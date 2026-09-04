"""Calculate organ risk indices from protein expression data.

The script balances AD and NC samples within each ``Range`` stratum with
SMOTE, fits a calibrated random forest for each organ, and writes mean and
standard-deviation predictions across bootstrap iterations.

Example
-------
python organ-risk-index.py \
    --protein protein_adjusted.csv \
    --organ protein_to_organ.tsv \
    --metadata sample_metadata.csv \
    --output results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate calibrated organ risk indices with SMOTE and random forests."
    )
    parser.add_argument("--protein", "-p", type=Path, default=Path("protein_adjusted.csv"))
    parser.add_argument("--organ", "-o", type=Path, default=Path("protein_to_organ.tsv"))
    parser.add_argument(
        "--metadata", "--info", "-i", dest="metadata", type=Path,
        default=Path("sample_metadata.csv")
    )
    parser.add_argument("--output", "--result", "-r", dest="output", type=Path, default=Path("results"))
    parser.add_argument("--target-ad", type=int, default=None)
    parser.add_argument("--target-nc", type=int, default=None)
    parser.add_argument("--n-bootstrap", type=int, default=100)
    parser.add_argument("--n-jobs", type=int, default=4)
    return parser.parse_args()


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    protein = pd.read_csv(args.protein, index_col=0)
    mapping = pd.read_csv(args.organ, sep="\t")
    metadata = pd.read_csv(args.metadata, index_col=0)

    required_mapping = {"Organ", "Protein"}
    required_metadata = {"Group", "Range"}
    if not required_mapping.issubset(mapping.columns):
        raise ValueError("The mapping file must contain 'Organ' and 'Protein' columns.")
    if not required_metadata.issubset(metadata.columns):
        raise ValueError("The metadata file must contain 'Group' and 'Range' columns.")
    if protein.index.has_duplicates or metadata.index.has_duplicates:
        raise ValueError("Sample IDs must be unique in both input tables.")
    if protein.columns.has_duplicates:
        raise ValueError("Protein IDs must be unique in the expression matrix.")

    common = protein.index.intersection(metadata.index)
    if common.empty:
        raise ValueError("No shared sample IDs were found between expression and metadata files.")
    protein = protein.loc[common].apply(pd.to_numeric, errors="coerce")
    if protein.isna().all(axis=1).any():
        bad = protein.index[protein.isna().all(axis=1)].tolist()
        raise ValueError(f"These samples contain no numeric protein values: {bad}")
    protein = protein.fillna(protein.median(numeric_only=True)).fillna(0.0).astype(float)
    metadata = metadata.loc[common].copy()
    metadata["Group"] = metadata["Group"].astype(str)
    metadata["Range"] = metadata["Range"].astype(str)
    if not set(metadata["Group"]).issubset({"AD", "NC"}):
        raise ValueError("Metadata Group values must be 'AD' or 'NC'.")
    return protein, mapping.astype({"Organ": str, "Protein": str}), metadata


def balance_by_range(
    protein: pd.DataFrame,
    metadata: pd.DataFrame,
    target_ad: int | None,
    target_nc: int | None,
) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    """Apply SMOTE within strata when both classes have enough neighbours."""
    protein_parts: list[pd.DataFrame] = []
    metadata_parts: list[pd.DataFrame] = []
    original_ids = set(protein.index.astype(str))

    for range_name, range_meta in metadata.groupby("Range", sort=False):
        range_protein = protein.loc[range_meta.index]
        labels = range_meta["Group"].map({"NC": 0, "AD": 1}).astype(int)
        counts = labels.value_counts().to_dict()
        n_ad, n_nc = counts.get(1, 0), counts.get(0, 0)
        default_target = max(n_ad, n_nc)
        desired_ad = max(n_ad, target_ad if target_ad is not None else default_target)
        desired_nc = max(n_nc, target_nc if target_nc is not None else default_target)

        can_smote = n_ad >= 2 and n_nc >= 2 and (desired_ad > n_ad or desired_nc > n_nc)
        if not can_smote:
            protein_parts.append(range_protein)
            metadata_parts.append(range_meta)
            continue

        sampling_strategy = {1: desired_ad, 0: desired_nc}
        k_neighbors = min(5, min(n_ad, n_nc) - 1)
        from imblearn.over_sampling import SMOTE

        sampler = SMOTE(
            sampling_strategy=sampling_strategy,
            k_neighbors=max(1, k_neighbors),
            random_state=42,
        )
        x_res, y_res = sampler.fit_resample(range_protein, labels)
        n_original = len(range_protein)
        synthetic_ids = [f"{range_name}_SMOTE_{i + 1}" for i in range(len(x_res) - n_original)]
        res_ids = list(range_protein.index.astype(str)) + synthetic_ids
        x_res_df = pd.DataFrame(x_res, columns=protein.columns, index=res_ids)
        y_res_series = pd.Series(y_res, index=res_ids).map({0: "NC", 1: "AD"})
        meta_res = pd.DataFrame({"Range": range_name, "Group": y_res_series}, index=res_ids)
        protein_parts.append(x_res_df)
        metadata_parts.append(meta_res)

    balanced_protein = pd.concat(protein_parts)
    balanced_metadata = pd.concat(metadata_parts)
    return balanced_protein, balanced_metadata, original_ids


def process_organ(
    organ: str,
    protein: pd.DataFrame,
    metadata: pd.DataFrame,
    organ_proteins: list[str],
    n_bootstrap: int,
) -> tuple[str, np.ndarray | None, np.ndarray | None, float | None, float | None]:
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.metrics import roc_auc_score

    features = [p for p in organ_proteins if p in protein.columns]
    if len(features) < 2:
        return organ, None, None, None, None

    x_values = protein[features].to_numpy(dtype=float)
    y_values = metadata.loc[protein.index, "Group"].map({"NC": 0, "AD": 1}).to_numpy()
    if len(np.unique(y_values)) < 2 or np.min(np.bincount(y_values)) < 2:
        return organ, None, None, None, None

    predictions: list[np.ndarray] = []
    aucs: list[float] = []
    for seed in range(n_bootstrap):
        try:
            model = RandomForestClassifier(
                n_estimators=50,
                max_depth=3,
                min_samples_split=2,
                random_state=seed,
            )
            model.fit(x_values, y_values)
            calibrated = CalibratedClassifierCV(model, method="sigmoid", cv=2)
            calibrated.fit(x_values, y_values)
            probabilities = calibrated.predict_proba(x_values)[:, 1]
            predictions.append(probabilities)
            aucs.append(roc_auc_score(y_values, probabilities))
        except (ValueError, RuntimeError):
            continue

    if not predictions:
        return organ, None, None, None, None
    prediction_array = np.vstack(predictions)
    return organ, prediction_array.mean(axis=0), prediction_array.std(axis=0), float(np.mean(aucs)), float(np.std(aucs))


def main() -> None:
    args = parse_args()
    from joblib import Parallel, delayed

    if args.n_bootstrap < 1 or args.n_jobs == 0:
        raise ValueError("n-bootstrap must be positive and n-jobs cannot be zero.")
    protein, mapping, metadata = load_inputs(args)
    balanced_protein, balanced_metadata, original_ids = balance_by_range(
        protein, metadata, args.target_ad, args.target_nc
    )
    args.output.mkdir(parents=True, exist_ok=True)

    proteins_by_organ = mapping.groupby("Organ")["Protein"].apply(list).to_dict()
    results = Parallel(n_jobs=args.n_jobs)(
        delayed(process_organ)(organ, balanced_protein, balanced_metadata, proteins, args.n_bootstrap)
        for organ, proteins in proteins_by_organ.items()
    )

    mean_values: dict[str, np.ndarray] = {}
    sd_values: dict[str, np.ndarray] = {}
    auc_rows: list[dict[str, float | str]] = []
    for organ, mean_pred, sd_pred, auc_mean, auc_sd in results:
        if mean_pred is None:
            continue
        mean_values[organ] = mean_pred
        sd_values[organ] = sd_pred
        auc_rows.append({"Organ": organ, "AUC_mean": auc_mean, "AUC_SD": auc_sd})

    mean_all = pd.DataFrame(mean_values, index=balanced_protein.index)
    sd_all = pd.DataFrame(sd_values, index=balanced_protein.index)
    original_index = [sample for sample in balanced_protein.index if str(sample) in original_ids]
    mean_original = mean_all.loc[original_index]
    sd_original = sd_all.loc[original_index]

    mean_original.to_csv(args.output / "organ_risk_index_mean.csv")
    sd_original.to_csv(args.output / "organ_risk_index_sd.csv")
    mean_all.to_csv(args.output / "organ_risk_index_mean_smote.csv")
    sd_all.to_csv(args.output / "organ_risk_index_sd_smote.csv")
    pd.DataFrame(auc_rows).to_csv(args.output / "organ_risk_index_auc.csv", index=False)
    print(f"Completed. Results written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
