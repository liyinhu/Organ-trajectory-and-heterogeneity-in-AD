# Data dictionary

All files under `data/example/` are synthetic and contain no study participants.

## Common files

### `protein_expression.csv` and `protein_adjusted.csv`

The first column is `Sample_ID`. Remaining columns contain numeric protein abundances, with one protein per column.

### `protein_time_series.csv`

The first column is `Feature_ID`. Remaining columns are ordered time points or longitudinal intervals.

### `sample_metadata.csv`

Only columns required by the selected script need to be present.

| Column | Description | Expected format |
| --- | --- | --- |
| `Sample_ID` | Unique sample identifier | Character |
| `Group` | Clinical group | `NC` or `AD` |
| `Sex` | Sex covariate | Study-specific factor |
| `Age` | Age | Numeric |
| `BMI` | Body mass index | Numeric |
| `TDI` | Deprivation/index covariate | Numeric |
| `HQ` | Highest qualification/education | Study-specific factor |
| `E4_status` | APOE epsilon-4 status | `E4_carrier` or `E4_nocarrier` |
| `Range` | Time/stage stratum | Study-specific factor |
| `Year` | Event or censoring time for Cox analysis | Positive numeric |

### `protein_to_organ.tsv`

A tab-separated table with columns `Organ` and `Protein`.

### `organ_risk_index.csv`

The first column is `Sample_ID`. Remaining numeric columns are organ-level risk indices, with one organ per column.

### `cluster_membership.csv`

The first column is `Sample_ID`. `Cluster` contains the subtype assignment; `Range` is optional metadata.

### `protein_network.csv`

Required columns are `Protein1`, `Protein2`, `Partial_r`, `P_value`, and `Cluster`.

### `organ_pathway_scores.csv`

The first column is `Pathway`, using names such as `Brain__GO_NEURON_DEATH`. Remaining columns are samples and contain numeric ssGSEA scores.

## Module-specific requirements

- `01-proteome`: protein matrices and metadata must share sample IDs.
- `02-organ-risk`: mapping protein IDs must match expression-matrix columns.
- `03-subtyping`: cluster IDs must match the values used in downstream scripts.
- `04-functional-analysis`: use human gene symbols for MSigDB overlap.
- `05-survival-analysis`: metadata must contain `Group` and `Year`; organ-index and metadata files must share sample IDs.

Never commit direct identifiers, raw clinical records, or unapproved individual-level study data.
