# Organ trajectory and heterogeneity in AD

<p align="center">
  <img src="figures/figure-1.png" width="950">
</p>

<p align="center">
Overview of the study.
</p>


# Overview

This repository contains the analysis code and example files for a plasma proteomics study of Alzheimer’s disease (AD). The project uses plasma protein profiles to derive organ-resolved risk indices and characterize their trajectories and heterogeneity during the 15 years preceding AD.


## Repository Structure

```text
.
|-- data/example/       Example input files organized by analysis module
|-- docs/               Data dictionary and R dependency list
|-- figures/            Figures displayed in this README
|-- scripts/            R and Python analysis scripts
|-- .gitignore
|-- LICENSE
|-- README.md
`-- Requirements.txt
```

The directories under `scripts/` and `data/example/` use the same module names:

| Module | Analysis |
| --- | --- |
| `01-proteome` | PERMANOVA, residual adjustment, differential analysis, and Mfuzz clustering |
| `02-organ-risk` | Organ risk index, LASSO feature selection, and linear models |
| `03-subtyping` | NMF subtyping, differential proteins, enrichment, networks, and driver scores |
| `04-functional-analysis` | Organ-specific gene sets, ssGSEA, and pathway-level linear models |
| `05-survival-analysis` | Cox proportional-hazards models |

## Example Data

| Module | Files provided |
| --- | --- |
| `01-proteome` | `protein_expression.csv`, `protein_adjusted.csv`, `protein_time_series.csv`, `sample_metadata.csv` |
| `02-organ-risk` | `protein_adjusted.csv`, `protein_to_organ.tsv`, `organ_risk_index.csv`, `sample_metadata.csv` |
| `03-subtyping` | `organ_risk_index.csv`, `protein_adjusted.csv`, `cluster_membership.csv`, `protein_network.csv`, `protein_to_organ.tsv`, `sample_metadata.csv` |
| `04-functional-analysis` | `protein_adjusted.csv`, `protein_to_organ.tsv`, `organ_pathway_scores.csv`, `sample_metadata.csv` |
| `05-survival-analysis` | `organ_risk_index.csv`, `sample_metadata.csv` |

The example files are synthetic and intended to demonstrate input schemas. They are not real study data.

## Installation

Python 3.10+ dependencies are listed in `Requirements.txt`:

```bash
python -m venv .venv
source .venv/bin/activate       # Windows: .venv\Scripts\activate
python -m pip install -r Requirements.txt
```

R 4.4.1+ package requirements are listed in `docs/r-packages.txt`. 

```r
install.packages("BiocManager")
packages <- scan("docs/r-packages.txt", what = character(), comment.char = "#", quiet = TRUE)
BiocManager::install(packages, ask = FALSE, update = FALSE)
```

- **Tested environment:** [Ubuntu 22.04.5 LTS], [Intel Xeon(R) Platinum 8383C CPU and 512GB RAM], [Python 3.12.4], [R 4.4.1]. 
- **Typical installation time:** approximately 30 minutes for Python dependencies and 15 minutes for the R packages on the computer above, with Python and R already installed.

## Usage

Run scripts from the matching example directory.

```bash
# Proteome differential analysis
cd data/example/01-proteome
Rscript ../../../scripts/01-proteome/02-differential-analysis.R

# Organ risk index
cd ../02-organ-risk
python ../../../scripts/02-organ-risk/01-organ-risk-index.py
```

- **Differential analysis:** writes `data/example/01-proteome/diff/differential_results.csv` with one row per protein feature and columns for the effect estimate, standard error, t statistic, P value and adjusted FDR. Expected runtime is approximately 1 minute.
- **Organ risk index:** writes five CSV files under `data/example/02-organ-risk/results/` with mean risk scores, score standard deviations and an organ-level AUC summary. Expected runtime is approximately 20 minutes with the script's default settings.

## Input Requirement

- Sample IDs must match exactly across files used in the same analysis.
- Clinical groups are encoded as `NC` and `AD`.
- Organ mappings are tab-separated files with columns `Organ` and `Protein`.
- Human gene symbols are required for overlap with MSigDB in functional analysis.
- Survival metadata requires `Group` and `Year`; `E4_status` and `Range` are optional interaction variables.

See `docs/data-dictionary.md` for column definitions.

## License

The code is released under the MIT License.

## Contact

For questions, collaborations, or bug reports, please open an issue or contact:

Xuhua CHEN, [xhchena@ust.hk](mailto:xhchena@ust.hk); Mingfeng GUAN, [mf.guan@siat.ac.cn](mailto:mf.guan@siat.ac.cn); Yinhu LI, [xjy005351@siat.ac.cn](mailto:xjy005351@siat.ac.cn).
