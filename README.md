# HNSCC methylation analysis

Analysis for the poster **Association Between HPV Status and DNA Methylation
in Head and Neck Squamous Cell Carcinoma**.

[View the poster](poster/HNSCC_Poster.pdf)

## Project structure

```text
archive/                     Original mixed R/Python research notebook
poster/                      Final poster
R/setup.R                    Shared R dependencies
scripts/R/                   Ordered R analysis stages
scripts/python/              Python data-preparation scripts
scripts/run_r_analysis.R     R analysis runner
scripts/run_python_step.py   Python step runner
data/reference/              CpG reference lists
data/local/                  Local data workspace (ignored by Git)
results/figures/             Analysis figures
results/tables/              Gene-level result tables
results/report/              Rendered HTML analysis
```

## Data

The data are not included in this repository. See
[`data/README.md`](data/README.md) for the expected input filenames and local
setup.

## Software setup

Python preprocessing requires Python 3 and the packages in
`requirements.txt`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Install the R dependencies from CRAN and Bioconductor:

```r
install.packages(c(
  "arrow", "data.table", "dplyr", "ggplot2", "ggrepel", "gridExtra",
  "htmltools", "kableExtra", "knitr", "magick", "pheatmap", "readr",
  "stringr", "survival", "survminer", "tidyr", "webshot", "writexl",
  "BiocManager"
))
BiocManager::install(c("ENmix", "wateRmelon"))
```

## Running the analysis

Run a Python preparation step in the local data workspace:

```bash
python scripts/run_python_step.py HNSCC_Analysis.py --data-dir data/local
python scripts/run_python_step.py remove_first_row.py --data-dir data/local
```

Run all R stages:

```bash
Rscript scripts/run_r_analysis.R data/local
```

Generate the Kaplan-Meier survival figure directly:

```bash
Rscript scripts/plot_survival.R
```

To stop after a particular R stage, pass its filename after the data
directory. All preceding stages run first:

```bash
Rscript scripts/run_r_analysis.R data/local \
  05_global_methylation_subgroups.R
```

The Python preparation sequence is:

1. `HNSCC_Analysis.py`
2. `remove_first_row.py`
3. `calculate_methyl.py`
4. `adjust_calculation.py`
5. `sum_append_clinical.py`
6. `remove_nan_beta_value.py`
7. `split_subdata.py`
8. `calculate_15th_percentile.py`
9. `calculate_num_patients_add_column.py`
10. `methyl_remove_genes.py`
11. `extract_first_row_new.py`
12. Epigenetic-clock branch: `methyl_clean_for_clock.py`,
    `methyl_impute_for_clock.py`, `methyl_filter_clock.py`, and
    `methyl_filter_clock_new.py`
13. Gene branch: `filter_genes.py`, `methyl_clean_genes.py`, and
    `methyl_clean_gene_by_group.py`

`modify_table1.py` is an optional figure-processing step.
