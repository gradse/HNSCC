# Data

The analysis data are not included in this repository.

Place working copies in `data/local/`, or pass another working directory to
the R and Python runners. Several preparation steps create or replace
intermediate files, so use a working copy rather than the source data folder.

The workflow begins with:

```text
clinical_All_CDEs.txt
methylation_level3_within_bioassay.txt
```

`HNSCC_Analysis.py` converts the text methylation input to
`methylation_level3_within_bioassay.parquet`. The remaining filenames and
execution order are documented in the project README.

`data/reference/` contains the CpG reference lists used by the three
epigenetic clocks.
