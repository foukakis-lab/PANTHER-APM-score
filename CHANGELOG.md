# Component changelogs

## APM analysis

- Added the main APM analysis (`workflows/apm_analysis/`), the RNA-seq-adapted pathologic SET2,3 calculation (`workflows/pathologic_set/`) and the whole-pipeline permutation test of the signature definition (`workflows/permutation_test/`).
- Added synthetic PANTHER and TEX example inputs in `data/` and the script that makes them (`workflows/apm_analysis/scripts/PANTHER_make_synthetic_data.R`, seed 2027).
- Added `renv.apm.lock` for the APM analysis environment (R 4.6.0, Bioconductor 3.23, GSVA 2.6.15).
- TEX probe sets are mapped to gene symbols with the GEO platform annotation GPL10379.

## RNA-seq QC and SSPBC module

The [QC/SSPBC module changelog](workflows/rnaseq_qc_and_subtyping/CHANGELOG.md) records the workflow's change history, scientific corrections, and validation.
