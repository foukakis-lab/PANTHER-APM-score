# Input requirements

All input data remain controlled and must be supplied through authorized local paths.

## Technical QC

The MultiQC data directory must contain exactly one `multiqc_star.txt` table with the `Sample`, `uniquely_mapped_percent`, and `total_reads` columns. The `total_reads` field is interpreted as STAR input read pairs or fragments for these paired-end data.

Set the following environment variables:

- `PANTHER_QCA_MULTIQC_DIR`: MultiQC data directory.
- `PANTHER_QCA_OUTPUT_DIR`: empty or new technical-QC output directory.
- `PANTHER_QCA_SAMPLE_ID_REGEX`: project-specific sample identifier pattern.

Optional expected counts can be supplied through `PANTHER_QCA_EXPECTED_SAMPLE_COUNT` and `PANTHER_QCA_EXPECTED_PATIENT_COUNT`.

## Expression QC

The expression-QC script requires:

- a tximport RDS object with `countsFromAbundance = lengthScaledTPM`;
- the technical-QC decision table produced by script 01;
- an empty or new output directory.

Set `PANTHER_TXIMPORT_RDS`, `PANTHER_STAGE_A_QC_FILE`, and `PANTHER_QC_OUTPUT_DIR`.

## SSPBC subtype and ROR

The subtype script requires gene-level FPKM and count RDS objects with the same unique sample identifiers, plus the ordered final main-analysis inclusion file produced by script 02. Set `PANTHER_SSPBC_FPKM_RDS`, `PANTHER_SSPBC_COUNTS_RDS`, `PANTHER_MAIN_ANALYSIS_INCLUDE_FILE`, and `PANTHER_SSPBC_OUTPUT_DIR`.

Required R packages and validation versions are recorded in the repository-level
[environment directory](../../../environment/package_versions.tsv). Scripts stop when a required package is unavailable; they do not install software.

## Invocation and configuration

Launch from the repository root using the module-relative script paths in the
[module README](../README.md). Environment variables are the implemented
configuration mechanism; no YAML file is read and no workflow command-line
arguments are parsed. Relative data/output paths resolve against the repository
root working directory, regardless of a private configuration file location.

Use fresh private output directories outside the import tree. Optional
`PANTHER_WRITE_LEGACY_QC_ALIASES` defaults to `false`; keep that setting for
the validated output file set. The example YAML contains placeholders only.
