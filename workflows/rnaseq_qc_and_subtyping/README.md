# PANTHER RNA-seq QC and SSPBC subtyping/ROR

Workflow author/maintainer: Emmanouil G. Sifakis (@sifakise).

## Purpose and scope

This module implements PANTHER STAR-based technical QC, expression-level QC,
SSPBC intrinsic subtyping, and SSPBC ROR calculation for the confirmed final
444-patient main-analysis set. It starts with 506 sequenced and quantified
RNA-seq samples;
it does not describe selection from all trial participants. It runs independently
of APM-specific analyses and does not source or require APM analysis code. It does
not reproduce all manuscript results or figures and is not generalized to other
cohorts. See [CITATION.md](CITATION.md) for the suggested module citation and
underlying method/package citations. See the local
[contribution record](CONTRIBUTORS.md) and [changelog](CHANGELOG.md).

For the corrected manuscript figure, see [Reproducing Supplementary Figure 13](docs/supplementary_figure_S13.md).

## Recorded environment

Use the repository-level [renv.lock](../../renv.lock),
[package versions](../../environment/package_versions.tsv),
[dependency provenance](../../environment/dependency_provenance.tsv), and
[session record](../../environment/sessionInfo.txt). Validation used R 4.6.0 and
the recorded installed package versions. All required packages must already be
available; these scripts stop on missing dependencies and do not install them.

Exact sspbc 1.0 and org.Hs.eg.db 3.23.1 source archives are now documented
and tested; their original historical archive identities remain unproven.
Reconstruction in an isolated Windows/R 4.6.0 library passed, retaining every
previously locked version and explicitly recording otel 0.2.0. The three stages
and both tests passed, with all 18 substantive outputs byte-identical to the
validated reference. See the [restoration instructions](../../environment/RESTORE.md)
for the resolved source/order issues and isolated execution setup. The copied
`.Rprofile` and historical session records remain unchanged. Restored-environment
validation used `Rscript --vanilla` with a private launcher selecting only the
fresh library and R's base/recommended sandbox; original user/site packages were
excluded. This validates the specified Windows/R environment and authorized
PANTHER inputs, not all operating systems or all APM analyses.

## Configuration and path resolution

Launch from the repository root, `PANTHER-APM-score`, using the
module-relative script paths below. Configuration is supplied through
the documented environment variables.

The scripts do not read a configuration file, parse YAML, or accept workflow
command-line arguments. [config/config.example.yml](config/config.example.yml)
is a placeholder planning template only; its location does not affect execution.
If a private launcher reads a local configuration, that launcher must set the
listed environment variables before invoking these scripts. A path to such a
configuration file is distinct from the data/output paths stored in it: the
workflow has no config-file-relative path resolution.

All relative data/output paths in environment variables are resolved against
the process working directory (the repository root), not the script directory,
module directory, or configuration-file directory. Prefer absolute authorized
local input paths and fresh private output directories outside this repository.
Do not place populated configurations or outputs in the repository.

| Stage | Required environment variable | Meaning |
| --- | --- | --- |
| 1 | `PANTHER_QCA_MULTIQC_DIR` | Directory containing exactly one `multiqc_star.txt` with `Sample`, `uniquely_mapped_percent`, and `total_reads` columns |
| 1 | `PANTHER_QCA_OUTPUT_DIR` | New or empty private technical-QC output directory |
| 1 | `PANTHER_QCA_SAMPLE_ID_REGEX` | Study sample pattern, `(XI-4066|XB-3945)-[0-9]{5}` |
| 2 | `PANTHER_TXIMPORT_RDS` | tximport RDS with `countsFromAbundance = "lengthScaledTPM"` |
| 2 | `PANTHER_STAGE_A_QC_FILE` | Stage 1 `multiqc_qc_calls_and_inclusion.csv` |
| 2 | `PANTHER_QC_OUTPUT_DIR` | Fresh private expression-QC output directory |
| 3 | `PANTHER_SSPBC_FPKM_RDS` | Nonnegative finite gene-by-sample FPKM matrix RDS, with source gene symbols as row names |
| 3 | `PANTHER_SSPBC_COUNTS_RDS` | Nonnegative finite count matrix RDS with identical gene rows and sample columns to the FPKM input |
| 3 | `PANTHER_MAIN_ANALYSIS_INCLUDE_FILE` | Stage 2 ordered `INCLUDE_main_analysis_FINAL.txt` |
| 3 | `PANTHER_SSPBC_OUTPUT_DIR` | Fresh private subtype/ROR output directory |

`PANTHER_QCA_EXPECTED_SAMPLE_COUNT` and `PANTHER_QCA_EXPECTED_PATIENT_COUNT`
are optional positive-integer input checks; use 506 for these validated inputs.
They check observed counts without selecting or padding samples.
`PANTHER_WRITE_LEGACY_QC_ALIASES` defaults to `false`; `1`, `true`, or `yes`
(case-insensitive) writes the existing legacy final inclusion-list aliases.
Keep it `false` to match the validated staging output file set.

## Execution order

After the isolated restored R library paths have been explicitly selected by
a private launcher and the required environment variables are set, run these
commands from the repository root. Each stage must complete successfully before
the next starts. The R executable must be the recorded R 4.6.0 runtime.

```powershell
Rscript --vanilla workflows/rnaseq_qc_and_subtyping/scripts/01_multiqc_technical_qc.R
Rscript --vanilla workflows/rnaseq_qc_and_subtyping/scripts/02_expression_qc_core.R
Rscript --vanilla workflows/rnaseq_qc_and_subtyping/scripts/03_intrinsic_subtyping_and_ror_sspbc.R
```

The stage 2 technical input must point to this run's stage 1 CSV. The stage 3
inclusion input must point to this run's stage 2 main-analysis list. Use fresh
output directories for all stages: stage 1 explicitly rejects a nonempty
output directory; stages 2 and 3 must also be given fresh locations to avoid
overwriting earlier outputs.

## Outputs and interpretation

| Stage | Principal outputs | Meaning |
| --- | --- | --- |
| 1 | `multiqc_qc_calls_and_inclusion.csv` | Ordered technical calls, reasons, and strict/main eligibility |
| 1 | `FAILED_samples.txt`, `FLAGGED_samples.txt` | Ordered technical FAIL and FLAG lists |
| 1 | `INCLUDE_primary_PASS_only.txt`, `INCLUDE_sensitivity_PASS_or_FLAG.txt` | Strict and main technical sets; historical filenames are preserved |
| 2 | `QC_decisions_per_sample.tsv` | Technical calls, expression-review metrics/flags, final inclusion, and reasons |
| 2 | `INCLUDE_strict_QC_FINAL.txt`, `INCLUDE_main_analysis_FINAL.txt` | Ordered strict and main final sets |
| 3 | `PANTHER_SSPBC_subtype_and_ROR_main_analysis.tsv` | Sample/patient IDs, four-class subtype call, subtype probability, and ROR |
| 3 | `PANTHER_SSPBC_gene_mapping_summary.tsv`, `PANTHER_SSPBC_model_gene_coverage.tsv` | Mapping and complete model-gene coverage checks |

Stages 1 and 2 also write `QC_sample_flow_summary.tsv`, `QC_thresholds.tsv`,
`QC_tool_versions.tsv`, and `QC_sessionInfo.txt` in their separate directories.
Stage 3 writes `PANTHER_SSPBC_tool_versions.tsv` and
`PANTHER_SSPBC_sessionInfo.txt`. All generated files remain private.

## Study-specific assumptions and safeguards

STAR `total_reads` represents post-trimming input read pairs/fragments for these
paired-end data, not individual mate-read counts. Technical thresholds are
FAIL below 10 million input fragments or 35% uniquely mapped, FLAG from 10 to
below 15 million or 35% to below 50%, otherwise PASS. A FAIL takes precedence.
These thresholds are specific to this study.

Expression QC uses tximport lengthScaledTPM-derived counts, not raw read counts.
Detected genes are counted at >= 10 on the full unrounded count matrix before
filtering. Gene filtering uses rounded counts and one explicit expression group;
PCA and Spearman sample correlation use at most 5,000 genes ranked by VST
variance. RLE uses all retained positive-variance VST genes. Robust covariance
uses seed 123 and a squared-distance chi-squared 99.5% review cutoff.

Technical calls and expression-review flags are distinct. Strict final inclusion
requires technical PASS and zero primary expression flags; main inclusion requires
technical PASS or FLAG and at most one sensitivity-level flag. The primary and
sensitivity tail thresholds differ, as recorded in each run's threshold table.
No claim of independence among QC metrics is made.

Technical output preserves normalized STAR row order. Expression QC aligns
technical calls by ID to the expression-column order. Subtyping uses the final
main-analysis inclusion-file order. Stage 3 requires exactly 444 unique samples
and 444 unique five-digit patient identifiers, with all included samples present
in both matrices. This safeguards the confirmed analysis set size; membership
comes from stage 2. It does not independently compare against a stored reference
membership list. Relocation validation compares that ordered list privately.

Subtype input is log2(FPKM + 1), passed through BreastSubtypeR mapping with
`RawCounts = FALSE`. SSPBC ROR calculation uses linear FPKM with the
`ssp.ror` model in the `sspbc` package. Supported alias-to-Entrez mappings
are checked against the recorded annotation package, and both models require
complete gene coverage in their actual inputs. In BreastSubtypeR 1.4.0,
`impute = FALSE` disables KNN replacement of NA matrix values; it does not create
absent model-gene rows. The explicit coverage checks remain necessary.

The validated staging reference produced 506 unique samples/patients,
365/79/62 technical PASS/FLAG/FAIL, 365/444 strict/main technical sets,
301/444 final strict/main sets, 444 subtype calls, and 444 finite ROR values.
These are comparison targets, never values to force. See
[method details](docs/qc_method_summary.md) and [sample flow](docs/sample_flow.md).
Validation on the
recorded inputs/environment is not a clinical diagnostic claim or a guarantee
for other cohorts or package versions.

## Existing validation tests

The two existing tests accept file-path arguments. From the repository root,
with the same environment and output variables used above:

```powershell
Rscript --vanilla workflows/rnaseq_qc_and_subtyping/tests/validate_sample_sets.R (Join-Path $env:PANTHER_QC_OUTPUT_DIR 'INCLUDE_strict_QC_FINAL.txt') (Join-Path $env:PANTHER_QC_OUTPUT_DIR 'INCLUDE_main_analysis_FINAL.txt')
Rscript --vanilla workflows/rnaseq_qc_and_subtyping/tests/validate_subtyping_outputs.R $env:PANTHER_MAIN_ANALYSIS_INCLUDE_FILE (Join-Path $env:PANTHER_SSPBC_OUTPUT_DIR 'PANTHER_SSPBC_subtype_and_ROR_main_analysis.tsv')
```

These check counts, uniqueness, subset membership, output completeness, and
subtype alignment. Exact comparison with the prior validated STAR outputs is
an additional private relocation check; the bundled tests alone do not establish
complete numerical equivalence.

## Controlled data

Inputs, pseudonymized sample/patient identifiers, individual QC decisions,
expression matrices, subtype probabilities, and ROR values are controlled
research data. Use approved storage and authorization. No real sample lists,
private paths/configurations, clinical data, or analytical results are included
here. Read the [controlled data notice](docs/controlled_data_notice.md) and
[input requirements](docs/input_requirements.md).

## Licensing

The original PANTHER QC/SSPBC module and its first-party supporting material are
available under the [MIT License](../../LICENSE).

Copyright (c) 2026 Emmanouil G. Sifakis

This notice applies to the original workflow contribution and its first-party
supporting material. It does not claim ownership of other contributors' analyses,
the underlying SSPBC models, external packages, or third-party bootstrap code.

Third-party material retains the notices documented in
[THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md), including the distributed
renv bootstrap. External dependencies retain their own licenses. The software
license does not provide permission to distribute controlled patient data; see
the [controlled data notice](docs/controlled_data_notice.md).

[CITATION.md](CITATION.md) and the module-specific [CITATION.cff.draft](CITATION.cff.draft)
retain draft module citation and release metadata. License selection and the
copyright notice are finalized for this contribution.
