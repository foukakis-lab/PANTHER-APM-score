# PANTHER-APM-score

Analysis code for the PANTHER translational study of an antigen processing
and presentation machinery (APM) score and differential benefit from
adjuvant dose-dense chemotherapy in high-risk early breast cancer.

The repository separates RNA-seq quality control and molecular
characterization workflows from study-specific APM-score analyses.

## Available components

The current branch contains the
[RNA-seq QC and subtyping module](workflows/rnaseq_qc_and_subtyping/README.md).
It implements technical QC using STAR/MultiQC, expression-level QC,
SSPBC intrinsic subtyping, and RNA-seq risk-of-recurrence (ROR) calculation.

The module documentation describes execution order, required inputs,
configuration, outputs, study-specific safeguards, and validation.
The module can be run independently of APM-specific analyses.

APM-score development, validation, and downstream statistical analysis
code are not yet included in this branch. The available module does
not reproduce every manuscript result or figure.

## Software environment

The recorded R environment is stored at the repository root in
`renv.lock`, `.Rprofile`, and `renv/`.

This environment supports the QC/subtyping module. Exact-version
environment reconstruction and end-to-end execution were validated on
Windows x64 with R 4.6.0, renv 1.2.4, and Bioconductor 3.23.

See the [restoration instructions](environment/RESTORE.md) and
[package versions](environment/package_versions.tsv) for the tested
environment, installation procedure, and scope. This validation does
not establish compatibility with other operating systems or with
APM analysis code not yet included.

## Data access

Controlled patient data, expression matrices, private configurations,
and individual-level analytical outputs are not distributed in this
repository. Execution requires authorized inputs supplied through the
documented configuration mechanism.

See the
[controlled data notice](workflows/rnaseq_qc_and_subtyping/docs/controlled_data_notice.md).

## Contributions and citation

Confirmed contributions are documented in [CONTRIBUTORS.md](CONTRIBUTORS.md).
The QC/subtyping module includes its own
[author and maintainer information](workflows/rnaseq_qc_and_subtyping/README.md)
and [citation guidance](workflows/rnaseq_qc_and_subtyping/CITATION.md).

Repository-wide citation and release metadata remain under preparation
in [CITATION.cff.draft](CITATION.cff.draft).

## License

The original RNA-seq QC and SSPBC/ROR workflow contribution and its
first-party supporting material are available under the
[MIT License](LICENSE).

Copyright (c) 2026 Emmanouil G. Sifakis

This notice applies to that contribution, not to other contributors'
analyses or third-party software. Third-party material retains the
notices documented in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md),
and external dependencies retain their own licenses.

The software license does not authorize distribution of controlled
patient data.
