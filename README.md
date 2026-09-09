# PANTHER-APM-score

Planned repository: `foukakis-lab/PANTHER-APM-score`.

This local import candidate contributes the validated PANTHER RNA-seq technical
QC, expression QC, and SSPBC intrinsic subtyping/ROR workflow. It currently
contains QC and SSPBC/ROR code, not the complete APM analysis code. These three
scripts do not reproduce every manuscript result or figure.

Start with the [RNA-seq QC and subtyping module](workflows/rnaseq_qc_and_subtyping/README.md)
for execution order, inputs, environment variables, output descriptions, and
study-specific safeguards. The module runs independently of APM-specific analyses.

The recorded R environment remains at the repository root: `renv.lock`,
`.Rprofile`, `renv/`, and [environment records](environment/package_versions.tsv).
Exact-version reconstruction from documented sources and one complete QC/SSPBC
run passed on Windows x64 with R 4.6.0, renv 1.2.4, and Bioconductor 3.23.
All 18 substantive outputs matched the validated reference byte-for-byte, with
the original working libraries excluded. All previously locked versions remain;
otel 0.2.0 is now explicit. See [restoration instructions](environment/RESTORE.md)
for isolation, source records, resolved installation-order issues, and the tested
scope. Historical session records remain unchanged.

Only code and approved documentation belong in this candidate. Controlled data,
private configuration, and all validation and analytical outputs remain outside
this tree. See the [controlled data notice](workflows/rnaseq_qc_and_subtyping/docs/controlled_data_notice.md).

[Contributions](CONTRIBUTORS.md) are limited to confirmed workflow contributions.
The [module citation](workflows/rnaseq_qc_and_subtyping/CITATION.md) is suggested
text with release/commit placeholders. [CITATION.cff.draft](CITATION.cff.draft)
remains a draft. This copy has not been imported into the planned repository.

## Licensing

The original PANTHER QC/SSPBC workflow contribution and its first-party supporting
material are available under the [MIT License](LICENSE).

Copyright (c) 2026 Emmanouil G. Sifakis

This notice applies to that contribution. It does not claim ownership of other
contributors' analyses, the underlying SSPBC models, external packages, or
third-party bootstrap code.

Third-party material retains the notices documented in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). External dependencies retain
their own licenses. The software license does not provide permission to
distribute controlled patient data; see the
[controlled data notice](workflows/rnaseq_qc_and_subtyping/docs/controlled_data_notice.md).

[CITATION.cff.draft](CITATION.cff.draft) remains a draft only because repository-wide
contributor and release metadata are unconfirmed. License selection and the
copyright notice are finalized for this contribution.
