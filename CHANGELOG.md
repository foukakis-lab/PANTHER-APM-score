# Changelog

## Unreleased license finalization

- Finalized the standard MIT LICENSE for the original QC/SSPBC contribution and first-party supporting material with Copyright (c) 2026 Emmanouil G. Sifakis.
- Added the approved copyright and MIT SPDX comment headers to the three module scripts and two first-party tests; executable expressions are unchanged.
- Archived LICENSE.draft and LICENSE_SELECTION_REQUIRED.md outside the candidate. Clarified licensing scope, external dependency terms, preserved third-party notices, and the controlled-data boundary.
- Retained CITATION.cff.draft only for unconfirmed repository-wide contributor and release metadata; license selection and the copyright notice are finalized.
- Preserved the promoted environment metadata, bootstrap files, analytical outputs, and prior scientific-correction and validation history.

## Earlier licensing preparation (superseded)

The following entries record preparation before the approved notice was supplied.

- Recorded MIT selection for original PANTHER analysis code and first-party documentation; prepared LICENSE.draft pending confirmed copyright notice wording.
- Preserved external dependency terms and documented the distributed renv bootstrap notice.
- Retained proposed first-party SPDX headers privately until copyright confirmation; no executable analysis or environment changes.

## Unreleased module organization

- Relocated the validated STAR QC and SSPBC/ROR workflow into `workflows/rnaseq_qc_and_subtyping/` for planned import into `foukakis-lab/PANTHER-APM-score`.
- Clarified script comments, environment-variable configuration, root-based invocation, contribution scope, and pending citation/license metadata.
- Preserved repository-level environment and renv files byte-for-byte.
- Comment edits preserve identical parsed R expressions; relocation requires no executable path changes. Full relocation execution status is recorded in the private validation report.

## Unreleased publication code curation

- Replaced the biotype-category count with the validated detected-gene metric based on tximport `lengthScaledTPM`-derived counts greater than or equal to 10.
- Updated the strict final set from 300 to 301 patients after the detected-gene correction.
- Confirmed no change to the ordered 444-patient main-analysis set.
- Selected STAR input read pairs or fragments as the sequencing-yield metric.
- Confirmed exact agreement of technical classifications and final memberships with the historical samtools scenario.
- Restricted SSPBC subtype and ROR analysis to the ordered final 444-patient set.
- Removed clinical, survival, figure, exploratory, package-installation, and private-path code from the three core scripts.
