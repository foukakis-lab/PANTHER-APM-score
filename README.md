# PANTHER-APM-score

Analysis code for the PANTHER translational study of an antigen processing and presentation machinery (APM) score and differential benefit from adjuvant dose-dense chemotherapy in high-risk early breast cancer.

The APM80 gene set, its modules and submodules, and all gene memberships were predefined independently of clinical outcomes. 
Selection of the final signature components is a separate analysis using adjusted BCRFS treatment-by-score interaction tests in PANTHER.

## Workflows

The public release uses the following workflow structure:

| Component | Directory |
| --- | --- |
| Main APM analysis | `workflows/apm_analysis/` |
| SET2,3 calculation | `workflows/pathologic_set/` |
| APM permutation analysis | `workflows/permutation_test/` |
| RNA-seq QC and SSPBC intrinsic subtyping | `workflows/rnaseq_qc_and_subtyping/` |

R scripts belong in the `scripts/` subdirectory of each workflow. 
The main analysis uses the SET2,3 script and runs the permutation analysis by default (about 30 minutes; it can be switched off in its settings).

The existing [QC and subtyping workflow](workflows/rnaseq_qc_and_subtyping/README.md) contains technical QC, expression-level QC, and SSPBC intrinsic subtyping.
It runs independently of the APM analysis and requires its own authorized inputs.

## Data

See the [data download and example-input guide](data/README.md) for source
links, expected filenames, and local folder locations.

Public annotation, MATADOR, SCAN-B, and TEX expression data are obtained from their original sources rather than redistributed in this repository.
Real PANTHER participant data and the non-public TEX clinical and sample mapping files are not distributed here.

Synthetic PANTHER clinical and expression data, and synthetic TEX clinical data and sample mapping (made with `workflows/apm_analysis/scripts/PANTHER_make_synthetic_data.R`), are included in the data/ folder for the public example.
They demonstrate code execution and do not reproduce the study results. 

The processed PANTHER examples are intended for the downstream APM analysis; 
they do not provide the sequencing and other inputs needed to reproduce the upstream QC and subtyping workflow. 
See the [controlled data notice](workflows/rnaseq_qc_and_subtyping/docs/controlled_data_notice.md) for the existing QC module.

## Software and execution

Use the repository root as the working directory. The shared `data/` directory is located beside `workflows/`; paths should be resolved from the repository root.

Follow each component's documentation for its required inputs, dependencies, execution order, and output locations. 
   ## Running the APM analysis
   1. Open R in the repository root and install the packages: `install.packages("renv")`, then `renv::restore(lockfile = "renv.apm.lock")`.
   2. Download the public inputs (annotation files, SCAN-B and MATADOR) into `data/` as described in the [data guide](data/README.md); the TEX GEO files are downloaded automatically on the first run. 
      The synthetic PANTHER and TEX files already in `data/` let the analysis run end to end; to reproduce the study results, replace them with the original files under the same names.
   3. Run `workflows/apm_analysis/scripts/PANTHER_RNAseq_APM_score_Whole_Analysis.R` from top to bottom. It uses the SET2,3 functions in `workflows/pathologic_set/` and writes all results to `results/APM/`. Saving tables as PNG needs Chrome or Microsoft Edge (set chrome_path if neither is found).
   4. Part 1 runs the whole-pipeline permutation test (5,000 permutations, about 30 minutes; longer with fewer cores) and writes it to results/APM/APM_discovery_permutation/. To skip it, set SIGDEF_RUN_PERMUTATION = FALSE; Part 1 then uses a saved summary there if one exists.

   The first run needs internet access to download GEO GSE56493 and GPL10379 (TEX) and the MCP-counter gene list.

[renv.lock](renv.lock) records the QC and subtyping environment; [renv.apm.lock](renv.apm.lock) records the APM analysis environment (R 4.6.0, GSVA 2.6.15).

## Documentation and citation

Contributor information and citation guidance are provided within each component. See the [component contribution records](CONTRIBUTORS.md).

Study-level analysis documentation and publication citation metadata will be added on publication. Cite the original public data sources and underlying methods when using them.

## Licensing

See [LICENSE](LICENSE) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
for the terms applicable to the included code. Component documentation identifies the scope of the corresponding attribution and licensing.
External datasets retain their original terms of use; the code license does not grant access to controlled participant data.
