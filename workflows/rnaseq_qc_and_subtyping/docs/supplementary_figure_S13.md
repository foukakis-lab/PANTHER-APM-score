# Reproducing Supplementary Figure 13

[04_supplementary_figure_S13.R](../scripts/04_supplementary_figure_S13.R) is the authoritative entry point for the corrected manuscript Figure S13. It applies reproducible, image-specific corrections to the reviewed historical four-panel figure. The original image is required; the script cannot reconstruct the complete figure from expression data alone. It does not run the QC, subtyping, APM or survival pipelines or decide which patients enter an analysis.

## Required local inputs

Set each environment variable to an authorized local file. Keep all required inputs archived outside Git. Preserve the historical plotting scripts and exports unchanged, and keep local settings and generated figure files outside the checkout.

| Variable | Input and purpose |
|---|---|
| `PANTHER_S13_ORIGINAL` | The reviewed 3,920 x 3,288 PNG, or a local DOCX containing it at `word/media/image14.png`. The renderer verifies the embedded image SHA-256, not the document's filename. |
| `PANTHER_S13_CORRECTED_REFERENCE` | The retained corrected plotting-run `Figure_QC_summary_NatureStyle_FULL_4panel_reordered.png`. Used only to verify both axes against the newly keyed annotation vectors. The historical export with the same filename is not interchangeable. |
| `PANTHER_S13_PRESEQ` | Original MultiQC `mqc_preseq_plot_1.txt` with sample curves and total-read grid in millions. |
| `PANTHER_TXIMPORT_RDS` | Full authorized tximport list with finite gene-by-sample `counts`, unique sample column names and `countsFromAbundance = "lengthScaledTPM"`. |
| `PANTHER_S13_DECISIONS` | Canonical corrected `QC_decisions_per_sample.tsv`: `sample`, `technical_qc_status`, `detected_genes_counts_ge10`, `median_sample_correlation`, `pca_review_flag`. |
| `PANTHER_S13_HISTORICAL_MASTER` | Original plotting-run `QC_MASTER_samples_FULL.tsv`: `sample`, `detected_features_call`, `library_complexity`. Used to verify the old Detect stripe and the available preseq estimates, never to supply a missing estimate. |
| `PANTHER_S13_OUTPUT_DIR` | An existing, empty, private directory outside the repository. |
| `PANTHER_S13_PYTHON` | Python executable for the [raster/PDF renderer](../scripts/s13_render.py). Defaults to `python` on PATH. |

A and D have different denominators. A takes the **preseq estimated distinct reads, in millions, at exactly 442 million total reads**. A sample without a finite value at that coordinate is absent from A only. There is no interpolation, duplicate-fraction fallback or imputation. The violin and box plots use the verified technical QC classes.

D counts genes with **full, unrounded lengthScaledTPM counts >= 10** before filtering or VST. The R script checks the canonical per-sample counts against the tximport matrix. Detect marks values **strictly below the type-7 10th percentile over the full cohort**; it does not use the 0.5th-percentile strong inclusion flag. All joins use exact sample IDs.

## Execution

Use a local R installation with ggplot2, ragg, DESeq2, edgeR, S4Vectors, SummarizedExperiment and matrixStats, and a Python installation with NumPy, Pillow, ReportLab and pypdf. The script installs nothing. Use the separate figure-rendering environment documented below; the core workflow lockfile is unchanged and does not specify this rendering environment.

Set the paths above in your shell, then run:

```powershell
Rscript --vanilla workflows/rnaseq_qc_and_subtyping/scripts/04_supplementary_figure_S13.R
```

The output directory receives `S13_corrected.png`, its lossless single-page PDF `S13_corrected.pdf`, and the synchronized `S13_caption.md`. The console reports aggregate preseq availability and Detect changes. No ordered sample IDs are written. Existing output files are not replaced.

## Verification and scope

The historical canvas is hash-locked to `8b6656daef49eeb13cfb2d71923f2ba209f8fab384eb745be2da41ab7ed5183b`; the corrected reference PNG to `9d4d25ab39100292275fc3def47a95d1b7258736ff3d89e5c8e977fdd459cb9a`. Pixel coordinates are specific to these images. A differently cropped, scaled or recompressed image is not supported.

To recover named leaves for the existing anonymous heatmap, the script evaluates the historical blind VST, selects up to 5,000 positive-variance genes in decreasing variance order, computes Spearman sample correlations and clusters the correlation profiles by complete linkage with Euclidean distance. It checks the keyed median correlations against retained records. This computation supplies the leaf permutation only: the heatmap body and dendrograms are copied from the original. Both the old and corrected annotation vectors are checked at every cell center on both axes. Panel B is preserved, and panel C receives only a thick solid 99.9% legend key.

The verified cohort gave 494 preseq estimates and 12 unavailable. Of 506 samples, 51 had fewer than 14,660 detected genes; 42 calls differed from the historical track. These are reference results, not counts enforced by modifying data. The final renderer verifies that pixels outside the four prescribed edit rectangles are unchanged and that the PDF image stream equals the PNG RGB pixels.

Reproduction was checked on Windows with R 4.6.2, Arial, ggplot2 4.0.3, ragg 1.5.2, DESeq2 1.52.0, edgeR 4.10.7, matrixStats 1.5.0, Python 3.12.14, NumPy 2.3.5, Pillow 13.0.0, ReportLab 4.4.9 and pypdf 6.10.0. Font or graphics-library differences can change rendering bytes on other systems. No change to the repository lockfile or historical scripts is required.
