# Suggested module citation

Draft suggested text; fill the release and commit placeholders only after they
are confirmed:

> Sifakis, Emmanouil G. PANTHER RNA-seq QC and SSPBC intrinsic subtyping/ROR
> workflow. Module `workflows/rnaseq_qc_and_subtyping` in
> `foukakis-lab/PANTHER-APM-score`. Release: [FUTURE_RELEASE].
> Commit: [FUTURE_COMMIT].

This module is available in `foukakis-lab/PANTHER-APM-score`.
Release, publication-date, DOI, and release-specific commit metadata
have not yet been finalized. Repository-wide citation metadata remain
in [CITATION.cff.draft](../../CITATION.cff.draft). Cite the associated
manuscript when its approved publication metadata are available.

## Underlying methods and packages

The workflow implementation is distinct from development of its
underlying models and packages. The citations below identify the
methods and software used by this module.

- Yang Q, Hartman J, Sifakis E (2025). BreastSubtypeR: a unified R/Bioconductor package for intrinsic molecular subtyping in breast cancer research. NAR Genomics and Bioinformatics 7(4), lqaf131. [doi:10.1093/nargab/lqaf131](https://doi.org/10.1093/nargab/lqaf131).
- Yang Q, Sifakis E (2026). BreastSubtypeR: Cohort-aware methods for intrinsic molecular subtyping of breast cancer. R package version 1.4.0. [Package DOI](https://doi.org/10.18129/B9.bioc.BreastSubtypeR).
- Staaf J, Hakkinen J, Hegardt C, et al. (2022). RNA sequencing-based
  single sample predictors of molecular subtype and risk of recurrence
  for clinical assessment of early-stage breast cancer. npj Breast
  Cancer 8, 94.
  [doi:10.1038/s41523-022-00465-3](https://doi.org/10.1038/s41523-022-00465-3).
- Vallon-Christersson J (2022). sspbc: Single-Sample Predictors for Breast Cancer (sspbc). R package version 1.0.
- Soneson C, Love MI, Robinson MD (2015). Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. F1000Research 4, 1521. [doi:10.12688/f1000research.7563.1](https://doi.org/10.12688/f1000research.7563.1).
- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology 15, 550. [doi:10.1186/s13059-014-0550-8](https://doi.org/10.1186/s13059-014-0550-8).
- Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK (2025). edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. Nucleic Acids Research 53(2), gkaf018. [doi:10.1093/nar/gkaf018](https://doi.org/10.1093/nar/gkaf018).
- Todorov V, Filzmoser P (2009). An Object-Oriented Framework for Robust Multivariate Analysis. Journal of Statistical Software 32(3), 1-47. [doi:10.18637/jss.v032.i03](https://doi.org/10.18637/jss.v032.i03).
- Pages H, Carlson M, Falcon S, Li N (2026). AnnotationDbi: Manipulation of SQLite-based annotations in Bioconductor. R package version 1.74.0. [Package record](https://bioconductor.org/packages/AnnotationDbi).
- Carlson M (2026). org.Hs.eg.db: Genome wide annotation for Human. R package version 3.23.1.

Reconstruction sources for sspbc 1.0 and org.Hs.eg.db 3.23.1 are
documented and have been tested. Identity of the original historical
source archives remains unproven. See the
[restoration instructions](../../environment/RESTORE.md).

AnnotationDbi and org.Hs.eg.db citation text is generated from installed package
metadata and should be reviewed when finalizing publication citations. Recorded
versions and provenance for the full environment remain in
[environment/](../../environment/package_versions.tsv).
