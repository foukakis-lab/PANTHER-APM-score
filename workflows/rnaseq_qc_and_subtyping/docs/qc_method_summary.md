# QC method summary

## Technical QC

Technical QC uses the MultiQC STAR table. Sequencing yield is the STAR `total_reads` value, interpreted for these paired-end data as post-trimming input read pairs or fragments. Uniquely mapped percentage is read from `uniquely_mapped_percent`.

Input-fragment values below 10 million are hard failures, values from 10 million to below 15 million are review flags, and values at least 15 million pass. Uniquely mapped percentages below 35 percent are hard failures, values from 35 percent to below 50 percent are review flags, and values at least 50 percent pass. The main-analysis technical set includes PASS and FLAG samples.

## Expression QC

The expression input is a tximport object created with `countsFromAbundance = lengthScaledTPM`. The detected-gene metric counts genes with these abundance-derived counts greater than or equal to 10; it is not a raw-read-count threshold.

Genes are filtered with an explicit single-level group in `edgeR::filterByExpr`. Variance-stabilizing transformation uses an intercept-only DESeq2 model with blind transformation. Four sample-level review dimensions are calculated:

- robust PCA distance;
- relative log expression;
- median sample-sample correlation;
- detected-gene count.

The strict QC set requires technical QC PASS and zero expression review flags. The main-analysis set requires technical QC PASS or FLAG and at most one expression review flag. Exclusion reasons retain all applicable technical and expression-level reasons.

## SSPBC intrinsic subtyping

The four-class subtype model returns LumA, LumB, Her2, or Basal. Subtyping uses `BreastSubtypeR::Mapping` with `RawCounts = FALSE` followed by `BreastSubtypeR::BS_sspbc` with `ssp.name = "ssp.subtype"`. The input is log2(FPKM + 1), as supported by the validated package interface.

Gene identifiers are mapped to Entrez identifiers using the validated shared
alias dictionary. Complete `ssp.subtype` model-gene coverage is required, and
output order is matched exactly to the main-analysis inclusion file. The output
contains sample and patient identifiers, subtype labels, and subtype probabilities.
