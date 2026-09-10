# Copyright (c) 2026 Emmanouil G. Sifakis
# SPDX-License-Identifier: MIT

# PANTHER expression QC and strict/main-analysis sample inclusion.
# Workflow author/maintainer: Emmanouil G. Sifakis.
# Inputs: lengthScaledTPM tximport RDS and stage 1 technical decision CSV.
# Outputs: expression-review decisions, ordered final sample sets, and thresholds.
# Execution: stage 2, after technical QC; run from the repository root.

options(warn = 1)

required_packages <- c(
  "DESeq2", "edgeR", "matrixStats", "readr", "rrcov",
  "S4Vectors", "SummarizedExperiment", "tximport"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

required_environment <- c(
  tximport_rds = "PANTHER_TXIMPORT_RDS",
  stage_a_qc_file = "PANTHER_STAGE_A_QC_FILE",
  output_directory = "PANTHER_QC_OUTPUT_DIR"
)
configuration <- vapply(
  required_environment,
  function(variable) Sys.getenv(variable, unset = ""),
  character(1)
)
if (any(!nzchar(trimws(configuration)))) {
  stop(
    "Required environment variables are unset: ",
    paste(
      unname(required_environment[!nzchar(trimws(configuration))]),
      collapse = ", "
    )
  )
}

tximport_path <- configuration[["tximport_rds"]]
stage_a_path <- configuration[["stage_a_qc_file"]]
output_dir <- configuration[["output_directory"]]
missing_inputs <- c(tximport_path, stage_a_path)[
  !file.exists(c(tximport_path, stage_a_path))
]
if (length(missing_inputs)) {
  stop("Missing required input files: ", paste(missing_inputs, collapse = "; "))
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Unable to create the output directory.")

normalize_sample <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("\\.bam$", "", x, ignore.case = TRUE)
  x <- sub("_Aligned(\\.sorted)?\\.out$", "", x)
  x <- sub("_aligned(\\.sorted)?\\.out$", "", x)
  x <- sub("_aligned\\.sortedbycoord\\.out$", "", x)
  x <- sub("\\.sorted$", "", x)
  x <- sub("\\.dedup$", "", x)
  x <- sub("_[Rr][12]$", "", x)
  x <- sub("_001$", "", x)
  x <- sub("\\.f(ast)?q(\\.gz)?$", "", x, ignore.case = TRUE)
  x <- sub("\\.lc_extrap$", "", x)
  sub("^.*?((XI-4066|XB-3945)-[0-9]{5})$", "\\1", x)
}

patient_id_from_sample <- function(x) {
  result <- sub("^.*-([0-9]{5})$", "\\1", x)
  result[!grepl("[0-9]{5}$", x)] <- NA_character_
  result
}

assert_identifiers <- function(x, label) {
  if (anyNA(x) || any(trimws(as.character(x)) == "")) {
    stop(label, " contains a missing identifier.")
  }
  if (anyDuplicated(x)) stop(label, " contains a duplicated identifier.")
}

collapse_reasons <- function(values) {
  values <- values[nzchar(values)]
  if (!length(values)) NA_character_ else paste(values, collapse = "; ")
}

txi <- readRDS(tximport_path)
if (!is.list(txi) || is.null(txi$counts)) stop("Invalid tximport object.")
if (!identical(as.character(txi$countsFromAbundance), "lengthScaledTPM")) {
  stop("The tximport object must record countsFromAbundance = lengthScaledTPM.")
}

cts <- as.matrix(txi$counts)
if (!is.numeric(cts) || any(!is.finite(cts)) || any(cts < 0)) {
  stop("The expression matrix must contain finite nonnegative values.")
}
colnames(cts) <- normalize_sample(colnames(cts))
assert_identifiers(colnames(cts), "Expression sample identifiers")
assert_identifiers(rownames(cts), "Expression gene identifiers")
sample_order <- colnames(cts)
patient_id <- patient_id_from_sample(sample_order)
assert_identifiers(patient_id, "Patient identifiers")

stage_a <- readr::read_csv(stage_a_path, show_col_types = FALSE, progress = FALSE)
required_stage_a_columns <- c("sample", "qc_status")
if (!all(required_stage_a_columns %in% names(stage_a))) {
  stop("Stage A decisions must contain sample and qc_status columns.")
}
stage_a$sample <- normalize_sample(stage_a$sample)
assert_identifiers(stage_a$sample, "Stage A sample identifiers")
if (!setequal(stage_a$sample, sample_order)) {
  stop("Stage A and expression sample sets differ.")
}
# Align technical calls to expression-column order by ID, never by row position.
stage_a <- stage_a[match(sample_order, stage_a$sample), , drop = FALSE]
if (!identical(stage_a$sample, sample_order)) {
  stop("Stage A sample order could not be aligned exactly.")
}
qc_status <- toupper(as.character(stage_a$qc_status))
if (anyNA(qc_status) || any(!qc_status %in% c("PASS", "FLAG", "FAIL"))) {
  stop("Stage A status values must be PASS, FLAG, or FAIL.")
}

# Number of genes with tximport length-scaled abundance-derived counts >= 10.
# This is not a raw read-count threshold. Calculate it on the full unrounded
# tximport count matrix, before gene filtering or variance transformation.
detected_genes_counts_gt0 <- colSums(cts > 0, na.rm = TRUE)
detected_genes_counts_ge10 <- colSums(cts >= 10, na.rm = TRUE)
if (
  any(!is.finite(detected_genes_counts_ge10)) ||
    any(detected_genes_counts_ge10 < 0) ||
    any(detected_genes_counts_ge10 > nrow(cts)) ||
    any(detected_genes_counts_ge10 != round(detected_genes_counts_ge10))
) {
  stop("Invalid detected-gene values.")
}

# Low detection is reviewed below the cohort 2% (strict) or 0.5% (main) quantile.
detected_primary_threshold <- unname(
  stats::quantile(detected_genes_counts_ge10, 0.02, na.rm = TRUE)
)
detected_sensitivity_threshold <- unname(
  stats::quantile(detected_genes_counts_ge10, 0.005, na.rm = TRUE)
)
detected_primary_flag <-
  detected_genes_counts_ge10 < detected_primary_threshold
detected_sensitivity_flag <-
  detected_genes_counts_ge10 < detected_sensitivity_threshold

# Retain genes using one explicit expression group across all samples.
cts_rounded <- round(cts)
expression_group <- factor(rep("all_samples", ncol(cts_rounded)))
keep_gene <- edgeR::filterByExpr(
  cts_rounded,
  group = expression_group,
  min.count = 10
)
if (!any(keep_gene)) stop("No genes passed expression filtering.")
retained_genes <- rownames(cts_rounded)[keep_gene]

# Preserve the tximport structure for DESeq2 construction.
txi_filtered <- txi
txi_filtered$counts <- cts
for (component in intersect(c("counts", "abundance", "length"), names(txi_filtered))) {
  component_matrix <- as.matrix(txi_filtered[[component]])
  colnames(component_matrix) <- normalize_sample(colnames(component_matrix))
  if (!setequal(colnames(component_matrix), sample_order)) {
    stop("A tximport matrix contains a different sample set: ", component)
  }
  if (!setequal(rownames(component_matrix), rownames(cts))) {
    stop("A tximport matrix contains a different gene set: ", component)
  }
  txi_filtered[[component]] <- component_matrix[
    retained_genes,
    sample_order,
    drop = FALSE
  ]
}

col_data <- S4Vectors::DataFrame(row.names = sample_order)
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi_filtered,
  colData = col_data,
  design = ~ 1
)
vst_object <- DESeq2::vst(dds, blind = TRUE)
vst_matrix <- SummarizedExperiment::assay(vst_object)
if (!identical(colnames(vst_matrix), sample_order)) {
  stop("VST sample order differs from the expression input order.")
}
if (any(!is.finite(vst_matrix))) stop("The VST matrix contains nonfinite values.")

positive_variance <- matrixStats::rowVars(vst_matrix) > 0
vst_matrix <- vst_matrix[positive_variance, , drop = FALSE]
if (!nrow(vst_matrix)) stop("No positive-variance genes remain after VST.")

# PCA and correlation use at most 5,000 genes with the highest VST variance,
# after expression filtering and removal of zero-variance genes.
top_index <- head(
  order(matrixStats::rowVars(vst_matrix), decreasing = TRUE),
  min(5000L, nrow(vst_matrix))
)
review_matrix <- vst_matrix[top_index, , drop = FALSE]

pca <- stats::prcomp(
  t(review_matrix),
  center = TRUE,
  scale. = FALSE,
  rank. = 10
)
component_count <- min(
  10L,
  ncol(pca$x),
  floor((nrow(pca$x) - 1L) / 2L)
)
pca_matrix <- pca$x[, seq_len(component_count), drop = FALSE]
if (!identical(rownames(pca_matrix), sample_order)) {
  stop("PCA sample order differs from the expression input order.")
}
if (nrow(pca_matrix) <= 2L * ncol(pca_matrix)) {
  stop("Too few samples for robust covariance estimation.")
}

# Fix the MCD subsampling seed. Review squared robust distances above the
# chi-squared 99.5% cutoff with one degree of freedom per retained component.
set.seed(123L)
robust_covariance <- rrcov::CovMcd(pca_matrix)
pca_distance_squared <- stats::mahalanobis(
  pca_matrix,
  center = robust_covariance@center,
  cov = robust_covariance@cov
)
pca_alpha <- 0.995
pca_threshold <- stats::qchisq(pca_alpha, df = component_count)
pca_review_flag <- pca_distance_squared > pca_threshold

# Calculate RLE summaries using every retained positive-variance gene.
# RLE uses the full retained VST matrix, not the PCA/correlation subset.
rle_values <- sweep(
  vst_matrix,
  1,
  matrixStats::rowMedians(vst_matrix),
  FUN = "-"
)
rle_median_absolute <- vapply(
  seq_len(ncol(rle_values)),
  function(index) stats::median(abs(rle_values[, index]), na.rm = TRUE),
  numeric(1)
)
names(rle_median_absolute) <- sample_order
rle_primary_threshold <- unname(
  stats::quantile(rle_median_absolute, 0.98, na.rm = TRUE)
)
rle_sensitivity_threshold <- unname(
  stats::quantile(rle_median_absolute, 0.995, na.rm = TRUE)
)
rle_primary_flag <- rle_median_absolute > rle_primary_threshold
rle_sensitivity_flag <- rle_median_absolute > rle_sensitivity_threshold

# Use the same finite variance-ranked VST matrix for PCA and sample correlation.
sample_correlation <- stats::cor(
  review_matrix,
  use = "everything",
  method = "spearman"
)
sample_correlation <- sample_correlation[
  match(sample_order, rownames(sample_correlation)),
  match(sample_order, colnames(sample_correlation)),
  drop = FALSE
]
if (
  !identical(rownames(sample_correlation), sample_order) ||
    !identical(colnames(sample_correlation), sample_order)
) {
  stop("Sample-correlation order differs from the expression input order.")
}
median_sample_correlation <- vapply(
  sample_order,
  function(sample_name) {
    values <- sample_correlation[, sample_name]
    stats::median(values[names(values) != sample_name], na.rm = TRUE)
  },
  numeric(1)
)
correlation_primary_threshold <- unname(
  stats::quantile(median_sample_correlation, 0.02, na.rm = TRUE)
)
correlation_sensitivity_threshold <- unname(
  stats::quantile(median_sample_correlation, 0.005, na.rm = TRUE)
)
correlation_primary_flag <-
  median_sample_correlation < correlation_primary_threshold
correlation_sensitivity_flag <-
  median_sample_correlation < correlation_sensitivity_threshold

primary_expression_flag_count <-
  as.integer(pca_review_flag) +
  as.integer(rle_primary_flag) +
  as.integer(correlation_primary_flag) +
  as.integer(detected_primary_flag)
sensitivity_expression_flag_count <-
  as.integer(pca_review_flag) +
  as.integer(rle_sensitivity_flag) +
  as.integer(correlation_sensitivity_flag) +
  as.integer(detected_sensitivity_flag)

# Expression-review flags remain distinct from technical PASS/FLAG/FAIL calls.
# Strict: technical PASS and no primary flags. Main: technical PASS or FLAG
# and at most one sensitivity-level flag; the two sets use distinct tail cutoffs.
strict_qc_include <-
  qc_status == "PASS" & primary_expression_flag_count == 0L
main_analysis_include <-
  qc_status %in% c("PASS", "FLAG") &
  sensitivity_expression_flag_count <= 1L

strict_reasons <- vapply(seq_along(sample_order), function(index) {
  collapse_reasons(c(
    if (qc_status[[index]] != "PASS") {
      paste0("Technical QC status: ", qc_status[[index]])
    } else "",
    if (pca_review_flag[[index]]) "PCA robust distance review flag" else "",
    if (rle_primary_flag[[index]]) "High median absolute RLE" else "",
    if (correlation_primary_flag[[index]]) {
      "Low median sample-sample correlation"
    } else "",
    if (detected_primary_flag[[index]]) "Low detected-gene count" else ""
  ))
}, character(1))

main_reasons <- vapply(seq_along(sample_order), function(index) {
  collapse_reasons(c(
    if (qc_status[[index]] == "FAIL") "Technical QC FAIL" else "",
    if (!qc_status[[index]] %in% c("PASS", "FLAG", "FAIL")) {
      "Technical QC status is not eligible"
    } else "",
    if (sensitivity_expression_flag_count[[index]] > 1L) {
      "More than one sensitivity-level expression review flag"
    } else "",
    if (pca_review_flag[[index]]) "PCA robust distance review flag" else "",
    if (rle_sensitivity_flag[[index]]) "Extreme median absolute RLE" else "",
    if (correlation_sensitivity_flag[[index]]) {
      "Very low median sample-sample correlation"
    } else "",
    if (detected_sensitivity_flag[[index]]) {
      "Very low detected-gene count"
    } else ""
  ))
}, character(1))

decisions <- data.frame(
  sample = sample_order,
  patient_id = patient_id,
  technical_qc_status = qc_status,
  detected_genes_counts_gt0 = as.numeric(detected_genes_counts_gt0),
  detected_genes_counts_ge10 = as.numeric(detected_genes_counts_ge10),
  pca_robust_distance_squared = unname(pca_distance_squared),
  pca_review_flag = unname(pca_review_flag),
  rle_median_absolute = unname(rle_median_absolute),
  rle_primary_flag = unname(rle_primary_flag),
  rle_sensitivity_flag = unname(rle_sensitivity_flag),
  median_sample_correlation = unname(median_sample_correlation),
  correlation_primary_flag = unname(correlation_primary_flag),
  correlation_sensitivity_flag = unname(correlation_sensitivity_flag),
  detected_genes_primary_flag = unname(detected_primary_flag),
  detected_genes_sensitivity_flag = unname(detected_sensitivity_flag),
  primary_expression_flag_count = primary_expression_flag_count,
  sensitivity_expression_flag_count = sensitivity_expression_flag_count,
  strict_qc_include = strict_qc_include,
  main_analysis_include = main_analysis_include,
  strict_exclusion_reasons = ifelse(strict_qc_include, NA_character_, strict_reasons),
  main_analysis_exclusion_reasons = ifelse(
    main_analysis_include,
    NA_character_,
    main_reasons
  ),
  stringsAsFactors = FALSE
)
if (!identical(decisions$sample, sample_order)) stop("Decision output order changed.")

count_patients <- function(include) {
  ids <- patient_id[include]
  c(
    samples = sum(include),
    unique_nonmissing_patients = length(unique(ids[!is.na(ids)])),
    missing_patient_ids = sum(is.na(ids)),
    duplicated_patient_ids = sum(duplicated(ids) | duplicated(ids, fromLast = TRUE))
  )
}
strict_counts <- count_patients(strict_qc_include)
main_counts <- count_patients(main_analysis_include)

flow_summary <- data.frame(
  metric = c(
    "input_samples",
    "unique_nonmissing_patient_ids",
    "missing_patient_ids",
    "duplicated_patient_ids",
    "technical_qc_PASS",
    "technical_qc_FLAG",
    "technical_qc_FAIL",
    "strict_qc_included_samples",
    "strict_qc_included_patients",
    "main_analysis_included_samples",
    "main_analysis_included_patients"
  ),
  value = c(
    length(sample_order),
    length(unique(patient_id[!is.na(patient_id)])),
    sum(is.na(patient_id)),
    sum(duplicated(patient_id) | duplicated(patient_id, fromLast = TRUE)),
    sum(qc_status == "PASS"),
    sum(qc_status == "FLAG"),
    sum(qc_status == "FAIL"),
    strict_counts[["samples"]],
    strict_counts[["unique_nonmissing_patients"]],
    main_counts[["samples"]],
    main_counts[["unique_nonmissing_patients"]]
  ),
  definition = c(
    "Expression samples entering QC",
    "Unique nonmissing patient identifiers",
    "Missing patient identifiers",
    "Patient identifiers represented by more than one sample",
    "Stage A PASS",
    "Stage A FLAG",
    "Stage A FAIL",
    "Technical QC PASS and no expression-level review flags",
    "Strict set; corresponds to the corrected former primary set",
    paste0(
      "Technical QC PASS or FLAG and no more than one ",
      "expression-level review flag"
    ),
    "Main-analysis set; corresponds to the former sensitivity set"
  ),
  stringsAsFactors = FALSE
)

thresholds <- data.frame(
  metric = c(
    "detected_genes_primary_quantile",
    "detected_genes_sensitivity_quantile",
    "RLE_primary_quantile",
    "RLE_sensitivity_quantile",
    "sample_correlation_primary_quantile",
    "sample_correlation_sensitivity_quantile",
    "PCA_robust_distance_alpha",
    "PCA_robust_distance_threshold",
    "filterByExpr_retained_genes",
    "VST_positive_variance_genes",
    "PCA_and_correlation_genes"
  ),
  probability = c(0.02, 0.005, 0.98, 0.995, 0.02, 0.005, 0.995, NA, NA, NA, NA),
  value = c(
    detected_primary_threshold,
    detected_sensitivity_threshold,
    rle_primary_threshold,
    rle_sensitivity_threshold,
    correlation_primary_threshold,
    correlation_sensitivity_threshold,
    pca_alpha,
    pca_threshold,
    length(retained_genes),
    nrow(vst_matrix),
    nrow(review_matrix)
  ),
  stringsAsFactors = FALSE
)

tool_versions <- data.frame(
  package = c("R", required_packages),
  version = c(
    as.character(getRversion()),
    vapply(
      required_packages,
      function(package_name) as.character(packageVersion(package_name)),
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)

readr::write_tsv(decisions, file.path(output_dir, "QC_decisions_per_sample.tsv"))
writeLines(
  sample_order[strict_qc_include],
  file.path(output_dir, "INCLUDE_strict_QC_FINAL.txt")
)
writeLines(
  sample_order[main_analysis_include],
  file.path(output_dir, "INCLUDE_main_analysis_FINAL.txt")
)
readr::write_tsv(flow_summary, file.path(output_dir, "QC_sample_flow_summary.tsv"))
readr::write_tsv(thresholds, file.path(output_dir, "QC_thresholds.tsv"))
readr::write_tsv(tool_versions, file.path(output_dir, "QC_tool_versions.tsv"))

session_connection <- file(file.path(output_dir, "QC_sessionInfo.txt"), open = "wt")
writeLines(capture.output(sessionInfo()), session_connection)
close(session_connection)

write_legacy_aliases <- tolower(Sys.getenv(
  "PANTHER_WRITE_LEGACY_QC_ALIASES",
  unset = "false"
)) %in% c("1", "true", "yes")
if (write_legacy_aliases) {
  writeLines(
    sample_order[strict_qc_include],
    file.path(output_dir, "INCLUDE_primary_FINAL.txt")
  )
  writeLines(
    sample_order[main_analysis_include],
    file.path(output_dir, "INCLUDE_sensitivity_FINAL.txt")
  )
}

cat("Input samples:", length(sample_order), "\n")
cat("Retained genes:", length(retained_genes), "\n")
cat("Strict QC included samples:", strict_counts[["samples"]], "\n")
cat("Main-analysis included samples:", main_counts[["samples"]], "\n")
