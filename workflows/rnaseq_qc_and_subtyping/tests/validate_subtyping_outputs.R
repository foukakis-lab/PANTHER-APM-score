# Copyright (c) 2026 Emmanouil G. Sifakis
# SPDX-License-Identifier: MIT

# Validate SSPBC output alignment and completeness.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript validate_subtyping_outputs.R INCLUDE_FILE OUTPUT_TSV")
}

include_ids <- trimws(readLines(args[[1L]], warn = FALSE, encoding = "ASCII"))
include_ids <- include_ids[nzchar(include_ids)]
results <- utils::read.delim(
  args[[2L]],
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "sample",
  "patient_id",
  "sspbc_intrinsic_subtype",
  "sspbc_subtype_probability",
  "sspbc_ror"
)
if (!all(required_columns %in% names(results))) {
  stop("SSPBC output is missing required columns.")
}
if (length(include_ids) != 444L || nrow(results) != 444L) {
  stop("Expected 444 inclusion identifiers and 444 output rows.")
}
if (anyDuplicated(include_ids) || anyDuplicated(results$sample)) {
  stop("Duplicated sample identifiers detected.")
}
if (anyDuplicated(results$patient_id) || anyNA(results$patient_id)) {
  stop("Missing or duplicated patient identifiers detected.")
}
if (!identical(include_ids, results$sample)) {
  stop("SSPBC output order differs from inclusion-file order.")
}
if (anyNA(results$sspbc_intrinsic_subtype) ||
    anyNA(results$sspbc_subtype_probability) ||
    anyNA(results$sspbc_ror)) {
  stop("SSPBC output contains missing subtype, probability, or ROR values.")
}
allowed_subtypes <- c("LumA", "LumB", "Her2", "Basal")
if (any(!results$sspbc_intrinsic_subtype %in% allowed_subtypes)) {
  stop("SSPBC output contains an unexpected subtype label.")
}
probability <- suppressWarnings(as.numeric(results$sspbc_subtype_probability))
ror <- suppressWarnings(as.numeric(results$sspbc_ror))
if (any(!is.finite(probability)) || any(probability < 0 | probability > 1)) {
  stop("SSPBC subtype probabilities must be finite and between zero and one.")
}
if (any(!is.finite(ror))) {
  stop("SSPBC ROR values must be finite.")
}

message("SSPBC output validation passed.")
