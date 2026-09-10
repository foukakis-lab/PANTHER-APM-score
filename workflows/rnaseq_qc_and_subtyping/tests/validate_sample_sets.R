# Copyright (c) 2026 Emmanouil G. Sifakis
# SPDX-License-Identifier: MIT

# Validate final sample-set files without loading participant metadata.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript validate_sample_sets.R STRICT_FILE MAIN_FILE")
}

read_ids <- function(path, label) {
  if (!file.exists(path)) stop(label, " file does not exist: ", path)
  values <- trimws(readLines(path, warn = FALSE, encoding = "ASCII"))
  values <- values[nzchar(values)]
  if (anyDuplicated(values)) stop(label, " contains duplicated identifiers.")
  values
}

strict_ids <- read_ids(args[[1L]], "Strict set")
main_ids <- read_ids(args[[2L]], "Main-analysis set")

derive_patient_ids <- function(sample_ids, label) {
  patient_ids <- sub("^.*-", "", sample_ids)
  if (any(!grepl("^[0-9]{5}$", patient_ids))) {
    stop(label, " contains an identifier without a five-digit patient suffix.")
  }
  if (anyDuplicated(patient_ids)) {
    stop(label, " contains duplicated patient identifiers.")
  }
  patient_ids
}

strict_patient_ids <- derive_patient_ids(strict_ids, "Strict set")
main_patient_ids <- derive_patient_ids(main_ids, "Main-analysis set")

if (length(strict_ids) != 301L) stop("Strict set does not contain 301 identifiers.")
if (length(main_ids) != 444L) stop("Main-analysis set does not contain 444 identifiers.")
if (length(strict_patient_ids) != 301L) stop("Strict set does not contain 301 patients.")
if (length(main_patient_ids) != 444L) stop("Main-analysis set does not contain 444 patients.")
if (!all(strict_ids %in% main_ids)) stop("Strict set is not a subset of the main-analysis set.")

message("Sample-set validation passed.")
