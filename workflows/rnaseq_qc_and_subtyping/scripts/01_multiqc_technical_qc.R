# Copyright (c) 2026 Emmanouil G. Sifakis
# SPDX-License-Identifier: MIT

# PANTHER RNA-seq technical QC: STAR input depth and unique mapping.
# Workflow author/maintainer: Emmanouil G. Sifakis.
# Inputs: MultiQC STAR table and required PANTHER_QCA_* environment variables.
# Outputs: technical decisions, ordered inclusion lists, thresholds, and counts.
# Execution: stage 1; run before expression QC from the repository root.

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("utils")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

required_environment <- c(
  "PANTHER_QCA_MULTIQC_DIR",
  "PANTHER_QCA_OUTPUT_DIR",
  "PANTHER_QCA_SAMPLE_ID_REGEX"
)
environment_values <- Sys.getenv(required_environment, unset = "")
missing_environment <- !nzchar(trimws(environment_values))
if (any(missing_environment)) {
  stop(
    "Missing required environment variables: ",
    paste(required_environment[missing_environment], collapse = ", "),
    call. = FALSE
  )
}

multiqc_dir <- normalizePath(
  environment_values[["PANTHER_QCA_MULTIQC_DIR"]],
  winslash = "/",
  mustWork = TRUE
)
output_dir <- environment_values[["PANTHER_QCA_OUTPUT_DIR"]]
sample_id_regex <- environment_values[["PANTHER_QCA_SAMPLE_ID_REGEX"]]

parse_optional_count <- function(name) {
  value <- trimws(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) {
    return(NA_integer_)
  }
  if (!grepl("^[1-9][0-9]*$", value)) {
    stop(name, " must be a positive integer", call. = FALSE)
  }
  as.integer(value)
}

expected_sample_count <- parse_optional_count(
  "PANTHER_QCA_EXPECTED_SAMPLE_COUNT"
)
expected_patient_count <- parse_optional_count(
  "PANTHER_QCA_EXPECTED_PATIENT_COUNT"
)

if (dir.exists(output_dir)) {
  existing_output <- list.files(
    output_dir,
    all.files = TRUE,
    no.. = TRUE,
    recursive = TRUE
  )
  if (length(existing_output)) {
    stop("Output directory is not empty: ", output_dir, call. = FALSE)
  }
}

resolve_exact_file <- function(directory, filename) {
  entries <- list.files(
    directory,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE,
    no.. = TRUE
  )
  matches <- entries[
    basename(entries) == filename &
      !dir.exists(entries)
  ]
  if (length(matches) != 1L) {
    stop(
      "Expected exactly one file named ",
      filename,
      " in the MultiQC directory; found ",
      length(matches),
      call. = FALSE
    )
  }
  normalizePath(matches, winslash = "/", mustWork = TRUE)
}

star_path <- resolve_exact_file(multiqc_dir, "multiqc_star.txt")

read_source_table <- function(path) {
  utils::read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
}

require_columns_once <- function(data, columns, label) {
  counts <- vapply(
    columns,
    function(column) sum(names(data) == column),
    integer(1)
  )
  if (any(counts != 1L)) {
    bad <- paste(columns[counts != 1L], counts[counts != 1L], sep = "=")
    stop(
      label,
      " must contain each required column exactly once: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
}

normalize_sample <- function(values, pattern) {
  vapply(
    as.character(values),
    function(value) {
      if (is.na(value)) {
        return(NA_character_)
      }
      matches <- regmatches(
        trimws(value),
        gregexpr(pattern, trimws(value), perl = TRUE)
      )[[1L]]
      if (length(matches) != 1L) {
        return(NA_character_)
      }
      matches[[1L]]
    },
    character(1),
    USE.NAMES = FALSE
  )
}

assert_unique_ids <- function(ids, label) {
  if (anyNA(ids) || any(!nzchar(ids))) {
    stop(label, " contains missing or ambiguous normalized IDs", call. = FALSE)
  }
  if (anyDuplicated(ids)) {
    stop(label, " contains duplicated normalized IDs", call. = FALSE)
  }
}

parse_numeric_strict <- function(values, label) {
  text_values <- trimws(as.character(values))
  missing <- is.na(values) |
    !nzchar(text_values) |
    toupper(text_values) == "NA"
  output <- rep(NA_real_, length(values))
  output[!missing] <- suppressWarnings(as.numeric(text_values[!missing]))
  invalid <- !missing & !is.finite(output)
  if (any(invalid)) {
    stop("Nonfinite or nonnumeric values in ", label, call. = FALSE)
  }
  output
}

qc_call_low <- function(value, fail_threshold, flag_threshold) {
  output <- rep(NA_character_, length(value))
  output[!is.na(value) & value < fail_threshold] <- "FAIL"
  output[
    !is.na(value) &
      value >= fail_threshold &
      value < flag_threshold
  ] <- "FLAG"
  output[!is.na(value) & value >= flag_threshold] <- "PASS"
  output
}

combine_hard_calls <- function(uniquely_mapped_call, depth_call) {
  hard_fail <-
    (!is.na(uniquely_mapped_call) & uniquely_mapped_call == "FAIL") |
    (!is.na(depth_call) & depth_call == "FAIL")
  hard_missing <- is.na(uniquely_mapped_call) | is.na(depth_call)
  hard_flag <-
    (!is.na(uniquely_mapped_call) & uniquely_mapped_call == "FLAG") |
    (!is.na(depth_call) & depth_call == "FLAG")
  output <- rep("PASS", length(uniquely_mapped_call))
  output[hard_flag] <- "FLAG"
  output[hard_missing] <- "FLAG"
  output[hard_fail] <- "FAIL"
  output
}

make_reasons <- function(
  uniquely_mapped_percentage,
  star_input_fragments_million,
  uniquely_mapped_call,
  depth_call
) {
  vapply(
    seq_along(uniquely_mapped_percentage),
    function(index) {
      reasons <- character()
      if (is.na(uniquely_mapped_call[[index]])) {
        reasons <- c(reasons, "Aligned=MISSING")
      } else if (uniquely_mapped_call[[index]] %in% c("FLAG", "FAIL")) {
        reasons <- c(
          reasons,
          sprintf(
            "Aligned=%s (%.1f%%)",
            uniquely_mapped_call[[index]],
            uniquely_mapped_percentage[[index]]
          )
        )
      }
      if (is.na(depth_call[[index]])) {
        reasons <- c(reasons, "Depth=MISSING")
      } else if (depth_call[[index]] %in% c("FLAG", "FAIL")) {
        reasons <- c(
          reasons,
          sprintf(
            "Depth=%s (%.1f M STAR input fragments)",
            depth_call[[index]],
            star_input_fragments_million[[index]]
          )
        )
      }
      paste(reasons, collapse = "; ")
    },
    character(1)
  )
}

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "NA",
    fileEncoding = "ASCII"
  )
}

star <- read_source_table(star_path)
require_columns_once(
  star,
  c("Sample", "uniquely_mapped_percent", "total_reads"),
  "STAR table"
)

# Extract one unambiguous study sample ID per row; retain STAR input order.
star_sample <- normalize_sample(star$Sample, sample_id_regex)
assert_unique_ids(star_sample, "STAR table")

uniquely_mapped_percentage <- parse_numeric_strict(
  star$uniquely_mapped_percent,
  "STAR uniquely_mapped_percent"
)
star_total_reads <- parse_numeric_strict(
  star$total_reads,
  "STAR total_reads"
)

if (
  any(!is.finite(uniquely_mapped_percentage)) ||
    any(uniquely_mapped_percentage < 0) ||
    any(uniquely_mapped_percentage > 100)
) {
  stop(
    "STAR uniquely_mapped_percent must be finite percentages from 0 to 100",
    call. = FALSE
  )
}
if (any(!is.finite(star_total_reads)) || any(star_total_reads < 0)) {
  stop("STAR total_reads must be finite and nonnegative", call. = FALSE)
}

if (
  !is.na(expected_sample_count) &&
    length(star_sample) != expected_sample_count
) {
  stop(
    "Observed sample count differs from configured expectation: ",
    length(star_sample),
    " versus ",
    expected_sample_count,
    call. = FALSE
  )
}

patient_id <- sub("^.*-", "", star_sample)
assert_unique_ids(patient_id, "Derived patient IDs")
if (
  !is.na(expected_patient_count) &&
    length(patient_id) != expected_patient_count
) {
  stop(
    "Observed patient count differs from configured expectation: ",
    length(patient_id),
    " versus ",
    expected_patient_count,
    call. = FALSE
  )
}

# These study-specific cutoffs define technical calls, before expression review.
uniquely_mapped_fail_threshold <- 35
uniquely_mapped_flag_threshold <- 50
star_input_fail_threshold_million <- 10
star_input_flag_threshold_million <- 15

# MultiQC STAR total_reads is the STAR Number of input reads field.
# Paired-end evidence supports one count per post-trim read pair or fragment.
# Individual mate-read totals would use different units and cannot be substituted.
star_total_reads_million <- star_total_reads / 1e6
uniquely_mapped_call <- qc_call_low(
  uniquely_mapped_percentage,
  uniquely_mapped_fail_threshold,
  uniquely_mapped_flag_threshold
)
depth_call <- qc_call_low(
  star_total_reads_million,
  star_input_fail_threshold_million,
  star_input_flag_threshold_million
)
qc_status <- combine_hard_calls(uniquely_mapped_call, depth_call)
# Stage A strict eligibility requires PASS; main eligibility also permits FLAG.
# Expression-level review applies separate rules in stage 2.
strict_qc_include <- qc_status == "PASS"
main_analysis_include <- qc_status %in% c("PASS", "FLAG")
qc_reasons <- make_reasons(
  uniquely_mapped_percentage,
  star_total_reads_million,
  uniquely_mapped_call,
  depth_call
)

qc_decisions <- data.frame(
  sample = star_sample,
  patient_id = patient_id,
  uniq_map_pct = uniquely_mapped_percentage,
  star_total_reads = star_total_reads,
  star_total_reads_million = star_total_reads_million,
  depth_million = star_total_reads_million,
  uniq_map_call = uniquely_mapped_call,
  depth_call = depth_call,
  qc_status = qc_status,
  strict_qc_include = strict_qc_include,
  main_analysis_include = main_analysis_include,
  include_primary = strict_qc_include,
  include_sensitivity = main_analysis_include,
  qc_reasons = qc_reasons,
  stringsAsFactors = FALSE
)

if (!identical(qc_decisions$sample, star_sample)) {
  stop("Output sample order differs from STAR input order", call. = FALSE)
}

summary_table <- data.frame(
  metric = c(
    "total_samples",
    "unique_nonmissing_patient_ids",
    "missing_patient_ids",
    "duplicated_patient_ids",
    "depth_source_records",
    "missing_depth_values",
    "technical_qc_pass",
    "technical_qc_flag",
    "technical_qc_fail",
    "strict_qc_included",
    "main_analysis_included"
  ),
  value = c(
    nrow(qc_decisions),
    length(unique(patient_id)),
    sum(is.na(patient_id) | !nzchar(patient_id)),
    sum(duplicated(patient_id)),
    sum(!is.na(star_total_reads)),
    sum(is.na(star_total_reads)),
    sum(qc_status == "PASS"),
    sum(qc_status == "FLAG"),
    sum(qc_status == "FAIL"),
    sum(strict_qc_include),
    sum(main_analysis_include)
  ),
  stringsAsFactors = FALSE
)

threshold_table <- data.frame(
  metric = c(
    "uniquely_mapped_percentage",
    "star_total_reads_million"
  ),
  source_filename = c(
    "multiqc_star.txt",
    "multiqc_star.txt"
  ),
  source_column = c(
    "uniquely_mapped_percent",
    "total_reads"
  ),
  unit = c(
    "percent",
    "million_post_trim_input_fragments"
  ),
  fail_below = c(
    uniquely_mapped_fail_threshold,
    star_input_fail_threshold_million
  ),
  flag_below = c(
    uniquely_mapped_flag_threshold,
    star_input_flag_threshold_million
  ),
  stringsAsFactors = FALSE
)

tool_versions <- data.frame(
  component = c("R", "utils"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    as.character(utils::packageVersion("utils"))
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

utils::write.csv(
  qc_decisions,
  file.path(output_dir, "multiqc_qc_calls_and_inclusion.csv"),
  row.names = FALSE,
  na = "NA",
  fileEncoding = "ASCII"
)
writeLines(
  qc_decisions$sample[qc_decisions$qc_status == "FAIL"],
  file.path(output_dir, "FAILED_samples.txt"),
  useBytes = TRUE
)
writeLines(
  qc_decisions$sample[qc_decisions$qc_status == "FLAG"],
  file.path(output_dir, "FLAGGED_samples.txt"),
  useBytes = TRUE
)
writeLines(
  qc_decisions$sample[qc_decisions$strict_qc_include],
  file.path(output_dir, "INCLUDE_primary_PASS_only.txt"),
  useBytes = TRUE
)
writeLines(
  qc_decisions$sample[qc_decisions$main_analysis_include],
  file.path(output_dir, "INCLUDE_sensitivity_PASS_or_FLAG.txt"),
  useBytes = TRUE
)
write_tsv(
  summary_table,
  file.path(output_dir, "QC_sample_flow_summary.tsv")
)
write_tsv(
  threshold_table,
  file.path(output_dir, "QC_thresholds.tsv")
)
write_tsv(
  tool_versions,
  file.path(output_dir, "QC_tool_versions.tsv")
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(output_dir, "QC_sessionInfo.txt"),
  useBytes = TRUE
)

message("Technical QC completed.")
message("Samples: ", nrow(qc_decisions))
message("Unique nonmissing patient IDs: ", length(unique(patient_id)))
message("Strict QC included: ", sum(strict_qc_include))
message("Main analysis included: ", sum(main_analysis_include))
