#!/usr/bin/env Rscript
# SPDX-License-Identifier: MIT
# Reproduce the corrected four-panel S13 from authorized study inputs.
# From the repository root: Rscript --vanilla workflows/rnaseq_qc_and_subtyping/scripts/04_supplementary_figure_S13.R
# Input contract and tested dependencies: ../docs/supplementary_figure_S13.md

stop_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

input_path <- function(variable) {
  value <- Sys.getenv(variable, unset = "")
  stop_if(!nzchar(value), paste("Set", variable))
  stop_if(!file.exists(value) || dir.exists(value), paste("Missing file for", variable))
  normalizePath(value, winslash = "/", mustWork = TRUE)
}

sample_table <- function(path, required) {
  # Character keys preserve leading zeros; do not log sample names.
  x <- read.delim(path, colClasses = "character", check.names = FALSE)
  stop_if(!all(c("sample", required) %in% names(x)), "Missing required table columns.")
  stop_if(anyNA(x$sample) || any(!nzchar(x$sample)) || anyDuplicated(x$sample) > 0L,
          "Expected one nonmissing row per sample.")
  x
}

align <- function(x, ids) {
  stop_if(!setequal(x$sample, ids), "Sample-ID sets differ.")
  x <- x[match(ids, x$sample), , drop = FALSE]
  stopifnot(identical(x$sample, ids))
  x
}

preseq_at_442 <- function(path, ids) {
  # Each MultiQC curve block supplies its own total-read grid, in millions.
  grid <- numeric()
  estimates <- numeric()
  seen <- character()
  for (line in readLines(path, warn = FALSE)) {
    if (!nzchar(trimws(line))) next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (!nzchar(trimws(fields[1L]))) {
      grid <- suppressWarnings(as.numeric(fields[-1L]))
      next
    }
    id <- sub("\\.lc_extrap$", "", trimws(fields[1L]))
    if (!id %in% ids) next  # e.g. the perfect-library reference curve
    stop_if(id %in% seen, "Duplicate preseq curve for a cohort sample.")
    seen <- c(seen, id)
    depth <- which(grid == 442)
    stop_if(length(depth) > 1L, "Repeated 442-million coordinate.")
    if (!length(depth) || length(fields) < depth + 1L) next
    y <- suppressWarnings(as.numeric(fields[depth + 1L]))
    if (is.finite(y)) estimates[id] <- y
  }
  stop_if(!length(estimates) || any(estimates < 0), "Invalid preseq input at the required depth.")
  estimates
}

named_leaves <- function(counts, decisions) {
  # Recover IDs only. The submitted heatmap, dendrograms and axes are never redrawn.
  ids <- decisions$sample
  x <- round(counts)
  keep <- edgeR::filterByExpr(x, group = rep(1L, length(ids)), min.count = 10)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    x[keep, , drop = FALSE], S4Vectors::DataFrame(row.names = ids), ~ 1
  )
  rm(x)
  m <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  rv <- matrixStats::rowVars(m)
  top <- head(order(rv, decreasing = TRUE), min(5000L, sum(rv > 0)))
  m <- m[top, , drop = FALSE]
  stopifnot(all(is.finite(m)), identical(colnames(m), ids))
  # With no missing values, this equals pairwise-complete Spearman.
  cc <- stats::cor(m, method = "spearman", use = "everything")
  med <- vapply(seq_along(ids), function(i) median(cc[-i, i]), numeric(1))
  stopifnot(max(abs(med - as.numeric(decisions$median_sample_correlation))) < 1e-10)
  hc <- stats::hclust(stats::dist(cc), method = "complete")
  leaves <- hc$labels[hc$order]
  stopifnot(identical(leaves, ids[hc$order]))
  leaves
}

draw_preseq <- function(values, decisions, file) {
  qc <- decisions$technical_qc_status[match(names(values), decisions$sample)]
  classes <- c("PASS", "FLAG", "FAIL")
  stopifnot(all(qc %in% classes))
  a <- data.frame(class = factor(qc, levels = classes), estimate = unname(values))
  p <- ggplot2::ggplot(a, ggplot2::aes(class, estimate, fill = class)) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.8) +
    ggplot2::geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
    ggplot2::scale_fill_manual(values = c(PASS="#777777", FLAG="#E69F00", FAIL="#CC79A7"), drop = FALSE) +
    ggplot2::labs(x = "MultiQC class",
      y = "Preseq-estimated distinct reads (millions)\nat 442 million total reads") +
    ggplot2::theme_bw(base_family = "Arial", base_size = 12) +
    ggplot2::theme(legend.position = "none", plot.margin = ggplot2::margin(4, 7, 4, 3))
  ragg::agg_png(file, width = 1890, height = 1210, units = "px", res = 218, background = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
}

main <- function() {
  required <- c("ggplot2", "ragg", "DESeq2", "edgeR", "S4Vectors", "SummarizedExperiment", "matrixStats")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  stop_if(length(missing) > 0L, paste("Required packages are unavailable:", paste(missing, collapse = ", ")))
  helper <- "workflows/rnaseq_qc_and_subtyping/scripts/s13_render.py"
  stop_if(!file.exists(helper), "Run from the repository root with the renderer present.")
  py <- Sys.getenv("PANTHER_S13_PYTHON", unset = "python")
  if (!file.exists(py)) py <- Sys.which(py)
  stop_if(!nzchar(py) || !file.exists(py), "Set PANTHER_S13_PYTHON to an available Python executable.")
  out <- Sys.getenv("PANTHER_S13_OUTPUT_DIR", unset = "")
  stop_if(!nzchar(out) || !dir.exists(out), "Provide an existing empty private output directory.")
  out <- normalizePath(out, winslash = "/")
  root <- normalizePath(".", winslash = "/")
  # Case-insensitive guard also covers Windows paths with different capitalization.
  stop_if(startsWith(paste0(tolower(out), "/"), paste0(tolower(root), "/")),
          "Use an output directory outside the repository.")
  stop_if(length(list.files(out, all.files = TRUE, no.. = TRUE)) > 0L, "Output directory is not empty.")
  canvas <- input_path("PANTHER_S13_ORIGINAL")
  reference <- input_path("PANTHER_S13_CORRECTED_REFERENCE")
  preseq_file <- input_path("PANTHER_S13_PRESEQ")
  tx <- readRDS(input_path("PANTHER_TXIMPORT_RDS"))
  stopifnot(identical(tx$countsFromAbundance, "lengthScaledTPM"))
  counts <- tx$counts
  rm(tx)
  ids <- colnames(counts)
  stopifnot(is.matrix(counts), all(is.finite(counts)), !is.null(ids),
            !anyNA(ids), all(nzchar(ids)), !anyDuplicated(ids))
  decisions <- align(sample_table(input_path("PANTHER_S13_DECISIONS"),
    c("technical_qc_status", "detected_genes_counts_ge10", "median_sample_correlation", "pca_review_flag")), ids)
  historical <- align(sample_table(input_path("PANTHER_S13_HISTORICAL_MASTER"),
    c("detected_features_call", "library_complexity")), ids)
  stopifnot(all(historical$detected_features_call %in% c("PASS","FLAG")),
            all(decisions$pca_review_flag %in% c("TRUE","FALSE")),
            all(decisions$technical_qc_status %in% c("PASS","FLAG","FAIL")))

  # Panel D: canonical values, checked on the full unrounded matrix. This is
  # the descriptive type-7 lower decile, not the main-inclusion strong flag.
  detected <- as.numeric(decisions$detected_genes_counts_ge10)
  stopifnot(all(colSums(counts >= 10) == detected))
  cutoff <- unname(stats::quantile(detected, 0.10, type = 7))
  flag <- detected < cutoff

  # Panel A: preseq only. No missing estimate restricts any other panel.
  values <- preseq_at_442(preseq_file, ids)
  matched <- match(names(values), historical$sample)
  stopifnot(all(abs(as.numeric(historical$library_complexity[matched]) - values) < 1e-9))
  leaves <- named_leaves(counts, decisions)
  ii <- match(leaves, ids)
  stopifnot(identical(ids[ii], leaves), !anyDuplicated(ii), length(ii) == length(ids))
  tracks <- data.frame(
    historical = as.integer(historical$detected_features_call[ii] == "FLAG"),
    corrected = as.integer(flag[ii]),
    pca = as.integer(decisions$pca_review_flag[ii] == "TRUE"),
    qc = match(decisions$technical_qc_status[ii], c("PASS","FLAG","FAIL")) - 1L
  )
  # Temporary renderer inputs contain plot pixels and ordered integer classes only.
  panel <- tempfile(fileext = ".png")
  track_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(panel, track_file)), add = TRUE)
  draw_preseq(values, decisions, panel)
  utils::write.table(tracks, track_file, sep = "\t", row.names = FALSE, quote = FALSE)
  argv <- c("-B", helper, canvas, reference, panel, track_file, out,
            as.character(length(values)), format(cutoff, scientific = FALSE, trim = TRUE))
  status <- system2(py, args = vapply(argv, shQuote, character(1)))
  stop_if(status != 0L, "S13 renderer failed; no candidate should be promoted.")
  message(sprintf("Preseq: %d/%d available. Detect: %d flags; %d changes. No inclusion decisions modified.",
                  length(values), length(ids), sum(flag), sum(tracks$historical != tracks$corrected)))
}
main()
