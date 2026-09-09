# Copyright (c) 2026 Emmanouil G. Sifakis
# SPDX-License-Identifier: MIT

# PANTHER SSPBC intrinsic subtyping and RNA-seq ROR.
# Workflow author/maintainer: Emmanouil G. Sifakis.
# Inputs: matching gene-by-sample FPKM/count RDS matrices and final main list.
# Outputs: ordered subtype calls, probabilities, ROR, and mapping coverage.
# Execution: stage 3, after expression QC; run from the repository root.

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

required_packages <- c(
  "AnnotationDbi",
  "BreastSubtypeR",
  "org.Hs.eg.db",
  "S4Vectors",
  "sspbc",
  "SummarizedExperiment"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

fpkm_file <- Sys.getenv("PANTHER_SSPBC_FPKM_RDS", unset = "")
counts_file <- Sys.getenv("PANTHER_SSPBC_COUNTS_RDS", unset = "")
include_file <- Sys.getenv("PANTHER_MAIN_ANALYSIS_INCLUDE_FILE", unset = "")
output_dir <- Sys.getenv("PANTHER_SSPBC_OUTPUT_DIR", unset = "")

configured_paths <- c(fpkm_file, counts_file, include_file, output_dir)
if (any(!nzchar(configured_paths))) {
  stop(
    "Set PANTHER_SSPBC_FPKM_RDS, PANTHER_SSPBC_COUNTS_RDS, ",
    "PANTHER_MAIN_ANALYSIS_INCLUDE_FILE, and PANTHER_SSPBC_OUTPUT_DIR."
  )
}
if (any(!file.exists(c(fpkm_file, counts_file, include_file)))) {
  stop("One or more required input files do not exist.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(
    x,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA",
    fileEncoding = "ASCII"
  )
}

normalize_sample_id <- function(x) {
  normalized <- trimws(as.character(x))
  normalized[!nzchar(normalized)] <- NA_character_
  normalized
}

patient_id_from_sample <- function(x) {
  sub("^.*-", "", x, perl = TRUE)
}

extract_one_column <- function(x, expected_samples, label) {
  if (is.matrix(x) || is.data.frame(x)) {
    if (ncol(x) != 1L || is.null(rownames(x))) {
      stop(label, " must contain one named value per sample.")
    }
    values <- x[, 1L]
    names(values) <- rownames(x)
  } else {
    values <- x
  }
  if (is.null(names(values)) || anyDuplicated(names(values))) {
    stop(label, " names are missing or duplicated.")
  }
  result_index <- match(expected_samples, names(values))
  if (anyNA(result_index)) {
    stop(
      label,
      " is missing samples: ",
      paste(expected_samples[is.na(result_index)], collapse = ", ")
    )
  }
  result <- values[result_index]
  names(result) <- expected_samples
  result
}

load_models <- function(package_name) {
  model_environment <- new.env(parent = emptyenv())
  utils::data(
    list = "sspbc.models",
    package = package_name,
    envir = model_environment
  )
  model_environment[["sspbc.models"]]
}

required_model_genes <- function(package_name, model_name) {
  model <- load_models(package_name)[[model_name]]
  pair_members <- unlist(strsplit(model$all.pairs, "<", fixed = TRUE))
  sort(unique(sub("^e", "", pair_members)))
}

# These source symbols are aliases for model-required Entrez genes.
manual_entrez <- c(
  C2orf54 = "79919",
  SGOL1 = "151648",
  C2orf40 = "84417",
  C20orf194 = "25943",
  LHFP = "10186",
  COL4A3BP = "10087",
  C10orf32 = "119032",
  GLTSCR1L = "23506",
  ZMYM6NB = "100506144",
  MARCH8 = "220972",
  FAM63A = "55793",
  CCDC115 = "84317",
  FAM122A = "116224",
  FAM83D = "81610"
)

fpkm <- readRDS(fpkm_file)
counts <- readRDS(counts_file)
if (!is.matrix(fpkm) || !is.numeric(fpkm)) {
  stop("The FPKM input must be a numeric matrix.")
}
if (!is.matrix(counts) || !is.numeric(counts)) {
  stop("The count input must be a numeric matrix.")
}
if (any(!is.finite(counts)) || any(counts < 0)) {
  stop("Count values must be finite and nonnegative.")
}
if (is.null(rownames(fpkm)) || anyNA(rownames(fpkm)) ||
    any(!nzchar(rownames(fpkm))) || anyDuplicated(rownames(fpkm))) {
  stop("FPKM gene identifiers are missing or duplicated.")
}
if (!identical(rownames(fpkm), rownames(counts))) {
  stop("FPKM and count gene rows differ.")
}
if (!identical(colnames(fpkm), colnames(counts))) {
  stop("FPKM and count sample columns differ.")
}
if (any(!is.finite(fpkm)) || any(fpkm < 0)) {
  stop("FPKM values must be finite and nonnegative.")
}

fpkm_samples <- normalize_sample_id(colnames(fpkm))
counts_samples <- normalize_sample_id(colnames(counts))
if (anyNA(fpkm_samples) || anyDuplicated(fpkm_samples)) {
  stop("FPKM sample identifiers are missing or duplicated after normalization.")
}
if (anyNA(counts_samples) || anyDuplicated(counts_samples)) {
  stop("Count sample identifiers are missing or duplicated after normalization.")
}
if (!identical(fpkm_samples, counts_samples)) {
  stop("FPKM and count sample order differs after normalization.")
}

# The preceding QC stage supplies membership. This exact-444 safeguard checks
# the confirmed analysis set size and unique patients; it does not select samples.
main_samples <- normalize_sample_id(readLines(include_file, warn = FALSE))
if (length(main_samples) != 444L || anyNA(main_samples) ||
    anyDuplicated(main_samples)) {
  stop("The inclusion list must contain exactly 444 unique normalized samples.")
}
main_patients <- patient_id_from_sample(main_samples)
if (anyNA(main_patients) || any(!nzchar(main_patients)) ||
    any(!grepl("^[0-9]{5}$", main_patients)) || anyDuplicated(main_patients)) {
  stop("The inclusion list must contain exactly 444 unique patient identifiers.")
}

# Match by sample ID to preserve inclusion-file order in both expression inputs.
fpkm_index <- match(main_samples, fpkm_samples)
counts_index <- match(main_samples, counts_samples)
missing_from_fpkm <- main_samples[is.na(fpkm_index)]
missing_from_counts <- main_samples[is.na(counts_index)]
message(
  "Samples missing from the FPKM matrix: ",
  if (length(missing_from_fpkm) == 0L) "none" else
    paste(missing_from_fpkm, collapse = ", ")
)
message(
  "Samples missing from the count matrix: ",
  if (length(missing_from_counts) == 0L) "none" else
    paste(missing_from_counts, collapse = ", ")
)
if (length(missing_from_fpkm) > 0L || length(missing_from_counts) > 0L) {
  stop("The expression matrices do not contain the complete inclusion list.")
}

not_in_main_analysis <- fpkm_samples[!fpkm_samples %in% main_samples]
message(
  "Expression samples outside the main-analysis list: ",
  length(not_in_main_analysis)
)

fpkm_main <- fpkm[, fpkm_index, drop = FALSE]
counts_main <- counts[, counts_index, drop = FALSE]
colnames(fpkm_main) <- main_samples
colnames(counts_main) <- main_samples
if (!identical(colnames(fpkm_main), main_samples) ||
    !identical(colnames(counts_main), main_samples)) {
  stop("Expression matrices do not preserve inclusion-list order.")
}

annotation_database <- get(
  "org.Hs.eg.db",
  envir = asNamespace("org.Hs.eg.db")
)
automatic_entrez <- AnnotationDbi::mapIds(
  annotation_database,
  keys = rownames(fpkm_main),
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
# Every listed alias must resolve to its recorded Entrez ID in this annotation
# version; stop on unsupported mappings rather than substituting other genes.
alias_entrez <- AnnotationDbi::mapIds(
  annotation_database,
  keys = names(manual_entrez),
  column = "ENTREZID",
  keytype = "ALIAS",
  multiVals = "first"
)
if (anyNA(alias_entrez) ||
    !identical(unname(alias_entrez[names(manual_entrez)]), unname(manual_entrez))) {
  stop("A required manual Entrez mapping is not supported by the annotation package.")
}

feature_annotation <- data.frame(
  probe = rownames(fpkm_main),
  SYMBOL = rownames(fpkm_main),
  ENTREZID = unname(automatic_entrez),
  stringsAsFactors = FALSE
)
manual_index <- match(names(manual_entrez), feature_annotation$SYMBOL)
if (anyNA(manual_index)) {
  stop(
    "Required mapping symbols are missing from the FPKM matrix: ",
    paste(names(manual_entrez)[is.na(manual_index)], collapse = ", ")
  )
}
feature_annotation$ENTREZID[manual_index] <- unname(manual_entrez)

subtype_model_genes <- required_model_genes(
  "BreastSubtypeR",
  "ssp.subtype"
)
ror_model_genes <- required_model_genes("sspbc", "ssp.ror")

sample_metadata <- data.frame(
  sample = main_samples,
  PatientID = main_patients,
  stringsAsFactors = FALSE
)

# Mapping expects log2-normalized FPKM when RawCounts is FALSE.
subtype_input <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = log2(fpkm_main + 1)),
  rowData = S4Vectors::DataFrame(feature_annotation),
  colData = S4Vectors::DataFrame(sample_metadata)
)
# In BreastSubtypeR 1.4.0, impute = FALSE disables KNN replacement of NA
# matrix values. It does not supply absent model-gene rows; coverage is checked
# explicitly before prediction. This setting is separate from input scaling.
mapped <- BreastSubtypeR::Mapping(
  se_obj = subtype_input,
  RawCounts = FALSE,
  method = "max",
  impute = FALSE,
  verbose = FALSE
)
subtype_expression <- SummarizedExperiment::assay(mapped$se_SSP)
# Require complete subtype model-gene coverage in the actual mapped input.
missing_subtype_genes <- setdiff(subtype_model_genes, rownames(subtype_expression))
if (length(missing_subtype_genes) > 0L) {
  stop(
    "Missing ssp.subtype model genes: ",
    paste(missing_subtype_genes, collapse = ", ")
  )
}

subtype_prediction <- BreastSubtypeR::BS_sspbc(
  se_obj = mapped$se_SSP,
  ssp.name = "ssp.subtype"
)
subtype_call <- as.character(extract_one_column(
  subtype_prediction$cl,
  main_samples,
  "Subtype call"
))
subtype_probability <- as.numeric(extract_one_column(
  subtype_prediction$prob,
  main_samples,
  "Subtype probability"
))
names(subtype_probability) <- main_samples
allowed_subtypes <- c("LumA", "LumB", "Her2", "Basal")
if (anyNA(subtype_call) || any(!subtype_call %in% allowed_subtypes)) {
  stop("The subtype output contains missing or unexpected classes.")
}
if (any(!is.finite(subtype_probability))) {
  stop("The subtype output contains non-finite probabilities.")
}

# ROR uses nonnegative linear StringTie FPKM values.
# Keep this scale separate from the log2(FPKM + 1) subtype input above.
# When several rows map to one Entrez ID, retain the row with the largest sum.
mapped_rows <- !is.na(feature_annotation$ENTREZID)
ror_annotation <- feature_annotation[mapped_rows, , drop = FALSE]
ror_expression <- fpkm_main[mapped_rows, , drop = FALSE]
ror_annotation$row_sum <- rowSums(ror_expression)
ror_order <- order(ror_annotation$ENTREZID, -ror_annotation$row_sum)
ror_annotation <- ror_annotation[ror_order, , drop = FALSE]
ror_expression <- ror_expression[ror_order, , drop = FALSE]
unique_entrez <- !duplicated(ror_annotation$ENTREZID)
ror_annotation <- ror_annotation[unique_entrez, , drop = FALSE]
ror_expression <- ror_expression[unique_entrez, , drop = FALSE]
rownames(ror_expression) <- ror_annotation$ENTREZID
if (anyDuplicated(rownames(ror_expression))) {
  stop("Duplicate Entrez identifiers remain in the ROR matrix.")
}
# Require complete ROR model-gene coverage after duplicate Entrez resolution.
missing_ror_genes <- setdiff(ror_model_genes, rownames(ror_expression))
if (length(missing_ror_genes) > 0L) {
  stop(
    "Missing ssp.ror model genes: ",
    paste(missing_ror_genes, collapse = ", ")
  )
}

# sspbc 1.0 resolves these package data objects during applySSP execution.
utils::data("Gene.ID.ann", package = "sspbc", envir = .GlobalEnv)
utils::data("sspbc.models", package = "sspbc", envir = .GlobalEnv)
utils::data("sspbc.models.fullname", package = "sspbc", envir = .GlobalEnv)

ror_prediction <- sspbc::applySSP(
  gex = as.matrix(ror_expression),
  id = rownames(ror_expression),
  id.type = "EntrezGene",
  ssp.name = "ssp.ror",
  full.out = TRUE
)
ror_class <- as.character(extract_one_column(
  ror_prediction$cl,
  main_samples,
  "ROR class"
))
ror_score <- as.numeric(sub("^c0*", "", ror_class))
names(ror_score) <- main_samples
if (anyNA(ror_score) || any(!is.finite(ror_score))) {
  stop("The ROR output contains missing or non-finite values.")
}

result <- data.frame(
  sample = main_samples,
  patient_id = main_patients,
  sspbc_intrinsic_subtype = subtype_call,
  sspbc_subtype_probability = subtype_probability,
  sspbc_ror = ror_score,
  stringsAsFactors = FALSE
)
if (nrow(result) != 444L || anyDuplicated(result$sample) ||
    anyDuplicated(result$patient_id) || !identical(result$sample, main_samples)) {
  stop("Final SSPBC output alignment validation failed.")
}
write_tsv(
  result,
  file.path(output_dir, "PANTHER_SSPBC_subtype_and_ROR_main_analysis.tsv")
)

automatic_present_ids <- unique(na.omit(as.character(automatic_entrez)))
final_present_ids <- unique(na.omit(as.character(feature_annotation$ENTREZID)))
mapping_summary <- data.frame(
  category = c(
    rep("input", 8L),
    rep("mapping", 7L),
    rep("model_coverage", 6L),
    rep("manual_mapping", length(manual_entrez))
  ),
  item = c(
    "fpkm_gene_rows",
    "fpkm_expression_samples",
    "count_expression_samples",
    "main_analysis_samples",
    "main_analysis_unique_patients",
    "samples_missing_from_fpkm",
    "samples_missing_from_counts",
    "expression_samples_outside_main_analysis",
    "automatically_mapped_gene_rows",
    "unmapped_gene_rows_before_manual_mapping",
    "duplicate_entrez_rows_before_manual_mapping",
    "mapped_gene_rows_after_manual_mapping",
    "unmapped_gene_rows_after_manual_mapping",
    "duplicate_entrez_rows_after_manual_mapping",
    "manual_mapping_count",
    "ssp_subtype_required_genes",
    "ssp_subtype_required_genes_present",
    "ssp_subtype_required_genes_missing",
    "ssp_ror_required_genes",
    "ssp_ror_required_genes_present",
    "ssp_ror_required_genes_missing",
    names(manual_entrez)
  ),
  value = as.character(c(
    nrow(fpkm_main),
    length(fpkm_samples),
    length(counts_samples),
    length(main_samples),
    length(unique(main_patients)),
    length(missing_from_fpkm),
    length(missing_from_counts),
    length(not_in_main_analysis),
    sum(!is.na(automatic_entrez)),
    sum(is.na(automatic_entrez)),
    sum(duplicated(automatic_entrez[!is.na(automatic_entrez)])),
    sum(!is.na(feature_annotation$ENTREZID)),
    sum(is.na(feature_annotation$ENTREZID)),
    sum(duplicated(
      feature_annotation$ENTREZID[!is.na(feature_annotation$ENTREZID)]
    )),
    length(manual_entrez),
    length(subtype_model_genes),
    length(subtype_model_genes) - length(missing_subtype_genes),
    length(missing_subtype_genes),
    length(ror_model_genes),
    length(ror_model_genes) - length(missing_ror_genes),
    length(missing_ror_genes),
    unname(manual_entrez)
  )),
  detail = c(
    rep("", 21L),
    paste0("alias_to_entrez;annotation_alias=", unname(alias_entrez))
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  mapping_summary,
  file.path(output_dir, "PANTHER_SSPBC_gene_mapping_summary.tsv")
)

all_model_genes <- unique(c(subtype_model_genes, ror_model_genes))
current_symbols <- AnnotationDbi::mapIds(
  annotation_database,
  keys = all_model_genes,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
manual_symbol_by_entrez <- setNames(names(manual_entrez), manual_entrez)
coverage <- rbind(
  data.frame(
    model = "ssp.subtype",
    entrez_id = subtype_model_genes,
    current_symbol = unname(current_symbols[subtype_model_genes]),
    manual_source_symbol = unname(
      manual_symbol_by_entrez[subtype_model_genes]
    ),
    automatic_mapping_present = subtype_model_genes %in%
      automatic_present_ids,
    final_mapping_present = subtype_model_genes %in% final_present_ids,
    model_input_present = subtype_model_genes %in%
      rownames(subtype_expression),
    stringsAsFactors = FALSE
  ),
  data.frame(
    model = "ssp.ror",
    entrez_id = ror_model_genes,
    current_symbol = unname(current_symbols[ror_model_genes]),
    manual_source_symbol = unname(manual_symbol_by_entrez[ror_model_genes]),
    automatic_mapping_present = ror_model_genes %in%
      automatic_present_ids,
    final_mapping_present = ror_model_genes %in% final_present_ids,
    model_input_present = ror_model_genes %in% rownames(ror_expression),
    stringsAsFactors = FALSE
  )
)
write_tsv(
  coverage,
  file.path(output_dir, "PANTHER_SSPBC_model_gene_coverage.tsv")
)

tool_versions <- data.frame(
  tool = c("R", required_packages),
  version = c(
    as.character(getRversion()),
    vapply(
      required_packages,
      function(package_name) as.character(utils::packageVersion(package_name)),
      character(1)
    )
  ),
  stringsAsFactors = FALSE
)
write_tsv(
  tool_versions,
  file.path(output_dir, "PANTHER_SSPBC_tool_versions.tsv")
)
writeLines(
  capture.output(utils::sessionInfo()),
  con = file.path(output_dir, "PANTHER_SSPBC_sessionInfo.txt"),
  useBytes = TRUE
)

message("SSPBC subtype calls: ", nrow(result))
message("SSPBC ROR values: ", sum(!is.na(result$sspbc_ror)))
message("SSPBC output order matches the inclusion list: TRUE")
