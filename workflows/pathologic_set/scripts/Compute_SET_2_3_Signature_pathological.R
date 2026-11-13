# ==============================================================================
# RNA-seq-adapted pathologic SET2,3 calculation
#
# Purpose:
#   Calculate SETER/PR, continuous RNA4 risk score, pathologic BPI, and
#   SET2,3 using RNA-seq expression data and surgical pathology information.
#
# Required function inputs:
#   gene_exp_data
#       Gene-expression matrix with gene symbols as rownames and sample IDs
#       as column names. Values are assumed to be log2-scale RNA-seq
#       expression values, for example log2(TMM-normalized TPM + 1).
#
#   clinical_data
#       Clinical data frame for the cohort to be scored.
#
# Required clinical columns:
#   sample_id
#       Sample identifiers matching gene_exp_data column names.
#
#   tumdiam
#       Pathologic invasive tumour diameter in millimetres.
#
#   posnodes
#       Number of pathologically involved regional lymph nodes.
#
# Required whole-cohort reference objects in the global environment:
#   tmm_normalized_tpm_symbols
#       RNA-seq expression matrix for the whole PANTHER cohort.
#
#   all_clinical
#       Clinical data frame for the whole PANTHER cohort with sample_id.
#
# Why whole-cohort objects are required:
#   PGR high/low status is derived from the RNA4-scaled PGR distribution using
#   a two-component Gaussian mixture model. The threshold is estimated in the
#   whole PANTHER cohort, where the relevant bimodal expression structure is
#   retained, and then applied unchanged to the target cohort.
#
# Components directly following Bossuyt et al. 2021:
#   SET_ER_PR:
#       mean expression of 18 ER/PR-related genes minus mean expression of 10 reference genes.
#
#   RNA4 scaling:
#       RNA4 gene expression - mean expression of 10 reference genes + 10.
#
#   pT risk score:
#       0 if pathologic tumour diameter <1 cm;
#       3 if pathologic tumour diameter >3.5 cm;
#       otherwise linearly scaled from 0 to 3 in increments of 0.1.
#
#   pN risk score:
#       5 if >10 involved nodes;
#       otherwise linearly scaled from 0 to 5 in increments of 0.5.
#
#   RNA4 risk score:
#       Continuous score ranging from 0 to 3, determined by AURKA expression
#       in the context of PGR status.
#
#   BPI:
#       BPI = (11 - (pT risk + pN risk + RNA4 risk)) * 4 / 11.
#
#   SET2,3:
#       SET2,3 = 0.5056 * BPI + 0.7542 * SET_ER_PR.
#
# Explicit RNA-seq adaptations used here:
#   1. PGR high/low threshold:
#       The published fixed QGP PGR threshold cannot be directly transferred
#       to uncalibrated RNA-seq values. Therefore, RNA4-scaled PGR status is
#       defined using a two-component Gaussian mixture model fitted in the
#       whole PANTHER RNA-seq cohort:
#
#           threshold = mean(higher-expression component) - 2 * SD
#
#   2. Continuous AURKA anchors:
#       Bossuyt et al. report QGP anchors:
#
#           PGR-high: 7.0 to 10.0 -> RNA4 risk 0 to 3
#           PGR-low:  7.0 to  8.5 -> RNA4 risk 0 to 3
#
#       In PANTHER RNA-seq, RNA4-scaled AURKA values span approximately
#       4 to 10 over the principal observed range. The RNA-seq adaptation
#       maps the PGR-high interval to 4.0 to 10.0 and preserves the published
#       halfway location of the PGR-low upper anchor:
#
#           PGR-high: 4.0 to 10.0 -> RNA4 risk 0 to 3
#           PGR-low:  4.0 to  7.0 -> RNA4 risk 0 to 3
#
# Important interpretation:
#   This is an exploratory RNA-seq-adapted pathologic SET2,3 calculation.
#   It is not a formally QGP-calibrated implementation because the published
#   RNA-seq-to-QGP calibration coefficients are unavailable.
#
# References:
#   Bossuyt V, et al. Clin Chem. 2021;67:1240-1248.
#   Lau R, et al. Clin Chem. 2020;66:934-945.
#   Du L, et al. Ann Oncol. 2021;32:642-651.
#   Speers CW, et al. J Clin Oncol. 2023;41:1841-1848.
# ==============================================================================

set23_signature_pathological <- function(gene_exp_data, clinical_data) {
  
  # --------------------------------------------------------------------------
  # 0. Package requirements
  # --------------------------------------------------------------------------
  
  if (!requireNamespace("mclust", quietly = TRUE)) { stop("Package 'mclust' is required. Please install it before running this function.") }
  
  if (!requireNamespace("tibble", quietly = TRUE)) { stop("Package 'tibble' is required. Please install it before running this function.") }
  
  # --------------------------------------------------------------------------
  # 1. Internal constants: genes and published/adapted score parameters
  # --------------------------------------------------------------------------
  
  erpr_genes <- c("ABAT", "ADCY1", "AZGP1", "CA12", "CD2", "CD3D", "DNAJC12", "ESR1", "KCNE4", "MAPT", "MRPS30", "NAT1", "NPY1R", "PDZK1", "QDPR", "SCUBE2", "SLC39A6", "STC2")
  ref_genes <- c("AK2", "APPBP2", "ATP5J2", "DARS", "LDHA", "TRIM2", "UBE2Z", "UGP2", "VDAC2", "WIPF2")
  rna4_genes <- c("ESR1", "PGR", "ERBB2", "AURKA")
  
  # Bossuyt et al. reported coefficients
  coefficient_BPI <- 0.5056
  coefficient_SET_ER_PR <- 0.7542
  
  # Bossuyt et al. reported QGP continuous AURKA anchors
  published_QGP_AURKA_anchors <- list(
    PGR_high_lower = 7.0,
    PGR_high_upper = 10.0,
    PGR_low_lower = 7.0,
    PGR_low_upper = 8.5)
  
  # Published relative position of the PGR-low upper anchor within the
  # PGR-high AURKA scaling interval:
  #
  #   (8.5 - 7.0) / (10.0 - 7.0) = 0.5
  pgr_low_relative_upper_position <- (published_QGP_AURKA_anchors$PGR_low_upper - published_QGP_AURKA_anchors$PGR_low_lower) / (published_QGP_AURKA_anchors$PGR_high_upper - published_QGP_AURKA_anchors$PGR_high_lower)
  
  # Fixed exploratory RNA-seq-adapted AURKA anchors
  RNAseq_AURKA_anchors <- list(
    PGR_high_lower = 4.0,
    PGR_high_upper = 10.0,
    PGR_low_lower = 4.0,
    PGR_low_upper = 4.0 + pgr_low_relative_upper_position * (10.0 - 4.0))
  
  stopifnot(isTRUE(all.equal(pgr_low_relative_upper_position, 0.5)), isTRUE(all.equal(RNAseq_AURKA_anchors$PGR_low_upper, 7.0)))
  
  # --------------------------------------------------------------------------
  # 2. Validate and retrieve whole-cohort reference objects
  #
  # These are required only to estimate the fixed RNA-seq PGR threshold.
  # --------------------------------------------------------------------------
  
  if (!exists("tmm_normalized_tpm_symbols", inherits = TRUE)) { stop("Whole-cohort RNA-seq object 'tmm_normalized_tpm_symbols' was not found in the environment.") }
  
  if (!exists("all_clinical", inherits = TRUE)) { stop("Whole-cohort clinical object 'all_clinical' was not found in the environment.") }
  
  whole_expr <- as.matrix(get("tmm_normalized_tpm_symbols", inherits = TRUE))
  
  whole_clin <- get("all_clinical", inherits = TRUE)
  
  if (is.null(rownames(whole_expr))) { stop("The whole-cohort RNA-seq matrix 'tmm_normalized_tpm_symbols' must have gene symbols as rownames.") }
  
  if (is.null(colnames(whole_expr))) { stop("The whole-cohort RNA-seq matrix 'tmm_normalized_tpm_symbols' must have sample IDs as column names.") }
  
  if (!"sample_id" %in% colnames(whole_clin)) { stop("'all_clinical' must contain a 'sample_id' column.") }
  
  if (anyDuplicated(whole_clin$sample_id) > 0) { stop("Duplicated sample_id values detected in 'all_clinical'.") }
  
  if (anyDuplicated(colnames(whole_expr)) > 0) { stop("Duplicated sample IDs detected in 'tmm_normalized_tpm_symbols'.") }
  
  required_whole_genes <- unique(c(ref_genes, "PGR"))
  missing_whole_genes <- setdiff(required_whole_genes, rownames(whole_expr))
  
  if (length(missing_whole_genes) > 0) { stop("Missing required genes in 'tmm_normalized_tpm_symbols': ", paste(missing_whole_genes, collapse = ", ")) }
  
  # --------------------------------------------------------------------------
  # 3. Validate and align the target cohort to be scored
  # --------------------------------------------------------------------------
  
  expr <- as.matrix(gene_exp_data)
  clin <- clinical_data
  
  if (is.null(rownames(expr))) { stop("gene_exp_data must have gene symbols as rownames.") }
  
  if (is.null(colnames(expr))) { stop("gene_exp_data must have sample IDs as column names.") }
  
  required_clin_cols <- c("sample_id", "tumdiam", "posnodes")
  missing_clin_cols <- setdiff(required_clin_cols, colnames(clin))
  
  if (length(missing_clin_cols) > 0) { stop("Missing required columns in clinical_data: ", paste(missing_clin_cols, collapse = ", ")) }
  
  if (anyDuplicated(clin$sample_id) > 0) { stop("Duplicated sample_id values detected in clinical_data.") }
  
  if (anyDuplicated(colnames(expr)) > 0) { stop("Duplicated sample IDs detected in gene_exp_data.") }
  
  required_target_genes <- unique(c(erpr_genes, ref_genes, rna4_genes))
  missing_target_genes <- setdiff(required_target_genes, rownames(expr))
  
  if (length(missing_target_genes) > 0) { stop("Missing required genes in gene_exp_data: ", paste(missing_target_genes, collapse = ", ")) }
  
  target_samples <- clin$sample_id[clin$sample_id %in% colnames(expr)]
  
  if (length(target_samples) == 0) { stop("No overlapping sample IDs between gene_exp_data and clinical_data.") }
  
  clin <- clin[match(target_samples, clin$sample_id), , drop = FALSE]
  expr <- expr[, match(target_samples, colnames(expr)), drop = FALSE]
  
  if (!identical(colnames(expr), clin$sample_id)) { stop("Target clinical and RNA-seq samples could not be aligned correctly.") }
  
  rownames(clin) <- clin$sample_id
  
  # For faithful calculation, all required expression measurements must be
  # available for every scored sample; do not silently average incomplete sets.
  if (anyNA(expr[required_target_genes, , drop = FALSE])) { stop("Missing expression values detected among required SET2,3 genes in the target cohort.") }
  
  # --------------------------------------------------------------------------
  # 4. Align the whole PANTHER cohort used for PGR threshold derivation
  # --------------------------------------------------------------------------
  
  whole_samples <- whole_clin$sample_id[whole_clin$sample_id %in% colnames(whole_expr)]
  
  if (length(whole_samples) == 0) { stop("No overlapping sample IDs between 'all_clinical' and 'tmm_normalized_tpm_symbols'.") }
  
  whole_expr <- whole_expr[ , match(whole_samples, colnames(whole_expr)), drop = FALSE]
  
  if (anyNA(whole_expr[required_whole_genes, , drop = FALSE])) { stop("Missing expression values detected among reference/PGR genes in the whole-cohort threshold derivation set.") }
  
  # --------------------------------------------------------------------------
  # 5. Calculate SET_ER_PR in the target cohort
  #
  # Bossuyt et al.:
  #   SET_ER_PR = mean expression of 18 target genes - mean expression of 10 reference genes
  #
  # No +2 or +10 constant is added to SET_ER_PR.
  # --------------------------------------------------------------------------
  
  target_ref_average <- colMeans(expr[ref_genes, , drop = FALSE])
  
  target_erpr_average <- colMeans(expr[erpr_genes, , drop = FALSE])
  
  clin$SET_ER_PR <- as.numeric(target_erpr_average - target_ref_average) # + 2.5
  
  # --------------------------------------------------------------------------
  # 6. Calculate RNA4-scaled expression in the target cohort
  #
  # Bossuyt et al.:
  #   RNA4-scaled gene expression = individual gene expression - mean(reference genes) + 10
  # --------------------------------------------------------------------------
  
  target_rna4_scaled <- sweep(expr[rna4_genes, , drop = FALSE], 2, target_ref_average, FUN = "-") + 10
  
  clin$RNA4_ESR1_scaled <- as.numeric(target_rna4_scaled["ESR1", ])
  clin$RNA4_PGR_scaled <- as.numeric(target_rna4_scaled["PGR", ])
  clin$RNA4_ERBB2_scaled <- as.numeric(target_rna4_scaled["ERBB2", ])
  clin$RNA4_AURKA_scaled <- as.numeric(target_rna4_scaled["AURKA", ])
  
  # --------------------------------------------------------------------------
  # 7. Derive the RNA-seq-adapted PGR threshold in the whole PANTHER cohort
  #
  # Bossuyt et al. used a fixed QGP-scale PGR threshold (>5.1).
  # Since uncalibrated RNA-seq values are not on the validated QGP scale,
  # this implementation derives an RNA-seq-specific threshold from the
  # higher PGR expression component in the whole PANTHER cohort.
  # --------------------------------------------------------------------------
  
  whole_ref_average <- colMeans(whole_expr[ref_genes, , drop = FALSE])
  
  whole_PGR_scaled <- as.numeric(whole_expr["PGR", ] - whole_ref_average + 10)
  
  define_high_peak_threshold_gmm <- function(x) {
    
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    
    if (length(x) < 20) { stop("Too few non-missing values to estimate the whole-cohort PGR threshold robustly.") }
    
    # Mclust() calls mclustBIC() by name in the calling environment, so make it visible here (mclust need not be attached)
    mclustBIC <- mclust::mclustBIC
    fit <- mclust::Mclust(data = x, G = 2, verbose = FALSE)
    
    if (is.null(fit) || is.null(fit$parameters$mean)) { stop("Two-component Gaussian mixture model fitting failed for PGR.") }
    
    component_means <- as.numeric(fit$parameters$mean)
    component_vars <- as.numeric(fit$parameters$variance$sigmasq)
    
    # mclust returns one common variance when the fitted two-component
    # model uses equal component variance.
    if (length(component_vars) == 1) {
      component_vars <- rep(component_vars, 2)
    }
    
    if (length(component_means) != 2 || length(component_vars) != 2 || any(!is.finite(component_means)) || any(!is.finite(component_vars))) { stop("Unexpected or invalid parameter structure from PGR GMM.") }
    
    high_component <- which.max(component_means)
    
    threshold <- component_means[high_component] - 2 * sqrt(component_vars[high_component])
    
    list(
      threshold = as.numeric(threshold),
      model = fit,
      component_means = component_means,
      component_sd = sqrt(component_vars),
      component_proportions = as.numeric(fit$parameters$pro),
      high_component = high_component)
  }
  
  whole_PGR_gmm <- define_high_peak_threshold_gmm(whole_PGR_scaled)
  PGR_threshold_whole_cohort <- whole_PGR_gmm$threshold
  
  # --------------------------------------------------------------------------
  # 8. Apply the frozen whole-cohort PGR threshold to the target cohort
  # --------------------------------------------------------------------------
  
  clin$RNA4_PGR_high <- ifelse(is.na(clin$RNA4_PGR_scaled), NA, clin$RNA4_PGR_scaled > PGR_threshold_whole_cohort)
  
  clin$RNA4_PGR_group <- factor(
    ifelse(is.na(clin$RNA4_PGR_high), NA,
      ifelse(clin$RNA4_PGR_high, "PGR-high", "PGR-low")), levels = c("PGR-low", "PGR-high"))
  
  pgr_high <- !is.na(clin$RNA4_PGR_high) & clin$RNA4_PGR_high
  pgr_low <- !is.na(clin$RNA4_PGR_high) & !clin$RNA4_PGR_high
  
  if (sum(pgr_low) == 0) { stop("No PGR-low samples were identified in the target cohort after application of the frozen whole-cohort RNA-seq PGR threshold.") }
  
  if (sum(pgr_high) == 0) { stop("No PGR-high samples were identified in the target cohort after application of the frozen whole-cohort RNA-seq PGR threshold.") }
  
  # --------------------------------------------------------------------------
  # 9. Calculate continuous RNA4 risk score from AURKA and PGR context
  #
  # Bossuyt et al. continuous QGP structure:
  #   PGR-high: 7.0 to 10.0 -> RNA4 risk 0 to 3
  #   PGR-low:  7.0 to  8.5 -> RNA4 risk 0 to 3
  #
  # RNA-seq adaptation used here:
  #   PGR-high: 4.0 to 10.0 -> RNA4 risk 0 to 3
  #   PGR-low:  4.0 to  7.0 -> RNA4 risk 0 to 3
  #
  # The PGR-low upper anchor remains halfway through the PGR-high
  # scaling interval, preserving the published PGR-dependent scoring
  # geometry while adapting the lower AURKA floor to RNA-seq.
  # --------------------------------------------------------------------------
  
  scale_to_0_3 <- function(x, low_anchor, high_anchor) {
    
    if (!is.finite(low_anchor) || !is.finite(high_anchor) || high_anchor <= low_anchor) { stop("Invalid AURKA scaling anchors.") }
    
    unbounded_score <- ((x - low_anchor) / (high_anchor - low_anchor)) * 3
    
    pmin(pmax(unbounded_score, 0), 3)
  }
  
  clin$RNA4_risk_score <- NA_real_
  
  clin$RNA4_risk_score[pgr_high] <- scale_to_0_3(
    x = clin$RNA4_AURKA_scaled[pgr_high],
    low_anchor = RNAseq_AURKA_anchors$PGR_high_lower,
    high_anchor = RNAseq_AURKA_anchors$PGR_high_upper)
  
  clin$RNA4_risk_score[pgr_low] <- scale_to_0_3(
    x = clin$RNA4_AURKA_scaled[pgr_low],
    low_anchor = RNAseq_AURKA_anchors$PGR_low_lower,
    high_anchor = RNAseq_AURKA_anchors$PGR_low_upper)
  
  # --------------------------------------------------------------------------
  # 10. Calculate pathologic pT risk score
  #
  # Bossuyt et al.:
  #   - score = 0 for tumour diameter <1.0 cm
  #   - score = 3 for tumour diameter >3.5 cm
  #   - otherwise linearly scaled from 0 to 3 in increments of 0.1
  #
  # clinical_data$tumdiam is assumed to be in millimetres.
  # --------------------------------------------------------------------------
  
  tumor_size_mm <- suppressWarnings(as.numeric(clin$tumdiam))
  tumor_size_cm <- tumor_size_mm / 10
  
  if (any(!is.na(tumor_size_mm) & tumor_size_mm < 0)) { stop("Negative values detected in tumdiam.") }
  
  clin$pT_risk_score <- ifelse(is.na(tumor_size_cm), NA_real_,
                               ifelse(tumor_size_cm < 1.0, 0,
                                      ifelse(tumor_size_cm > 3.5, 3,
                                             round(((tumor_size_cm - 1.0) / (3.5 - 1.0)) * 3, digits = 1))))
  
  # --------------------------------------------------------------------------
  # 11. Calculate pathologic pN risk score
  #
  # Bossuyt et al.:
  #   - score = 5 for >10 involved nodes
  #   - otherwise linearly scaled from 0 to 5 in increments of 0.5
  #     for 0 to 10 positive nodes
  # --------------------------------------------------------------------------
  
  positive_nodes <- suppressWarnings(as.numeric(clin$posnodes))
  
  if (any(!is.na(positive_nodes) & positive_nodes < 0)) { stop("Negative values detected in posnodes.") }
  
  if (any(!is.na(positive_nodes) & positive_nodes %% 1 != 0)) { warning("Non-integer values detected in posnodes. Verify that posnodes represents the number of involved lymph nodes.") }
  
  clin$pN_risk_score <- ifelse(is.na(positive_nodes), NA_real_,
    ifelse(positive_nodes > 10, 5, round(positive_nodes * 0.5, digits = 1)))
  
  # --------------------------------------------------------------------------
  # 12. Calculate pathologic BPI
  #
  # Bossuyt et al.:
  #   BPI = [11 - (pT risk + pN risk + RNA4 risk)] * 4 / 11
  #
  # Higher BPI indicates more indolent baseline prognosis.
  # --------------------------------------------------------------------------
  
  clin$BPI_RNAseq_adapted <- (11 - (clin$pT_risk_score + clin$pN_risk_score + clin$RNA4_risk_score)) * 4 / 11
  
  # Alias retained for compatibility with downstream code.
  clin$BPI <- clin$BPI_RNAseq_adapted
  
  # --------------------------------------------------------------------------
  # 13. Calculate SET2,3
  #
  # Bossuyt et al.:
  #   SET2,3 = 0.5056 * BPI + 0.7542 * SET_ER_PR
  #
  # The name SET_2_3_RNAseq_adapted is preferred for reporting because
  # the RNA-seq PGR threshold and AURKA anchors are adapted rather than
  # formally calibrated to QGP.
  # --------------------------------------------------------------------------
  
  clin$SET_2_3_RNAseq_adapted <- (coefficient_BPI * clin$BPI_RNAseq_adapted + coefficient_SET_ER_PR * clin$SET_ER_PR)
  
  # Alias retained for compatibility with existing analysis scripts.
  clin$SET_2_3 <- clin$SET_2_3_RNAseq_adapted
  
  # Optional exploratory postsurgical high/low category.
  # Speers et al. used >2.10 as high for pathologic SET2,3 after formal calibration. 
  # Since the present score is RNA-seq-adapted but not formally QGP-calibrated, this category should be interpreted cautiously.
  clin$SET_2_3_high_2.10_exploratory <- ifelse(is.na(clin$SET_2_3_RNAseq_adapted), NA,
                                               clin$SET_2_3_RNAseq_adapted > 2.10)
  
  # --------------------------------------------------------------------------
  # 14. Supporting metadata for reproducibility
  # --------------------------------------------------------------------------
  
  RNA4_floor_ceiling_summary <- data.frame(
    PGR_group = c("PGR-low", "PGR-high"),
    N = c(sum(pgr_low), sum(pgr_high)),
    RNA4_score_zero_N = c(sum(clin$RNA4_risk_score[pgr_low] == 0, na.rm = TRUE), sum(clin$RNA4_risk_score[pgr_high] == 0, na.rm = TRUE)),
    RNA4_score_three_N = c(sum(clin$RNA4_risk_score[pgr_low] == 3, na.rm = TRUE), sum(clin$RNA4_risk_score[pgr_high] == 3, na.rm = TRUE)))
  
  RNA4_floor_ceiling_summary$RNA4_score_zero_percentage <- round(RNA4_floor_ceiling_summary$RNA4_score_zero_N / RNA4_floor_ceiling_summary$N * 100, 1)
  
  RNA4_floor_ceiling_summary$RNA4_score_three_percentage <- round(RNA4_floor_ceiling_summary$RNA4_score_three_N / RNA4_floor_ceiling_summary$N * 100, 1)
  
  meta <- list(
    implementation = list(
      name = "RNA-seq-adapted pathologic SET2,3",
      interpretation = paste(
        "Exploratory RNA-seq adaptation of the continuous pathologic",
        "SET2,3 calculation. The SETER/PR, RNA4 scaling, pT/pN scoring,",
        "BPI formula, and SET2,3 weighting follow Bossuyt et al. 2021.",
        "The PGR threshold and AURKA anchors are RNA-seq adaptations",
        "because formal RNA-seq-to-QGP calibration coefficients are unavailable.")),
    
    target_cohort = list(
      used_samples = target_samples,
      N = length(target_samples),
      PGR_low_N = sum(pgr_low),
      PGR_high_N = sum(pgr_high)),
    
    genes = list(
      SET_ER_PR_target_genes = erpr_genes,
      reference_genes = ref_genes,
      RNA4_genes = rna4_genes),
    
    Bossuyt_reported_calculation = list(
      SET_ER_PR = "mean(18 target genes) - mean(10 reference genes)",
      RNA4_scaling = "gene expression - mean(10 reference genes) + 10",
      BPI = "(11 - (pT_risk_score + pN_risk_score + RNA4_risk_score)) * 4 / 11",
      SET_2_3 = "0.5056 * BPI + 0.7542 * SET_ER_PR",
      QGP_AURKA_anchors = published_QGP_AURKA_anchors),
    
    PGR_threshold_derivation = list(
      status = "RNA-seq-adapted",
      cohort = "Whole PANTHER cohort",
      expression_object = "tmm_normalized_tpm_symbols",
      clinical_object = "all_clinical",
      N = length(whole_samples),
      method = paste("Two-component Gaussian mixture model on RNA4-scaled PGR; threshold = mean(higher-expression component) - 2 SD"),
      threshold = PGR_threshold_whole_cohort,
      selected_model = whole_PGR_gmm$model$modelName,
      component_means = whole_PGR_gmm$component_means,
      component_sd = whole_PGR_gmm$component_sd,
      component_proportions = whole_PGR_gmm$component_proportions),
    
    AURKA_scaling = list(
      status = "RNA-seq-adapted",
      rationale = paste("Published QGP AURKA scoring geometry was preserved after mapping the approximate lower RNA-seq AURKA floor from 7.0 to 4.0."),
      published_QGP_anchors = published_QGP_AURKA_anchors,
      published_PGR_low_relative_upper_position = pgr_low_relative_upper_position,
      RNAseq_adapted_anchors = RNAseq_AURKA_anchors,
      floor_ceiling_summary = RNA4_floor_ceiling_summary),
    
    pathologic_scoring = list(
      tumor_diameter_input_units = "millimetres",
      pT_rule = paste("0 if <1.0 cm; 3 if >3.5 cm; otherwise linear 0-3 scaling rounded to 0.1"),
      pN_rule = paste("5 if >10 involved nodes; otherwise 0.5 risk units per positive node")),
    
    caveat = paste(
      "Speers et al. applied SET2,3 to whole-transcriptome RNA-seq after calibration to QGP in an independent matched-sample cohort.",
      "The required transformation coefficients were not reported; therefore this implementation should not be described as formally QGP-calibrated."))
  
  return(list(clinical = tibble::as_tibble(clin), meta = meta))
}
