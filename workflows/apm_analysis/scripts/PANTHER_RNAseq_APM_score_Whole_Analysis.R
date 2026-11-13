###### RNA-Seq Data Analysis: APM score (PANTHER, SCAN-B, MATADOR, TEX) ######

# Run top to bottom in one session. Layout (follows the manuscript):
#   Packages, settings and constants
#   Functions (every reusable function of the script, grouped by use)
#   PANTHER: analysis inputs (clinical table, expression matrix); SET2,3; APM gene sets (BC360, APM80); APM signature (definition, scoring, APM_score / APM_group);
#            Table 1 and Supplementary Table 2; continuous APM score; APM groups; intrinsic subtypes and 21-gene RS; 14-gene IgG signature and immune deconvolution
#   Validation cohorts: SCAN-B, MATADOR, TEX
#   Median follow-up; main and supplementary figures; Supplementary Table 5; session info

###### Load required packages ######

# Attach packages for gene-set scoring, survival analysis, tables and figures.
packages = c("circlize", "ComplexHeatmap", "GSVA", "MCPcounter", "broom", "dplyr", "forcats", "forestploter", "genefu", "ggplot2", "ggpubr", 
             "grid", "gridExtra", "gt", "gtsummary", "lmtest", "readxl", "survival", "survminer", "tibble", "gtable", "patchwork", "GEOquery")

lapply(packages, library, character.only = TRUE)

# Installed but not attached (used as pkg::fun only): Biobase, curl, effsize, ggtext, mclust, openxlsx, preprocessCore, purrr, rlang, rstatix, splines, tidyr, here; 
# webshot2 (needed by gt::gtsave() for PNG files); future and future.apply (permutation script)

###### Settings: folders, Cox adjustment and plot colors ######

# Folders (see data/README.md):
#   data_dir     input data, in data/ under the repository root
#   code_dir     the repository root, found by here::i_am(); results are written to code_dir/results/APM
#   chrome_path  optional: path to the Chrome executable, only if gt::gtsave() cannot find Chrome when saving tables as PNG
here::i_am("workflows/apm_analysis/scripts/PANTHER_RNAseq_APM_score_Whole_Analysis.R")   # finds the repository root
code_dir = here::here()
data_dir = here::here("data")
chrome_path = ""

output_dir = file.path(code_dir, "results", "APM")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

adjust_vars = c("ERpos", "newHER2", "tumdiam", "posnodes", "histgrade")
adjust_vars_no_clin_sub = c("tumdiam", "posnodes", "histgrade")
adjust_vars_no_her2_sub = c("ERpos", "tumdiam", "posnodes", "histgrade")
adjust_vars_no_er_sub = c("newHER2", "tumdiam", "posnodes", "histgrade")

adjust_terms = paste("strata(ERpos)", "newHER2", "scale(tumdiam)", "scale(posnodes)", "histgrade", sep = " + ")
adjust_terms_no_clin_sub = paste("scale(tumdiam)", "scale(posnodes)", "histgrade", sep = " + ")
adjust_terms_no_her2_sub = paste("strata(ERpos)", "scale(tumdiam)", "scale(posnodes)", "histgrade", sep = " + ")
adjust_terms_no_er_sub = paste("newHER2", "scale(tumdiam)", "scale(posnodes)", "histgrade", sep = " + ")

group_col_treat = c('FEC-D' = 'dodgerblue3' , 'tddEC-D' = 'maroon3')
group_col_apm = c('APM-low' = '#3B6EA8', 'APM-high' = 'seagreen4')

###### ----- Functions & Helpers -----
###### External: SET2,3 index (set23_signature_pathological) ######

# Load external functions for pathological SET2,3 calculation.
source(file.path(code_dir, "workflows", "pathologic_set", "scripts", "Compute_SET_2_3_Signature_pathological.R"))

###### Formatting and table helpers ######

# Format p-values while preserving precision for small values.
format_p_regular = function(p, digits = 3) {
  if (is.null(p) || length(p) == 0) return(NA_character_)
  
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.0001 ~ "<0.0001",
    p < 0.001 ~ formatC(p, format = "f", digits = 4),   # e.g. 0.0003 instead of 0.000
    TRUE ~ formatC(p, format = "f", digits = digits))
}

# "HR (low–high)" with two decimals and an en dash, as in the KM plots
format_hr_ci = function(est, low, high, digits = 2) {
  f = function(x) formatC(x, format = "f", digits = digits)
  paste0(f(est), " (", f(low), "\u2013", f(high), ") ")
}

# Export a UTF-8 CSV file, creating its parent directory if needed.
write_table = function(x, path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  utils::write.csv(x,path,row.names=FALSE,na="",fileEncoding="UTF-8")
}

###### Continuous treatment-effect plots (HR across a continuous score) ######

# CORE PLOT: one biomarker + one endpoint + one dataset/subgroup
plot_hr_vs_continuous = function(data, time_col, status_col, continuous_var, continuous_label = NULL, treatment_col = "TREAT", treat_levels = c("FEC-D", "tddEC-D"),
                                 adjust_vars = NULL, adjust_terms = NULL, df_spline = 2, n_grid = 200, endpoint_label = NULL, subtype_label = NULL, y_limits = c(-0.10, 4)) {
  # symbols
  time_sym   = rlang::sym(time_col)
  status_sym = rlang::sym(status_col)
  cont_sym   = rlang::sym(continuous_var)
  treat_sym  = rlang::sym(treatment_col)

  # required columns
  required_cols = unique(c(time_col, status_col, continuous_var, treatment_col, adjust_vars))

  # filter + prepare
  d = data %>%
    dplyr::filter(.data[[treatment_col]] %in% treat_levels) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.x))) %>%
    dplyr::mutate(time_years = .data[[time_col]] / 365.25, "{treatment_col}" := factor(.data[[treatment_col]], levels = treat_levels))

  if (nrow(d) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = "No data available"))
  }

  # precompute scaled covariates if used
  # This avoids problems with scale() inside model.matrix(newdata)
  if ("tumdiam" %in% adjust_vars) {
    tumdiam_scaled = scale(d$tumdiam)
    d$tumdiam_z = as.numeric(tumdiam_scaled)
  }

  if ("posnodes" %in% adjust_vars) {
    posnodes_scaled = scale(d$posnodes)
    d$posnodes_z = as.numeric(posnodes_scaled)
  }

  # labels
  if (is.null(subtype_label)) { subtype_label = "All" }
  if (is.null(endpoint_label)) { endpoint_label = gsub("^time", "", time_col) }
  if (is.null(continuous_label)) { continuous_label = rlang::as_name(cont_sym) }

  # adjustment terms
  adjust_term = if (!is.null(adjust_terms) && nzchar(adjust_terms)) {
    paste0(" + ", adjust_terms)
  } else if (length(adjust_vars) > 0) {
    paste0(" + ", paste(adjust_vars, collapse = " + "))
  } else {
    ""
  }

  # replace inline scale() with precomputed z variables
  adjust_term_model = adjust_term %>%
    gsub("scale\\(tumdiam\\)", "tumdiam_z", .) %>%
    gsub("scale\\(posnodes\\)", "posnodes_z", .)

  # fit spline interaction Cox model
  fml = stats::as.formula(paste0("survival::Surv(time_years, ", status_col, ") ~ ", treatment_col, " * splines::ns(", continuous_var, ", df = ", df_spline, ")", adjust_term_model))

  cox_model = survival::coxph(fml, data = d, control = survival::coxph.control(iter.max = 50))
  
  # Test departure from a linear treatment-by-score interaction.
  d$trt_x_score = as.numeric(d[[treatment_col]] == treat_levels[2]) * d[[continuous_var]]
  fml_lin_int = stats::as.formula(paste0("survival::Surv(time_years, ", status_col, ") ~ ", treatment_col, " + splines::ns(", continuous_var, ", df = ", df_spline, ") + trt_x_score", adjust_term_model))
  cox_lin_int = survival::coxph(fml_lin_int, data = d, control = survival::coxph.control(iter.max = 50))
  ll_full = stats::logLik(cox_model); ll_red = stats::logLik(cox_lin_int)
  p_nonlinear = stats::pchisq(2 * (as.numeric(ll_full) - as.numeric(ll_red)), df = attr(ll_full, "df") - attr(ll_red, "df"), lower.tail = FALSE)
  message(continuous_label, " / ", endpoint_label, ": P non-linear interaction = ", format_p_regular(p_nonlinear))

  # interaction p-values
  interaction_stats = get_linear_interaction_pvalues(
    data = d,
    time_col = time_col,
    status_col = status_col,
    continuous_var = continuous_var,
    treatment_col = treatment_col,
    treat_levels = treat_levels,
    adjust_vars = adjust_vars,
    adjust_terms = adjust_terms)

  p_label_unadjusted = paste0("P-inter., linear = ", format_p_regular(interaction_stats$p_unadjusted))
  p_label_adjusted = paste0("P-inter., linear adj. = ", format_p_regular(interaction_stats$p_adjusted))

  # prediction grid
  range_vals = range(d[[rlang::as_name(cont_sym)]], na.rm = TRUE)
  pred_range = seq(range_vals[1], range_vals[2], length.out = n_grid)

  newdata = expand.grid(
    cont_val = pred_range,
    treat_val = treat_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)

  colnames(newdata) = c(
    rlang::as_name(cont_sym),
    rlang::as_name(treat_sym))

  newdata[[rlang::as_name(treat_sym)]] = factor(
    newdata[[rlang::as_name(treat_sym)]],
    levels = treat_levels)

  # Add reference values for adjustment variables
  if (length(adjust_vars) > 0) {
    for (v in adjust_vars) {

      # tumdiam and posnodes enter the model as z variables
      if (v %in% c("tumdiam", "posnodes")) {
        next
      }

      if (is.factor(d[[v]])) {
        reference_value = names(which.max(table(d[[v]])))
        newdata[[v]] = factor(reference_value, levels = levels(d[[v]]))

      } else if (is.character(d[[v]])) {
        reference_value = names(which.max(table(d[[v]])))
        newdata[[v]] = factor(reference_value, levels = sort(unique(d[[v]])))

      } else {
        newdata[[v]] = median(d[[v]], na.rm = TRUE)
      }
    }
  }

  # reference values for scaled continuous covariates
  # 0 corresponds to the mean of the covariate in the fitted dataset
  if ("tumdiam" %in% adjust_vars) { newdata$tumdiam_z = 0 }

  if ("posnodes" %in% adjust_vars) { newdata$posnodes_z = 0 }

  beta = stats::coef(cox_model)
  V    = stats::vcov(cox_model)

  if (anyNA(beta)) {
    stop("The spline Cox model has non-estimable coefficients.")
  }

  ## Reuse the fitted Cox model's spline basis and factor contrasts
  X = stats::model.matrix(cox_model, data = newdata)

  stopifnot(all(names(beta) %in% colnames(X)))

  X = X[, names(beta), drop = FALSE]

  # remove intercept if present
  if ("(Intercept)" %in% colnames(X)) {
    X = X[, colnames(X) != "(Intercept)", drop = FALSE]
  }

  # align prediction matrix to Cox coefficients
  X = X[, names(beta), drop = FALSE]

  idx1 = which(newdata[[rlang::as_name(treat_sym)]] == treat_levels[1])
  idx2 = which(newdata[[rlang::as_name(treat_sym)]] == treat_levels[2])

  stopifnot(length(idx1) == length(idx2))

  X1 = X[idx1, , drop = FALSE]
  X2 = X[idx2, , drop = FALSE]

  C = X2 - X1

  diff_lp = as.numeric(C %*% beta)
  se_diff = sqrt(rowSums((C %*% V) * C))

  hr_df = data.frame(
    x_var    = pred_range,
    diff_lp  = diff_lp,
    se_diff  = se_diff,
    HR       = exp(diff_lp),
    HR_lower = exp(diff_lp - 1.96 * se_diff),
    HR_upper = exp(diff_lp + 1.96 * se_diff))

  hr_axis_label = if (length(adjust_vars) > 0) {
    "Hazard ratio (HR), adjusted"
  } else {
    "Hazard ratio (HR)"
  }

  # hide HR and CI lines outside plotting area
  hr_df = hr_df %>% dplyr::mutate(HR_plot = dplyr::if_else(HR >= y_limits[1] & HR <= y_limits[2], HR, NA_real_),
                                  HR_lower_plot = dplyr::if_else(HR_lower >= y_limits[1] & HR_lower <= y_limits[2], HR_lower, NA_real_),
                                  HR_upper_plot = dplyr::if_else(HR_upper >= y_limits[1] & HR_upper <= y_limits[2], HR_upper, NA_real_))

  ggplot2::ggplot(hr_df, ggplot2::aes(x = x_var)) +
    ggplot2::geom_line(ggplot2::aes(y = HR_plot), color = "black", linewidth = 1.2, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = HR_lower_plot), linetype = "dashed", linewidth = 0.5, color = "gray50", na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = HR_upper_plot), linetype = "dashed", linewidth = 0.5, color = "gray50", na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 1, linetype = "solid", color = "gray") +
    ggplot2::geom_segment(data = d %>% dplyr::filter(!!treat_sym == treat_levels[1]), ggplot2::aes(x = !!cont_sym, xend = !!cont_sym, y = -0.05, yend = 0.00), inherit.aes = FALSE, color = "dodgerblue3", alpha = 0.6) +
    ggplot2::geom_segment(data = d %>% dplyr::filter(!!treat_sym == treat_levels[2]), ggplot2::aes(x = !!cont_sym, xend = !!cont_sym, y = -0.10, yend = -0.05), inherit.aes = FALSE, color = "maroon3", alpha = 0.6) +
    ggplot2::annotate("segment", x = range_vals[2] - 0.65, xend = range_vals[2] - 0.55, y = 0.8, yend = 0.8, color = "dodgerblue3", linewidth = 2) +
    ggplot2::annotate("text", x = range_vals[2] - 0.5, y = 0.8, label = treat_levels[1], hjust = 0, size = 4) +
    ggplot2::annotate("segment", x = range_vals[2] - 0.65, xend = range_vals[2] - 0.55, y = 0.6, yend = 0.6, color = "maroon3", linewidth = 2) +
    ggplot2::annotate("text", x = range_vals[2] - 0.5, y = 0.6, label = treat_levels[2], hjust = 0, size = 4) +
    # P-int values
    ggplot2::annotate("text", x = range_vals[2]-2.2, y = y_limits[2] - 0.15, label = p_label_unadjusted, hjust = 0, vjust = 1, size = 5) +
    ggplot2::annotate("text", x = range_vals[2]-2.2, y = y_limits[2] - 0.45, label = p_label_adjusted, hjust = 0, vjust = 1, size = 5) +
    ggplot2::annotate("text", x = range_vals[2]-2.2, y = y_limits[2] - 0.75, label = paste0("P-inter., non-linearity = ", format_p_regular(p_nonlinear)), hjust = 0, vjust = 1, size = 5) +
    ggplot2::coord_cartesian(xlim = range_vals, ylim = y_limits) +
    ggplot2::labs(x = continuous_label, y = hr_axis_label, title = paste0("tddEC-D vs. FEC-D across continuous ", continuous_label, " for ", endpoint_label)) +
    ggplot2::theme_bw(base_size = 14)+
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5))
}

# SAFE WRAPPER: one biomarker + one endpoint
make_plot_hr_single_endpoint = function(data, biomarker, biomarker_label = NULL, endpoint, subtype_label = "All", adjust_vars = NULL, adjust_terms = NULL, y_limits = c(-0.10, 4), min_events = 5) {

  time_col   = paste0("time", endpoint)
  status_col = paste0(endpoint, ".status")

  if (!all(c(time_col, status_col, biomarker) %in% colnames(data))) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(subtype_label, ": missing columns for ", endpoint)))
  }

  status_vec = data[[status_col]]

  if (all(is.na(status_vec)) || sum(status_vec == 1, na.rm = TRUE) < min_events) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::labs(title = paste0(subtype_label, ": ", endpoint, " (too few events)")))
  }

  p = plot_hr_vs_continuous(
    data = data,
    time_col = time_col,
    status_col = status_col,
    continuous_var = biomarker,
    continuous_label = biomarker_label,
    endpoint_label = endpoint,
    subtype_label = subtype_label,
    adjust_vars = adjust_vars,
    adjust_terms = adjust_terms,
    y_limits = y_limits)

  p
}

###### Interaction tests (continuous and median split) ######

# Test continuous treatment-by-score interactions using unadjusted and adjusted Cox LRTs in the same complete cases.
get_linear_interaction_pvalues = function(data, time_col, status_col, continuous_var, treatment_col = "TREAT", treat_levels = c("FEC-D", "tddEC-D"), adjust_vars = NULL, adjust_terms = NULL) {

  required_cols = unique(c(
    time_col,
    status_col,
    continuous_var,
    treatment_col,
    adjust_vars))

  missing_cols = setdiff(required_cols, colnames(data))

  if (length(missing_cols) > 0) { stop("Missing variables in data: ", paste(missing_cols, collapse = ", ")) }

  # Use the same complete-case analysis population for both
  # unadjusted and adjusted interaction tests
  d = data %>%
    dplyr::filter(.data[[treatment_col]] %in% treat_levels) %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::all_of(required_cols), ~ !is.na(.x))) %>%
    dplyr::mutate(time_years = .data[[time_col]] / 365.25, "{treatment_col}" := factor(.data[[treatment_col]], levels = treat_levels))

  if (nrow(d) == 0) { stop("No complete-case observations available for interaction models.") }

  surv_term = paste0("survival::Surv(time_years, ", status_col, ")")

  adjust_term = if (!is.null(adjust_terms) && nzchar(adjust_terms)) {
    paste0(" + ", adjust_terms)
  } else if (length(adjust_vars) > 0) {
    paste0(" + ", paste(adjust_vars, collapse = " + "))
  } else {
    ""
  }

  # Helper: likelihood-ratio test from two nested Cox models
  get_lrt_p = function(reduced_model, full_model) {

    ll_reduced = stats::logLik(reduced_model)
    ll_full = stats::logLik(full_model)

    chisq = 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))

    df_test = attr(ll_full, "df") - attr(ll_reduced, "df")

    stats::pchisq(q = chisq, df = df_test, lower.tail = FALSE)
  }

  # --------------------------------
  # Unadjusted linear interaction LRT
  # --------------------------------

  fml_unadj_reduced = stats::as.formula(paste0(surv_term, " ~ ", treatment_col, " + ", continuous_var))

  fml_unadj_full = stats::as.formula(paste0(surv_term, " ~ ", treatment_col, " * ", continuous_var))

  model_unadj_reduced = survival::coxph(fml_unadj_reduced, data = d )

  model_unadj_full = survival::coxph(fml_unadj_full, data = d)

  p_unadjusted = get_lrt_p(reduced_model = model_unadj_reduced, full_model = model_unadj_full)

  # ------------------------------
  # Adjusted linear interaction LRT
  # ------------------------------

  if (length(adjust_vars) > 0) {

    fml_adj_reduced = stats::as.formula(paste0(surv_term, " ~ ", treatment_col, " + ", continuous_var, adjust_term))

    fml_adj_full = stats::as.formula(paste0(surv_term, " ~ ", treatment_col, " * ", continuous_var, adjust_term))

    model_adj_reduced = survival::coxph(fml_adj_reduced, data = d)

    model_adj_full = survival::coxph(fml_adj_full, data = d)

    p_adjusted = get_lrt_p(reduced_model = model_adj_reduced, full_model = model_adj_full)

  } else {

    model_adj_reduced = NULL
    model_adj_full = NULL
    p_adjusted = NA_real_
  }

  list(
    p_unadjusted = as.numeric(p_unadjusted),
    p_adjusted = as.numeric(p_adjusted),
    n = nrow(d),
    events = sum(d[[status_col]] == 1, na.rm = TRUE),
    model_unadj_reduced = model_unadj_reduced,
    model_unadj_full = model_unadj_full,
    model_adj_reduced = model_adj_reduced,
    model_adj_full = model_adj_full)
}

# Split a biomarker at its supplied-population median; estimate treatment HRs within groups and test the adjusted interaction.
get_median_cutoff_pint = function(data, biomarker, endpoint, adjust_vars, adjust_terms) {

  sig_cat = paste0(biomarker, "_median_cat")

  time_col   = paste0("time", endpoint)
  status_col = paste0(endpoint, ".status")

  # Cohort median over every scored patient, fixed before any endpoint or covariate
  # exclusion (the same rule as APM_group: High >= median)
  cutoff_value = median(data[[biomarker]], na.rm = TRUE)

  required_cols = unique(c(biomarker, "TREAT", adjust_vars, time_col, status_col))

  d = data %>% dplyr::filter(dplyr::if_all(dplyr::all_of(required_cols), ~ !is.na(.x))) %>%
    dplyr::mutate("{sig_cat}" := factor(ifelse(.data[[biomarker]] >= cutoff_value, "High", "Low"), levels = c("Low", "High")),
                  TREAT = factor(TREAT, levels = c("FEC-D", "tddEC-D")))

  # Interaction models
  f_base = stats::as.formula(paste0("survival::Surv(", time_col, "/365.25, ", status_col, ") ~ ", "TREAT + ", sig_cat, " + ", adjust_terms))
  f_int = stats::as.formula(paste0("survival::Surv(", time_col, "/365.25, ", status_col, ") ~ ", "TREAT * ", sig_cat, " + ", adjust_terms))

  res.cox_base = survival::coxph(f_base, data = d)
  res.cox_int  = survival::coxph(f_int, data = d)

  mv_interaction_lr = lmtest::lrtest(res.cox_base, res.cox_int)
  mv_interaction_p_value = mv_interaction_lr$`Pr(>Chisq)`[2]

  # Helper: treatment HR within each median group
  get_group_hr = function(group_name) {

    d_group = d %>% dplyr::filter(.data[[sig_cat]] == group_name)

    f_group = stats::as.formula(paste0("survival::Surv(", time_col, "/365.25, ", status_col, ") ~ ", "TREAT + ", adjust_terms))

    fit = survival::coxph(f_group, data = d_group)

    res = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      dplyr::filter(term == "TREATtddEC-D")

    tibble::tibble(
      group = group_name,
      n = nrow(d_group),
      n_FEC_D = sum(d_group$TREAT == "FEC-D", na.rm = TRUE),
      n_tddEC_D = sum(d_group$TREAT == "tddEC-D", na.rm = TRUE),
      events = sum(d_group[[status_col]] == 1, na.rm = TRUE),
      HR = res$estimate,
      HR_low = res$conf.low,
      HR_high = res$conf.high,
      p_value = res$p.value)
  }

  group_results = dplyr::bind_rows(get_group_hr("Low"), get_group_hr("High"))

  group_results %>%
    dplyr::mutate(
      biomarker = biomarker,
      endpoint = endpoint,
      cutoff = cutoff_value,
      N_total = nrow(d),
      p_interaction = mv_interaction_p_value,
      HR_CI = paste0(round(HR, 2), " (", round(HR_low, 2), "-", round(HR_high, 2), ")"),
      p_value_print = format_p_regular(p_value),
      p_interaction_print = format_p_regular(p_interaction)) %>%
    dplyr::select(
      biomarker,
      endpoint,
      cutoff,
      N_total,
      group,
      n,
      n_FEC_D,
      n_tddEC_D,
      events,
      HR,
      HR_low,
      HR_high,
      HR_CI,
      p_value,
      p_value_print,
      p_interaction,
      p_interaction_print)
}

# Per-arm rows for one patient group: header, then the APM-score HR from separate adjusted Cox models; p interaction (pooled-model LRT, complete cases)
apm_score_by_arm_rows = function(data, endpoint, adjust_terms, adjust_vars, label, label_col, score = "APM_score", spacer = hr_column) {
  need = unique(c("TREAT", score, adjust_vars, paste0("time", endpoint), paste0(endpoint, ".status")))
  d = droplevels(data[complete.cases(data[, need, drop = FALSE]) & data$TREAT %in% c("FEC-D", "tddEC-D"), , drop = FALSE])
  surv = paste0("Surv(time", endpoint, "/365.25, ", endpoint, ".status) ~ ")
  
  base = coxph(as.formula(paste0(surv, "TREAT + ", score, " + ", adjust_terms)), data = d, na.action = na.fail)
  full = coxph(as.formula(paste0(surv, "TREAT * ", score, " + ", adjust_terms)), data = d, na.action = na.fail)
  p_int = lmtest::lrtest(base, full)$`Pr(>Chisq)`[2L]
  
  rows = data.frame(label = label, N = '', est = NA, low = NA, hi = NA, se = NA, ` ` = spacer,
                    `HR (95% CI)` = '', `p value` = '', `p interaction` = '', check.names = FALSE)
  for (arm in c("FEC-D", "tddEC-D")) {
    fit = coxph(as.formula(paste0(surv, score, " + ", adjust_terms)), data = d[d$TREAT == arm, , drop = FALSE])
    r = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    r = r[r$term == score, , drop = FALSE]
    stopifnot(nrow(r) == 1L)
    rows = rbind(rows, data.frame(label = paste0(arm, " arm"), N = fit$n, est = r$estimate, low = r$conf.low, hi = r$conf.high, se = r$std.error, ` ` = spacer,
                                  `HR (95% CI)` = format_hr_ci(r$estimate, r$conf.low, r$conf.high),
                                  `p value` = format_p_regular(r$p.value, digits = 3),
                                  `p interaction` = if (arm == "FEC-D") format_p_regular(p_int, digits = 3) else '',
                                  check.names = FALSE))
  }
  names(rows)[1] = label_col
  rows
}

# Three-way interaction (heterogeneity) LRT: does the TREAT x score interaction differ across the levels of `group`?
# Reduced: TREAT * score + TREAT:g + score:g (+ g if group_main) + adjust_terms; full adds TREAT:score:g.
# g = indicators of the non-reference levels of `group` (first level = reference); same complete cases in both models.
# The products are built as numeric columns, so R's factor coding cannot alias them with strata(ERpos).
three_way_interaction_lrt = function(data, endpoint, group, adjust_terms, adjust_vars, group_main = TRUE, score = "APM_score") {
  need = unique(c("TREAT", score, group, adjust_vars, paste0("time", endpoint), paste0(endpoint, ".status")))
  d = droplevels(data[complete.cases(data[, need, drop = FALSE]) & data$TREAT %in% c("FEC-D", "tddEC-D"), , drop = FALSE])
  
  trt = as.numeric(d$TREAT == "tddEC-D")
  G = stats::model.matrix(~ d[[group]])[, -1, drop = FALSE]
  k = seq_len(ncol(G))
  for (j in k) {
    d[[paste0(".g", j)]]   = G[, j]                         # subgroup level j
    d[[paste0(".tg", j)]]  = trt * G[, j]                   # TREAT x subgroup
    d[[paste0(".sg", j)]]  = d[[score]] * G[, j]            # score x subgroup
    d[[paste0(".tsg", j)]] = trt * d[[score]] * G[, j]      # TREAT x score x subgroup
  }
  
  lower = c(if (group_main) paste0(".g", k), paste0(".tg", k), paste0(".sg", k))
  f0 = as.formula(paste0("Surv(time", endpoint, "/365.25, ", endpoint, ".status) ~ TREAT * ", score, " + ",
                         paste(lower, collapse = " + "), " + ", adjust_terms))
  f1 = update(f0, as.formula(paste(". ~ . +", paste0(".tsg", k, collapse = " + "))))
  
  base = coxph(f0, data = d, na.action = na.fail)
  full = coxph(f1, data = d, na.action = na.fail)
  stopifnot(!anyNA(coef(full)))
  lr = lmtest::lrtest(base, full)
  
  data.frame(levels = paste(levels(d[[group]]), collapse = " / "), n = full$n, events = full$nevent,
             df = lr$Df[2], chi2 = lr$Chisq[2], p = lr$`Pr(>Chisq)`[2])
}

###### Gene-set scoring ######

# Calculate GSVA enrichment scores using a Gaussian expression kernel.
score_signature_GSVA = function(expr, gene_sets, suffix = "_GSVA", min_size = 1) {

  ## Restrict each signature to genes present in the expression matrix
  gene_sets = lapply(
    gene_sets,
    function(g) {
      intersect(unique(g), rownames(expr))
    }
  )

  ## Remove signatures with too few available genes
  gene_sets = gene_sets[lengths(gene_sets) >= min_size]

  if (length(gene_sets) == 0) {
    stop("No gene sets left after filtering.")
  }

  ## Create GSVA parameter object
  par = GSVA::gsvaParam(
    exprData = expr,
    geneSets = gene_sets,
    minSize = min_size,
    kcdf = "Gaussian"
  )

  ## Calculate GSVA scores
  ## output = signatures x samples
  scores = GSVA::gsva(
    par,
    verbose = TRUE
  )

  ## Convert to samples x signatures
  df = as.data.frame(t(scores)) %>%
    tibble::rownames_to_column("sample_id")

  colnames(df)[colnames(df) != "sample_id"] =
    paste0(
      colnames(df)[colnames(df) != "sample_id"],
      suffix
    )

  df
}

# Standardize a score to mean zero and unit SD, rejecting nonfinite or invariant values.
z_score = function(x,label="score") {
  spread = stats::sd(x)
  if (any(!is.finite(x)) || !is.finite(spread) || spread < 1e-8)
    stop(label,": nonfinite or invariant score")
  as.numeric((x-mean(x))/spread)
}

###### APM signature: gene sets, scoring and groups (CONFIG is defined in Part 2) ######

# Resolve configured gene sets and validate their membership, identifiers and overlap.
resolve_signature_pair = function(config) {
  annotation = APM80_SUBGROUP_CATALOG
  if (anyDuplicated(annotation$subgroup_id))
    stop("The subgroup catalogue contains duplicate subgroup IDs.")
  empty_negative = function(requested_id = NA_character_) {
    selected = annotation[NA_integer_, , drop = FALSE]
    selected$module = "None"
    selected$module_label = "No negative gene set"
    selected$subgroup_number = NA_integer_
    selected$subgroup_id = "NONE"
    selected$subgroup_label = "No negative gene set"
    selected$n_genes = 0L
    selected$genes = ""
    selected$requested_id = requested_id
    list(id = "NONE", genes = character(0), annotation = selected,
         user_defined = FALSE, absent = TRUE)
  }
  resolve_one = function(field) {
    allow_empty = identical(field, "negative_subgroup_id")
    value = config[[field]]
    if (allow_empty && !length(value)) return(empty_negative())
    if (!is.character(value) || length(value) != 1L || is.na(value))
      stop("CONFIG$", field, " must be one subgroup label/ID or gene-vector object name",
           if (allow_empty) "; use NULL or character(0) for no negative set." else ". The positive set cannot be empty.")
    value = trimws(value)
    if (!nzchar(value)) {
      if (allow_empty) return(empty_negative(value))
      stop("The positive gene-set ID cannot be blank.")
    }
    matches = which(annotation$subgroup_id == value | annotation$subgroup_label == value)
    if (length(matches) > 1L)
      stop("Ambiguous APM80 subgroup label: ", value, ". Use its unique subgroup_id.")
    if (length(matches) == 1L) {
      selected = annotation[matches, , drop = FALSE]
      id = selected$subgroup_id[1L]
      if (!id %in% names(FIXED_APM80_SUBGROUPS)) stop("Unknown frozen APM80 subgroup: ", id)
      genes = unname(FIXED_APM80_SUBGROUPS[[id]])
      if (!length(genes) || anyDuplicated(genes) || !all(genes %in% panel_genes) ||
          !setequal(genes, strsplit(selected$genes[1L], ";", fixed = TRUE)[[1L]]))
        stop("Selected subgroup disagrees with the frozen APM80 annotation: ", id)
      user_defined = FALSE
    } else {
      if (!exists(value, envir = .GlobalEnv, inherits = FALSE))
        stop("Unknown subgroup or gene-set object: ", value,
             ". Define the gene vector before CONFIG; an unknown name is not treated as an empty set.")
      genes = get(value, envir = .GlobalEnv, inherits = FALSE)
      if (allow_empty && !length(genes) && is.null(dim(genes))) return(empty_negative(value))
      if (!is.character(genes) || !is.null(dim(genes)) || !length(genes) || anyNA(genes) ||
          any(!nzchar(trimws(genes))) || any(genes != trimws(genes)) || anyDuplicated(genes))
        stop(value, " must be a nonempty character vector of unique gene symbols without NA, blanks or surrounding whitespace.")
      genes = unname(genes)
      id = value
      selected = annotation[NA_integer_, , drop = FALSE]
      selected$module = "User-defined"
      selected$module_label = "User-defined gene sets"
      selected$subgroup_number = NA_integer_
      selected$subgroup_id = id
      selected$subgroup_label = value
      selected$n_genes = length(genes)
      selected$genes = paste(genes, collapse = ";")
      user_defined = TRUE
    }
    selected$requested_id = value
    list(id = id, genes = genes, annotation = selected, user_defined = user_defined, absent = FALSE)
  }
  parts = stats::setNames(lapply(c("positive_subgroup_id", "negative_subgroup_id"), resolve_one),
                          c("positive_subgroup_id", "negative_subgroup_id"))
  ids = vapply(parts, function(x) x$id, character(1))
  plus = parts[[1L]]$genes
  minus = parts[[2L]]$genes
  has_negative = length(minus) > 0L
  if (has_negative && ids[1L] == ids[2L]) stop("Positive and negative gene sets must be different.")
  overlap = intersect(plus, minus)
  if (length(overlap)) stop("Positive and negative gene sets overlap: ", paste(overlap, collapse = ", "))
  selected = do.call(rbind, lapply(parts, function(x) x$annotation))
  selected = cbind(role = c("positive", "negative"), selected)
  rownames(selected) = NULL
  signature_mode = if (has_negative) "positive_minus_negative" else "positive_only"
  id = if (has_negative) paste(ids[1L], ids[2L], sep = "__minus__") else paste0(ids[1L], "__positive_only")
  if (!any(vapply(parts, function(x) x$user_defined, logical(1)))) {
    pair_folder = if (has_negative) {
      sprintf("S%02d_minus_S%02d", match(ids[1L], annotation$subgroup_id), match(ids[2L], annotation$subgroup_id))
    } else sprintf("S%02d_positive_only", match(ids[1L], annotation$subgroup_id))
  } else {
    hash_file = tempfile(fileext = ".rds")
    on.exit(unlink(hash_file), add = TRUE)
    saveRDS(list(ids = unname(ids), plus = plus, minus = minus), hash_file, version = 2, compress = FALSE)
    pair_folder = paste0(if (has_negative) "USERPAIR_" else "USERSET_",
                         substr(unname(tools::md5sum(hash_file)), 1L, 12L))
  }
  list(plus = plus, minus = minus, signature_id = unname(id),
       signature_label = if (has_negative) paste(selected$subgroup_label, collapse = " minus ") else
         paste0(selected$subgroup_label[1L], " (positive only)"),
       positive_subgroup_id = unname(ids[1L]), negative_subgroup_id = unname(ids[2L]),
       positive_label = selected$subgroup_label[1L], negative_label = selected$subgroup_label[2L],
       positive_module_label = selected$module_label[1L], negative_module_label = selected$module_label[2L],
       selected_pair = selected, pair_folder = pair_folder, score_prefix = paste0(id, "_"),
       has_negative = has_negative, signature_mode = signature_mode,
       n_positive = length(plus), n_negative = length(minus))
}

# APM score for one cohort: z[GSVA(positive) - GSVA(negative)], or z[GSVA(positive)] when there is no negative set.
# GSVA: Gaussian kernel; background = every gene in expr (GSVA drops constant genes). 
# Also used for SCAN-B, MATADOR and TEX.
# The score is standardized over the samples (columns) of expr. Returns a vector named by sample ID.
score_signature = function(expr, plus, minus = character(0), label = "APM score") {
  genes = c(plus, minus)
  if (!length(plus)) stop(label, ": the positive gene set is empty.")
  if (!all(genes %in% rownames(expr)))
    stop(label, ": genes missing from the expression matrix: ", paste(setdiff(genes, rownames(expr)), collapse = ", "))
  spread = apply(expr[genes, , drop = FALSE], 1L, stats::sd)
  if (any(!is.finite(spread) | spread < 1e-8))
    stop(label, ": constant signature gene(s): ", paste(genes[!is.finite(spread) | spread < 1e-8], collapse = ", "))
  sets = list(POSITIVE = plus)
  if (length(minus)) sets$NEGATIVE = minus
  gsva = score_signature_GSVA(expr, sets, suffix = "", min_size = min(lengths(sets)))
  gsva = gsva[match(colnames(expr), gsva$sample_id), , drop = FALSE]
  raw = if (length(minus)) gsva$POSITIVE - gsva$NEGATIVE else gsva$POSITIVE
  stats::setNames(z_score(raw, label), colnames(expr))
}

# APM group from a score: APM-high = score >= cutoff (default: the median of the scores given)
apm_group = function(score, cutoff = stats::median(score, na.rm = TRUE)) {
  factor(ifelse(score >= cutoff, "APM-high", "APM-low"), levels = c("APM-low", "APM-high"))
}

###### Proportional-hazards checks ######

# Complete cases for the PH checks of the APM score (PANTHER endpoints)
prepare_ph_data = function(endpoint, data) {

  time_var   = paste0("time", endpoint)
  status_var = paste0(endpoint, ".status")

  vars_needed = c(
    time_var,
    status_var,
    "TREAT",
    "APM_score",
    "ERpos",
    "newHER2",
    "tumdiam",
    "posnodes",
    "histgrade"
  )

  missing_vars = setdiff(vars_needed, names(data))

  if (length(missing_vars) > 0) {
    stop(
      "Missing variable(s) in data for ",
      endpoint,
      ": ",
      paste(missing_vars, collapse = ", ")
    )
  }

  dat = data.frame(
    time      = as.numeric(data[[time_var]]) / 365.25,
    status    = as.numeric(data[[status_var]]),
    TREAT     = data[["TREAT"]],
    APM_score = as.numeric(data[["APM_score"]]),
    ERpos     = data[["ERpos"]],
    newHER2   = data[["newHER2"]],
    tumdiam   = as.numeric(data[["tumdiam"]]),
    posnodes  = as.numeric(data[["posnodes"]]),
    histgrade = data[["histgrade"]])

  dat = dat[
    complete.cases(dat) &
      dat$time > 0 &
      dat$status %in% c(0, 1),
  ]

  dat$TREAT = factor(dat$TREAT)
  dat$ERpos = factor(dat$ERpos)
  dat$newHER2 = factor(dat$newHER2)
  dat$histgrade = factor(dat$histgrade)

  dat
}

# Same model with ER as a covariate instead of a stratum: tests whether ER violates PH
check_ph_er_covariate = function(endpoint, data) {
  
  dat = prepare_ph_data(endpoint, data)
  
  fit = coxph(Surv(time, status) ~ TREAT * APM_score + ERpos + newHER2 + scale(tumdiam) + scale(posnodes) + histgrade, data = dat, x = TRUE, y = TRUE, model = TRUE)
  
  ph = cox.zph(fit, transform = "km")
  
  ph_table = as.data.frame(ph$table) %>%
    rownames_to_column("term") %>%
    mutate(endpoint = endpoint, model = "ER as covariate") %>%
    select(model, endpoint, term, chisq, df, p)
  
  list(
    endpoint = endpoint,
    model_fit = fit,
    ph_test = ph,
    table = ph_table)
}

# Assess proportional hazards in the adjusted treatment-by-APM model stratified by ER status.
check_ph_er_stratified = function(endpoint, data) {

  dat = prepare_ph_data(endpoint, data)

  cat("\n", endpoint, ": n =", nrow(dat), ", events =", sum(dat$status == 1), "\n")

  fit = coxph(Surv(time, status) ~ TREAT * APM_score + strata(ERpos) + newHER2 + scale(tumdiam) + scale(posnodes) + histgrade, data = dat, x = TRUE, y = TRUE, model = TRUE)

  ph = cox.zph(fit, transform = "km")

  ph_table = as.data.frame(ph$table) %>%
    rownames_to_column("term") %>%
    mutate(endpoint = endpoint, model = "ER-stratified") %>%
    select(model, endpoint, term, chisq, df, p)

  list(
    endpoint = endpoint,
    model_fit = fit,
    ph_test = ph,
    table = ph_table)
}

###### Kaplan-Meier curves ######

# S3 method used implicitly: lets ggplot2::ggsave() save the ggsurvplot objects returned by plot_km()
grid.draw.ggsurvplot = function(x){
  survminer:::print.ggsurvplot(x, newpage = FALSE)
}

# Format p-values for plot annotations.
p_lab = function(p) { f = format_p_regular(p); if (startsWith(f, "<")) paste("p", sub("<", "< ", f)) else paste("p =", f) }

# PANTHER (BCRFS, EFS, DDFS, OS) and SCAN-B (BCFi, RFi, DRFi, OS); saved to out_dir
plot_km = function(dataset, predictor, predictor_name, cohort, surv_endpoint = c("BCRFS", "EFS", "DDFS", "OS", "BCFi", "RFi", "DRFi"), colors, out_dir = output_dir) {

  surv_endpoint = match.arg(surv_endpoint)

  # ----- endpoint formula helper (always returns a real formula)
  surv_formula = function(group_var) {
    time_var = switch(surv_endpoint,
                      "BCRFS" = "timeBCRFS",
                      "EFS"   = "timeEFS",
                      "DDFS"  = "timeDDFS",
                      "OS"    = "timeOS",
                      "BCFi"  = "BCFi_days",
                      "RFi"   = "RFi_days",
                      "DRFi"  = "DRFi_days")

    status_var = switch(surv_endpoint,
                        "BCRFS" = "BCRFS.status",
                        "EFS"   = "EFS.status",
                        "DDFS"  = "DDFS.status",
                        "OS"    = "OS.status",
                        "BCFi"  = "BCFi_event",
                        "RFi"   = "RFi_event",
                        "DRFi"  = "DRFi_event")

    as.formula(paste0("Surv(", time_var, "/365.25, ", status_var, ") ~ ", group_var))
  }

  # Use 7-year estimates for BCFi, 10-year estimates for all other endpoints
  surv_time = ifelse(surv_endpoint == "BCFi", 7, 10)

  d = dataset
  stopifnot(predictor %in% names(d))

  # grouping var
  d$grp = factor(d[[predictor]])

  # KM
  fml_km = surv_formula("grp")
  fit    = survival::survfit(fml_km, data = d)
  fit$call$formula = fml_km

  K = if (is.null(fit$strata)) 1L else length(fit$strata)
  if (!(K %in% c(2L, 4L))) stop("This plot_km supports K = 2 or K = 4 groups; got K = ", K)

  # labels in USER-DEFINED order (from colors)
  if (is.null(names(colors))) { stop("colors must be a NAMED vector to enforce ordering") }

  labs = names(colors)

  # relevel grp explicitly to enforce order everywhere
  d$grp = factor(d$grp, levels = labs)

  # refit KM with enforced order
  fit = survival::survfit(fml_km, data = d)
  fit$call$formula = fml_km

  n_per = as.integer(fit$n)

  # align palette to labs (works if colors is named by group)
  pal = as.character(colors[labs])

  # Survival summary at selected timepoint:
  # 7 years for BCFi, 10 years for all other endpoints
  ss = summary(fit, times = surv_time, extend = TRUE)
  ss_group = if (!is.null(ss$strata)) sub("^grp=", "", ss$strata) else labs[1]
  surv_df = data.frame(
    grp   = ss_group,
    surv  = ss$surv,
    lower = ss$lower,
    upper = ss$upper,
    stringsAsFactors = FALSE)
  # ensure order matches labs
  surv_df = surv_df[match(labs, surv_df$grp), , drop = FALSE]
  surv_df$n = n_per

  # rich surv text for K=2 or K=4
  surv_text = paste0(mapply(function(g, n, s, lo, hi, col) {
    
    paste0( if (cohort != 'All' && !grepl("SCAN-B", cohort)) {
      paste0("<b style='color:", col, ";'>", cohort, " in ", g, "</b>") } else { paste0("<b style='color:", col, ";'>", g, "</b>")},
             "<span style='color:", col, ";'> (n = ", n, "): </span>", 
      formatC(s * 100, format = "f", digits = 1), "% (", formatC(lo * 100, format = "f", digits = 1), "–", formatC(hi * 100, format = "f", digits = 1), ")") },
    
    surv_df$grp, surv_df$n, surv_df$surv, surv_df$lower, surv_df$upper, pal, SIMPLIFY = TRUE), collapse = "<br>")

  # log-rank
  lr  = survival::survdiff(fml_km, data = d)
  pLR = 1 - stats::pchisq(lr$chisq, df = length(lr$n) - 1)
  pval_text = paste0("Log-rank, ", p_lab(pLR))

  # adjusted Cox
  if (grepl("SCAN-B", cohort)) {
    d$Size.mm_z = as.numeric(scale(d$Size.mm))
    fml_cox = update(fml_km, . ~ . + strata(ER) + HER2 + Size.mm_z + LN.spec + Grade_cat)
  } else {
    fml_cox = update(fml_km, . ~ . + strata(ERpos) + newHER2 + scale(tumdiam) + scale(posnodes) + histgrade)
  }

  cox = survival::coxph(fml_cox, data = d)

  if (K == 2L) {
    cs = summary(cox)
    ridx = grep("^grp", rownames(cs$coefficients))[1]
    hr       = formatC(cs$coefficients[ridx, "exp(coef)"], format = "f", digits = 2)
    ci_lower = formatC(cs$conf.int[ridx, "lower .95"], format = "f", digits = 2)
    ci_upper = formatC(cs$conf.int[ridx, "upper .95"], format = "f", digits = 2)
    pval_hr  = format_p_regular(cs$coefficients[ridx, "Pr(>|z|)"])
    hr_text   = paste0("HR-adj.: ", hr, " (", ci_lower, "–", ci_upper, "), ", p_lab(cs$coefficients[ridx, "Pr(>|z|)"]))
  } else { # K == 4
    lrt = drop1(cox, scope = ~ grp, test = "Chisq")
    p_global = lrt["grp", "Pr(>Chi)"]
    hr_text   = paste0("Cox adj. (global grp), ", p_lab(p_global))
  }
  
  # ggsurvplot
  plot = survminer::ggsurvplot(
    fit,
    data = d,
    size = 0.5, # KM curves size

    censor = TRUE,
    censor.shape = "+",
    censor.size = 3,

    linetype = rep("solid", K),
    ggtheme = ggplot2::theme_classic(),
    surv.plot.height = 0.75,

    conf.int = TRUE,
    conf.int.alpha = 0.05,

    pval = FALSE,
    pval.method = FALSE,

    risk.table = TRUE,
    risk.table.title = '',
    risk.table.height = 0.15,
    risk.table.fontsize = 3,
    risk.table.y.text = TRUE, # risk table risk group names

    tables.col = "strata",
    tables.theme = survminer::theme_survminer() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),

        axis.title.x = ggplot2::element_blank(),
        axis.text.x  = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),

        axis.text.y  = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()),

    break.time.by = 2,
    xlim = c(0, 10),
    ylim = c(0.5, 1),
    surv.scale = "percent",

    legend = "bottom",
    legend.title = "",
    legend.labs  = labs,
    palette      = unname(pal),

    xlab = "Follow-up time (years)",
    ylab = paste0(surv_endpoint, " probability"),

    font.x = c(12, face = "bold"),
    font.y = c(12, face = "bold"),
    font.title = c(18),
    font.tickslab = c(10),
    font.legend = c(8))

  plot$table = plot$table +
    ggplot2::labs(title = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = -5, r = 5, b = 5, l = 5))

  # Dashed horizontal lines at selected survival estimate + vertical line
  # BCFi: 7 years; all other endpoints: 10 years
  svals = surv_df$surv
  plot$plot = plot$plot + ggplot2::annotate("segment", x = surv_time, xend = surv_time, y = 0.5, yend = max(svals, na.rm = TRUE), linetype = "dashed", linewidth = 0.2, color = "gray70")
  for (s in svals) {
    plot$plot = plot$plot + ggplot2::annotate("segment", x = -1, xend = surv_time, y = s, yend = s, linetype = "dashed", linewidth = 0.2, color = "gray70")
  }

  # annotation blocks
  if (K == 2L) {

    plot$plot = plot$plot +
      ggplot2::annotate("text", x = 0.2, y = 0.72, label = paste0("At ", surv_time, " years:"), size = 4.5, hjust = 0, fontface = "bold") +
      ggtext::geom_richtext(data = data.frame(x = 0.1, y = 0.65, label = surv_text), ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 4.5, fill = NA, label.color = NA) +
      ggplot2::annotate("text", x = 0.2, y = 0.58, label = pval_text, size = 4.5, hjust = 0) +
      ggplot2::annotate("text", x = 0.2, y = 0.52, label = hr_text, size = 4.5, hjust = 0)

  } else { # K == 4

    plot$plot = plot$plot +
      ggplot2::annotate("text", x = 0.2, y = 0.73, label = paste0("At ", surv_time, " years:"), size = 4.5, hjust = 0, fontface = "bold") +
      ggtext::geom_richtext(data = data.frame(x = 0.1, y = 0.64, label = surv_text), ggplot2::aes(x = x, y = y, label = label), inherit.aes = FALSE, hjust = 0, size = 4.5, fill = NA, label.color = NA) +
      ggplot2::annotate("text", x = 0.2, y = 0.55, label = pval_text, size = 4.5, hjust = 0) +
      ggplot2::annotate("text", x = 0.2, y = 0.51, label = hr_text, size = 4.5, hjust = 0)
  }

  print(plot)

  ggplot2::ggsave(filename = paste0("KM_", cohort, "_", predictor, "_", surv_endpoint, ".png"),
                  plot = plot, width = 2500, height = 1500, units = "px", bg = "white",
                  path = out_dir)

  return(plot)
}

###### Table 1 ######

# Table 1 by treatment arm (one gtsummary table per cohort definition)
make_table1 = function(dat) {
  tbl_summary(
    dat,
    include = c(age, mensstat, optype, tumdiam_cat, posnodes_cat, histgrade, Clinical_Subtype),
    by = TREAT,
    missing = "ifany",
    label = list(
      age = "Age (years)",
      mensstat = "Menopausal Status",
      optype = "Operation Type",
      tumdiam_cat = "Tumor Diameter (cm)",
      posnodes_cat = "Positive Nodes",
      histgrade = "Histological Grade",
      Clinical_Subtype = "Clinical Subtype")) |>
    add_n() |>
    add_p() |>
    modify_header(label = "**Variable**") |>
    bold_labels()
}

###### 21-gene recurrence score (RS) ######

# 21-gene recurrence score: genefu's oncotypedx() with CTSL2 renamed CTSV (needs data(sig.oncotypedx))
OncotypeDX = function (data, annot, do.mapping = FALSE, mapping, do.scaling = TRUE, verbose = FALSE) {
  if (!exists('sig.oncotypedx')) data(sig.oncotypedx, envir = environment())
  sig2 = sig.oncotypedx[sig.oncotypedx[, 'group'] != 'reference', , drop = FALSE]
  dimnames(sig2)[[1]] = sig2[, 'symbol']
  gt = nrow(sig2)
  if (do.mapping) {
    gid1 = as.numeric(as.character(sig2[, 'EntrezGene.ID']))
    names(gid1) = dimnames(sig2)[[1]]
    gid2 = as.numeric(as.character(annot[, 'EntrezGene.ID']))
    names(gid2) = dimnames(annot)[[1]]
    rm.ix = is.na(gid1) | duplicated(gid1)
    gid1 = gid1[!rm.ix]
    rr = geneid.map(geneid1 = gid2, data1 = data, geneid2 = gid1, verbose = FALSE)
    gm = length(rr$geneid2)
    mymapping = c(mapped = gm, total = gt)
    if (length(rr$geneid1) != gt) {
      res = rep(NA, nrow(data))
      names(res) = dimnames(data)[[1]]
      if (verbose) {
        message(sprintf('probe candidates: %i/%i', gm, gt))
      }
      return(list(score = res, risk = res, mapping = mymapping, probe = NA))
    }
    gid1 = rr$geneid2
    gid2 = rr$geneid1
    data = rr$data1
    myprobe = cbind(probe = names(gid1), EntrezGene.ID = gid1, new.probe = names(gid2))
    dimnames(data)[[2]] = names(gid2) = names(gid1)
  }
  else {
    myprobe = NA
    data = data[, intersect(dimnames(sig2)[[1]], dimnames(data)[[2]])]
    gm = ncol(data)
    mymapping = c(mapped = gm, total = gt)
    if (nrow(sig2) != ncol(data)) {
      res = rep(NA, nrow(data))
      names(res) = dimnames(data)[[1]]
      if (verbose) {
        message(sprintf('probe candidates: %i/%i', ncol(data), gt))
      }
      return(list(score = res, risk = res, mapping = mymapping, probe = myprobe))
    }
  }
  dimnames(data)[[2]] = dimnames(sig2)[[1]] = sig2[, 'symbol']
  if (do.scaling) {
    data = apply(data, 2, function(x) {
      xx = (x - min(x, na.rm = TRUE))/(max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
      return(xx * 15)
    })
  }
  else if (max(data, na.rm = TRUE) > 20 || min(data, na.rm = TRUE) < -5) {
    warning('The max and min values of your data suggest it is not already scaled...\n\t\t\tIf it is please ignore this message, otherwise set `do.scaling=TRUE` to scale in the function call.')
  }
  cc.ix = complete.cases(data)
  rs = rs.unscaled = rsrisk = NULL
  for (i in 1:nrow(data)) {
    if (cc.ix[i]) {
      grb7.gs = 0.9 * data[i, 'GRB7'] + 0.1 * data[i, 'ERBB2']
      if (grb7.gs < 8) {
        grb7.gs = 8
      }
      er.gs = (0.8 * data[i, 'ESR1'] + 1.2 * data[i, 'PGR'] + data[i, 'BCL2'] + data[i, 'SCUBE2']) / 4
      proliferation.gs = (data[i, 'BIRC5'] + data[i, 'MKI67'] + data[i, 'MYBL2'] + data[i, 'CCNB1'] + data[i, 'AURKA']) / 5
      if (proliferation.gs < 6.5) {
        proliferation.gs = 6.5
      }
      invasion.gs = (data[i, 'CTSV'] + data[i, 'MMP11'])/2
      rsu = 0.47 * (grb7.gs) - 0.34 * (er.gs) + 1.04 * (proliferation.gs) + 0.1 * (invasion.gs) + 0.05 * data[i, 'CD68'] - 0.08 * data[i, 'GSTM1'] - 0.07 * data[i, 'BAG1']
      rsu2 = rsu
      if (rsu >= 0 & rsu <= 100) {
        rsu = 20 * (rsu - 6.7)
      }
      if (rsu < 0) {
        rsu = 0
      }
      if (rsu > 100) {
        rsu = 100
      }
      if (rsu < 18) {
        rsr = 0
      }
      if (rsu >= 18 & rsu < 31) {
        rsr = 0.5
      }
      if (rsu >= 31) {
        rsr = 1
      }
    }
    else {
      rsu = rsr = rsu2 = NA
    }
    rs.unscaled = c(rs.unscaled, rsu2)
    rs = c(rs, rsu)
    rsrisk = c(rsrisk, rsr)
  }
  names(rs) = names(rs.unscaled) = names(rsrisk) = dimnames(data)[[1]]
  return(list(score = rs, risk = rsrisk, mapping = mymapping, probe = myprobe))
}

###### Immune deconvolution boxplots ######

# Immune cell scores by APM group: boxplots, Wilcoxon (or t) tests, delta mean and Cliff's delta
plot_apm_immune_scores = function(data, score_cols = NULL, score_prefix = NULL, group_var = "APM_group", group_levels = c("APM-low", "APM-high"), colors = group_col_apm,
                                  x_label = "Cell Population", y_label = "Score", legend_title = "APM group", title = NULL, feature_order = NULL,
                                  test_method = c("wilcox", "t_test"), display_p = c("p", "p.adj"), p_adjust_method = "BH", y_multiplier = 1.20, dodge_width = 0.75,
                                  box_width = 0.65, point_size = 0.5, point_alpha = 0.4, base_size = 14, x_text_angle = 45) {

  # 1. Check inputs
  test_method = match.arg(test_method)
  display_p = match.arg(display_p)

  if (!group_var %in% colnames(data)) { stop("Grouping variable '", group_var, "' was not found in the dataset.") }

  if (is.null(score_cols) && is.null(score_prefix)) { stop("Provide either 'score_cols' or 'score_prefix'.") }

  if (!is.null(score_prefix)) { score_cols = grep(paste0("^", score_prefix), colnames(data), value = TRUE) }

  if (length(score_cols) == 0) { stop("No score columns were identified.") }

  missing_score_cols = setdiff(score_cols, colnames(data))

  if (length(missing_score_cols) > 0) { stop("The following score columns were not found: ", paste(missing_score_cols, collapse = ", ")) }

  if (!all(group_levels %in% unique(data[[group_var]]))) { stop("Not all specified group levels are present in '", group_var, "': ", paste(group_levels, collapse = ", ")) }


  # 2. Prepare long-format plotting data
  plot_df = data %>%
    dplyr::select(
      dplyr::all_of(c(group_var, score_cols))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(score_cols),
      names_to = "Cell_Population",
      values_to = "Score") %>%
    dplyr::filter(
      !is.na(.data[[group_var]]),
      !is.na(Score)) %>%
    dplyr::mutate(
      !!group_var := factor(
        .data[[group_var]],
        levels = group_levels),
      Cell_Population = if (!is.null(score_prefix)) {
        sub(paste0("^", score_prefix), "", Cell_Population)
      } else {
        Cell_Population
      },
      Cell_Population = gsub("_", " ", Cell_Population, fixed = TRUE))

  if (!is.null(feature_order)) {
    plot_df = plot_df %>% dplyr::mutate(Cell_Population = factor(Cell_Population, levels = feature_order))
  }

  test_formula = stats::as.formula(paste0("Score ~ ", group_var))

  # 3. Statistical tests + effect sizes
  if (test_method == "wilcox") {

    pvals = plot_df %>%
      dplyr::group_by(Cell_Population) %>%
      rstatix::wilcox_test(test_formula) %>%
      dplyr::ungroup()
  } else {
    pvals = plot_df %>%
      dplyr::group_by(Cell_Population) %>%
      rstatix::t_test(test_formula) %>%
      dplyr::ungroup()
  }

  pvals = pvals %>%
    rstatix::adjust_pvalue(method = p_adjust_method) %>%
    rstatix::add_xy_position(
      x = "Cell_Population",
      dodge = dodge_width,
      data = plot_df,
      formula = test_formula)

  # Δ mean: APM-high minus APM-low
  mean_stats = plot_df %>%
    dplyr::group_by(Cell_Population, .data[[group_var]]) %>%
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = .data[[group_var]],
      values_from = mean_score) %>%
    dplyr::mutate(
      delta_mean = .data[[group_levels[2]]] - .data[[group_levels[1]]])

  # Cliff's delta: direction follows group_levels[2] vs group_levels[1]
  cliff_stats = plot_df %>%
    dplyr::group_by(Cell_Population) %>%
    dplyr::group_modify(~{

      x = .x$Score[.x[[group_var]] == group_levels[2]]
      y = .x$Score[.x[[group_var]] == group_levels[1]]

      cd = effsize::cliff.delta(x, y)

      tibble::tibble(
        cliff_delta = unname(cd$estimate),
        cliff_magnitude = cd$magnitude)}) %>%
    dplyr::ungroup()

  # 4. Population-specific bracket positions + multiline labels
  y_positions = plot_df %>%
    dplyr::group_by(Cell_Population) %>%
    dplyr::summarise(
      y.position = dplyr::if_else(
        max(Score, na.rm = TRUE) > 0,
        y_multiplier * max(Score, na.rm = TRUE),
        max(Score, na.rm = TRUE) +
          (y_multiplier - 1) * diff(range(Score, na.rm = TRUE))), .groups = "drop")

  pvals = pvals %>%
    dplyr::select(-y.position) %>%
    dplyr::left_join(
      y_positions,
      by = "Cell_Population") %>%
    dplyr::left_join(
      mean_stats %>% dplyr::select(Cell_Population, delta_mean),
      by = "Cell_Population") %>%
    dplyr::left_join(
      cliff_stats %>% dplyr::select(Cell_Population, cliff_delta, cliff_magnitude),
      by = "Cell_Population") %>%
    dplyr::mutate(
      p_text = dplyr::case_when(
        display_p == "p" & p < 0.0001 ~ "p < 0.0001",
        display_p == "p" ~ paste0("p = ", format_p_regular(p)),
        display_p == "p.adj" & p.adj < 0.0001 ~ "P-adj. < 0.0001",
        display_p == "p.adj" ~ paste0("P-adj. = ", format_p_regular(p.adj))),

      p_label = paste0("\u0394mean = ", sprintf("%.2f", round(delta_mean, 2) + 0), "\n", p_text, "\n", "Cliff's \u03b4 = ", sprintf("%.2f", cliff_delta)))

  # 5. Plot
  p = ggplot(
    plot_df,
    aes(
      x = Cell_Population,
      y = Score,
      fill = .data[[group_var]])) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = box_width,
      outlier.shape = NA,
      alpha = 0.80,
      linewidth = 0.55) +
    geom_jitter(
      aes(color = .data[[group_var]]),
      position = position_jitterdodge(
        jitter.width = 0.15,
        dodge.width = dodge_width),
      size = point_size,
      alpha = point_alpha) +
    ggpubr::stat_pvalue_manual(
      pvals,
      label = "p_label",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.012,
      bracket.size = 0.45,
      size = 3.0,
      vjust = -0.25,
      inherit.aes = FALSE) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.18))) +
    labs(
      x = x_label,
      y = y_label,
      fill = legend_title,
      color = legend_title,
      title = title) +
    theme_bw(base_size = base_size) +
    theme(
      axis.text.x = element_text(
        angle = x_text_angle,
        hjust = ifelse(x_text_angle == 0, 0.5, 1),
        vjust = ifelse(x_text_angle == 0, 0.5, 1)),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank())

  # 6. Return plot, statistics and plotting data
  return(list(plot = p, statistics = pvals, data = plot_df))
}

###### Head-to-head Cox models ######

# Head-to-head Cox models (APM score, SET2,3, SET ER/PR): prognostic HR and treatment interaction per marker, cohort, endpoint and model, in identical complete-case patients within each comparison set.
get_biomarker_result = function(marker, cohort, data, endpoint, model_type = c("Univariate", "Adjusted"), matched_vars = character(0), 
                                adjustment_terms = character(0), adjustment_vars = character(0), standardize = FALSE) {
  
  model_type = match.arg(model_type)
  
  time_var  = paste0("time", endpoint)
  event_var = paste0(endpoint, ".status")
  
  if (model_type == "Univariate") {
    current_adjustment_terms = character(0)
    current_adjustment_vars  = character(0)
  } else {
    current_adjustment_terms = adjustment_terms
    current_adjustment_vars  = adjustment_vars
  }
  
  # matched_vars defines the comparison population;
  # adjustment_vars defines only the model-specific covariates
  vars_needed = unique(c(time_var, event_var, marker, matched_vars, current_adjustment_vars))
  
  missing_vars = setdiff(vars_needed, colnames(data))
  
  if (length(missing_vars) > 0) { stop("Missing variables in ", cohort, " for ", endpoint, ": ", paste(missing_vars, collapse = ", ")) }
  
  model_data = data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na()
  
  if (standardize) model_data[[marker]] = as.numeric(scale(model_data[[marker]]))
  
  surv_term = paste0("Surv(", time_var, "/365.25, ", event_var, ")")
  
  adjustment_string = if (length(current_adjustment_terms) > 0) {
    paste0(" + ", paste(current_adjustment_terms, collapse = " + "))
  } else { "" }
  
  model_formula = stats::as.formula(paste0(surv_term, " ~ ", marker, adjustment_string))
  
  cox_model = survival::coxph(formula = model_formula, data = model_data)
  
  biomarker_row = broom::tidy(
    cox_model,
    exponentiate = TRUE,
    conf.int = TRUE) %>%
    dplyr::filter(term == marker)
  
  if (nrow(biomarker_row) != 1) { stop("Expected exactly one biomarker coefficient for ", marker, " in ", cohort, " for ", endpoint, " (", model_type, "); found ", nrow(biomarker_row), ".") }
  
  biomarker_row %>%
    dplyr::transmute(
      cohort = cohort,
      marker = marker,
      endpoint = endpoint,
      model_type = model_type,
      adjustment = ifelse(model_type == "Univariate", "None", paste(current_adjustment_terms, collapse = " + ")),
      n = nrow(model_data),
      events = sum(model_data[[event_var]] == 1),
      biomarker_HR = estimate,
      HR_unit = if (standardize) "per SD in the analysed patients" else "per SD in all RNA-matched patients",
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value)
}

# Extract the treatment-by-marker interaction HR, confidence interval, Wald p-value and LRT p-value.
get_interaction_result = function(marker, cohort, data, endpoint, model_type = c("Unadjusted interaction", "Adjusted interaction"), matched_vars = character(0), 
                                  adjustment_terms = character(0), adjustment_vars = character(0), standardize = FALSE) {
  
  model_type = match.arg(model_type)
  
  time_var  = paste0("time", endpoint)
  event_var = paste0(endpoint, ".status")
  
  if (model_type == "Unadjusted interaction") {
    current_adjustment_terms = character(0)
    current_adjustment_vars  = character(0)
  } else {
    current_adjustment_terms = adjustment_terms
    current_adjustment_vars  = adjustment_vars
  }
  
  vars_needed = unique(c(
    time_var,
    event_var,
    "TREAT",
    marker,
    matched_vars,
    current_adjustment_vars))
  
  missing_vars = setdiff(vars_needed, colnames(data))
  
  if (length(missing_vars) > 0) { stop("Missing variables in ", cohort, " for ", endpoint, ": ", paste(missing_vars, collapse = ", ")) }
  
  model_data = data %>%
    dplyr::select(dplyr::all_of(vars_needed)) %>%
    tidyr::drop_na()
  
  if (standardize) model_data[[marker]] = as.numeric(scale(model_data[[marker]]))
  
  if (dplyr::n_distinct(model_data$TREAT) != 2) { stop("TREAT must contain exactly two treatment groups in ", cohort, " for ", endpoint, " (", model_type, ").") }
  
  surv_term = paste0("Surv(", time_var, "/365.25, ", event_var, ")")
  
  adjustment_string = if (length(current_adjustment_terms) > 0) {
    paste0(" + ", paste(current_adjustment_terms, collapse = " + "))
  } else { "" }
  
  base_formula = stats::as.formula(paste0(surv_term, " ~ TREAT + ", marker, adjustment_string))
  
  interaction_formula = stats::as.formula(paste0(surv_term, " ~ TREAT * ", marker, adjustment_string))
  
  cox_base = survival::coxph(formula = base_formula, data = model_data)
  
  cox_int = survival::coxph(formula = interaction_formula, data = model_data)
  
  lr_test = stats::anova(cox_base, cox_int, test = "LRT")
  
  interaction_lrt_p = lr_test[2, "Pr(>|Chi|)"]
  
  interaction_row = broom::tidy(
    cox_int,
    exponentiate = TRUE,
    conf.int = TRUE) %>%
    dplyr::filter(grepl(":", term) & grepl("TREAT", term) & grepl(marker, term))
  
  if (nrow(interaction_row) != 1) { stop("Expected exactly one treatment-by-marker interaction term for ", marker, " in ", cohort, " for ", endpoint, " (", model_type, "); found ", nrow(interaction_row), ".") }
  
  interaction_row %>%
    dplyr::transmute(
      cohort = cohort,
      marker = marker,
      endpoint = endpoint,
      model_type = model_type,
      adjustment = ifelse(model_type == "Unadjusted interaction", "None", paste(current_adjustment_terms, collapse = " + ")),
      n = nrow(model_data),
      events = sum(model_data[[event_var]] == 1),
      interaction_HR = estimate,
      HR_unit = if (standardize) "per SD in the analysed patients" else "per SD in all RNA-matched patients",
      CI_low = conf.low,
      CI_high = conf.high,
      Wald_p = p.value,
      LRT_interaction_p = interaction_lrt_p)
}

# Fit univariate and adjusted prognostic models across the specified analysis sets and endpoints.
run_biomarker_models = function(analysis_sets, endpoints) {
  
  purrr::map_dfr(endpoints, function(current_endpoint) {
    
    purrr::imap_dfr(analysis_sets, function(current_set, set_name) {
      
      purrr::map_dfr(c("Univariate", "Adjusted"), function(current_model_type) {
        
        get_biomarker_result(
          marker = current_set$marker,
          cohort = current_set$cohort,
          data = current_set$data,
          endpoint = current_endpoint,
          model_type = current_model_type,
          matched_vars = current_set$matched_vars,
          adjustment_terms = current_set$adjustment_terms,
          adjustment_vars = current_set$adjustment_vars,
          standardize = isTRUE(current_set$standardize)) }) }) }) }

# Fit unadjusted and adjusted treatment-interaction models across analysis sets and endpoints.
run_interaction_models = function(analysis_sets, endpoints) {
  
  purrr::map_dfr(endpoints, function(current_endpoint) {
    
    purrr::imap_dfr(analysis_sets, function(current_set, set_name) {
      
      purrr::map_dfr(c("Unadjusted interaction", "Adjusted interaction"), function(current_model_type) {
        
        get_interaction_result(
          marker = current_set$marker,
          cohort = current_set$cohort,
          data = current_set$data,
          endpoint = current_endpoint,
          model_type = current_model_type,
          matched_vars = current_set$matched_vars,
          adjustment_terms = current_set$adjustment_terms,
          adjustment_vars = current_set$adjustment_vars,
          standardize = isTRUE(current_set$standardize)) }) }) }) }

# Arrange prognostic results with univariate and adjusted estimates in adjacent columns.
make_biomarker_wide = function(results, endpoints) {
  
  results %>%
    dplyr::mutate(
      endpoint = factor(endpoint, levels = endpoints),
      model_type = factor(model_type, levels = c("Univariate", "Adjusted")),
      `HR (95% CI)` = paste0(sprintf("%.2f", biomarker_HR), " (", sprintf("%.2f", CI_low), "–", sprintf("%.2f", CI_high), ")"),
      `p-value` = p_value) %>%
    dplyr::select(
      cohort,
      marker,
      HR_unit,
      endpoint,
      model_type,
      adjustment,
      n,
      events,
      `HR (95% CI)`,
      `p-value`) %>%
    tidyr::pivot_wider(
      names_from = model_type,
      values_from = c(
        adjustment,
        n,
        events,
        `HR (95% CI)`,
        `p-value`),
      names_glue = "{model_type}_{.value}") %>%
    dplyr::arrange(cohort, marker, endpoint)
}

# Arrange interaction results with unadjusted and adjusted estimates in adjacent columns.
make_interaction_wide = function(results, endpoints) {
  
  results %>%
    dplyr::mutate(
      endpoint = factor(endpoint, levels = endpoints),
      model_type = factor(model_type, levels = c("Unadjusted interaction", "Adjusted interaction")),
      `Interaction HR (95% CI)` = paste0(sprintf("%.2f", interaction_HR), " (", sprintf("%.2f", CI_low), "–", sprintf("%.2f", CI_high), ")"),
      `Wald p-value` = Wald_p,
      `LRT interaction p-value` = LRT_interaction_p) %>%
    dplyr::select(
      cohort,
      marker,
      HR_unit,
      endpoint,
      model_type,
      adjustment,
      n,
      events,
      `Interaction HR (95% CI)`,
      `Wald p-value`,
      `LRT interaction p-value`) %>%
    tidyr::pivot_wider(
      names_from = model_type,
      values_from = c(
        adjustment,
        n,
        events,
        `Interaction HR (95% CI)`,
        `Wald p-value`,
        `LRT interaction p-value`),
      names_glue = "{model_type}_{.value}") %>%
    dplyr::arrange(cohort, marker, endpoint)
}

###### Gene-module scores by APM group ######

# Mean expression of gene modules by APM group (Wilcoxon tests, boxplots)
run_module_score_apm_comparison = function(expr_mat, pheno_df, gene_modules, cohort_name = "Cohort", sample_id_col = "sample_id", apm_group_col = "APM_group") {

  ## Make sure expression is matrix-like
  expr_mat = as.matrix(expr_mat)

  ## Check sample overlap
  expr_samples = colnames(expr_mat)
  pheno_samples = pheno_df[[sample_id_col]]

  n_overlap = sum(expr_samples %in% pheno_samples)
  message(cohort_name, ": ", n_overlap, "/", length(expr_samples), " expression samples found in phenotype data")

  if (n_overlap == 0) {
    stop("No matching sample IDs between expression matrix columns and phenotype sample_id column.")
  }

  ## Function to compute one module score
  compute_module_score_from_matrix = function(expr_mat, genes, score_name) {

    genes_present = intersect(genes, rownames(expr_mat))
    genes_missing = setdiff(genes, rownames(expr_mat))

    message(score_name, ": ", length(genes_present), "/", length(genes), " genes found")

    if (length(genes_missing) > 0) {
      message("Missing genes: ", paste(genes_missing, collapse = ", "))
    }

    if (length(genes_present) == 0) {
      score_df = data.frame(
        sample_id = colnames(expr_mat),
        score = NA_real_)
    } else {
      score_df = data.frame(
        sample_id = colnames(expr_mat),
        score = colMeans(expr_mat[genes_present, , drop = FALSE], na.rm = TRUE))
    }

    colnames(score_df)[2] = score_name

    return(score_df)
  }

  ## Compute all module scores
  module_scores_df = purrr::imap(
    gene_modules,
    ~ compute_module_score_from_matrix(
      expr_mat = expr_mat,
      genes = .x,
      score_name = .y)) %>%
    purrr::reduce(dplyr::left_join, by = "sample_id")

  ## Join to phenotype
  pheno_scored = pheno_df %>%
    dplyr::left_join(
      module_scores_df,
      by = setNames("sample_id", sample_id_col))

  ## Long-format data
  plot_df = pheno_scored %>%
    dplyr::select(
      dplyr::all_of(sample_id_col),
      dplyr::all_of(apm_group_col),
      dplyr::all_of(names(gene_modules))) %>%
    dplyr::rename(
      sample_id = dplyr::all_of(sample_id_col),
      APM_group = dplyr::all_of(apm_group_col)) %>%
    dplyr::filter(!is.na(APM_group)) %>%
    dplyr::mutate(
      APM_group = as.character(APM_group),
      APM_group = dplyr::case_when(
        APM_group %in% c("APM-low", "Low", "low") ~ "APM-low",
        APM_group %in% c("APM-high", "High", "high") ~ "APM-high",
        TRUE ~ APM_group),
      APM_group = factor(APM_group, levels = c("APM-low", "APM-high"))) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(names(gene_modules)),
      names_to = "Module",
      values_to = "Score") %>%
    dplyr::filter(
      !is.na(Score),
      !is.na(APM_group))

  ## Function for Wilcoxon test
  run_wilcox_one_module = function(df) {

    n_low = sum(df$APM_group == "APM-low")
    n_high = sum(df$APM_group == "APM-high")

    median_low = if (n_low > 0) {
      median(df$Score[df$APM_group == "APM-low"], na.rm = TRUE)
    } else {
      NA_real_
    }

    median_high = if (n_high > 0) {
      median(df$Score[df$APM_group == "APM-high"], na.rm = TRUE)
    } else {
      NA_real_
    }

    p_val = if (n_low > 0 & n_high > 0) {
      wilcox.test(
        Score ~ APM_group,
        data = df,
        exact = FALSE)$p.value
    } else {
      NA_real_
    }

    tibble::tibble(
      n_APM_low = n_low,
      n_APM_high = n_high,
      median_APM_low = median_low,
      median_APM_high = median_high,
      p_value = p_val)
  }

  ## Wilcoxon results
  wilcox_results = plot_df %>%
    split(.$Module) %>%
    purrr::imap_dfr(
      ~ run_wilcox_one_module(.x) %>%
        dplyr::mutate(Module = .y)) %>%
    dplyr::select(Module, everything()) %>%
    dplyr::mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      p_label = dplyr::case_when(
        is.na(p_value) ~ "p = NA",
        p_value < 0.001 ~ "p < 0.001",
        TRUE ~ paste0("p = ", signif(p_value, 2))),
      p_adj_label = dplyr::case_when(
        is.na(p_adj) ~ "adj. p = NA",
        p_adj < 0.001 ~ "adj. p < 0.001",
        TRUE ~ paste0("adj. p = ", signif(p_adj, 2))))

  ## P-value positions
  p_label_df = plot_df %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      y_position = max(Score, na.rm = TRUE) +
        0.10 * diff(range(Score, na.rm = TRUE)),
      .groups = "drop") %>%
    dplyr::left_join(
      wilcox_results %>% dplyr::select(Module, p_label),
      by = "Module") %>%
    dplyr::mutate(
      x_position = 1.5)

  ## Plot
  p = ggplot(
    plot_df,
    aes(
      x = APM_group,
      y = Score,
      fill = APM_group)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      alpha = 0.85) +
    geom_jitter(
      width = 0.15,
      size = 1.2,
      alpha = 0.5) +
    geom_text(
      data = p_label_df,
      aes(
        x = x_position,
        y = y_position,
        label = p_label),
      inherit.aes = FALSE,
      size = 4) +
    facet_wrap(
      ~ Module,
      scales = "free_y",
      ncol = 3) +
    labs(
      x = NULL,
      y = "Mean module expression",
      title = paste0(cohort_name, ": gene module scores by APM group")) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"))

  return(
    list(
      pheno_scored = pheno_scored,
      plot_df = plot_df,
      wilcox_results = wilcox_results,
      plot = p))
}

###### Expression and clinical data processing ######

# One row per gene symbol: keep the feature with the highest mean log2 expression
collapse_max_mean = function(expr, symbols) {
  expr = as.matrix(expr)
  stopifnot(nrow(expr) == length(symbols))
  ok = !is.na(symbols) & nzchar(symbols)
  expr = expr[ok, , drop = FALSE]; symbols = symbols[ok]
  o = order(symbols, -rowMeans(expr))
  keep = sort(o[!duplicated(symbols[o])])
  out = expr[keep, , drop = FALSE]
  rownames(out) = symbols[keep]
  out
}

###### Median follow-up ######

# Median follow-up by the reverse Kaplan-Meier method (censored patients become events).
# time_scale: "days" (PANTHER, SCAN-B) or "years" (MATADOR).
median_followup_reverse_km = function(data, endpoint, time_var = paste0("time", endpoint), status_var = paste0(endpoint, ".status"), time_scale = c("days", "years")) {

  time_scale = match.arg(time_scale)

  df = data %>%
    dplyr::select(
      time = dplyr::all_of(time_var),
      status = dplyr::all_of(status_var)) %>%
    dplyr::filter(!is.na(time), !is.na(status)) %>%
    dplyr::mutate(
      time_years = if (time_scale == "days") time / 365.25 else time,
      reverse_status = ifelse(status == 0, 1, 0))

  fit = survival::survfit(survival::Surv(time_years, reverse_status) ~ 1, data = df)

  q = quantile(fit, probs = c(0.25, 0.50, 0.75))$quantile

  tibble::tibble(
    endpoint = endpoint,
    N = nrow(df),
    events = sum(df$status == 1, na.rm = TRUE),
    censored = sum(df$status == 0, na.rm = TRUE),
    followup_Q1_years = unname(q[1]),
    median_followup_years = unname(q[2]),
    followup_Q3_years = unname(q[3]),
    followup_IQR = paste0(round(unname(q[1]), 1), "–", round(unname(q[3]), 1)),
    median_IQR = paste0(round(unname(q[2]), 1), " years (IQR ", round(unname(q[1]), 1), "–", round(unname(q[3]), 1), ")"))
}

###### Figure assembly ######

# Total width of a unit vector in cm
in_cm = function(u) convertWidth(sum(u), "cm", valueOnly = TRUE)

# Panel letter at the top-left of a gtable, with a small gap below
label_panel = function(g, letter, gap_cm = 0.5) {
  g = gtable_add_rows(g, unit(0.9, "cm"), pos = 0)
  g = gtable_add_grob(g, textGrob(letter, x = 0, hjust = 0, gp = gpar(fontsize = 18, fontface = "bold")), t = 1, l = 1, clip = "off", name = paste0("panel-label-", letter))
  gtable_add_rows(g, unit(gap_cm, "cm"), pos = -1)
}

# Left margin before a forest-plot grob, for panel letters
add_left_space = function(g, left_width = 1.5) {
  gridExtra::arrangeGrob(
    grid::nullGrob(),
    g,
    ncol = 2,
    widths = grid::unit.c(
      grid::unit(left_width, "cm"),
      grid::unit(1, "null")))
}

###### ----- Load data -----

# PANTHER analysis inputs, made from the raw data by the RNA-seq processing, QC and intrinsic-subtyping scripts (workflows/rnaseq_qc_and_subtyping; see data/README.md):
#   PANTHER_clinical.rds                  one row per randomized patient; sample_id is set for the 444 tumours that passed RNA-seq QC
#   PANTHER_expression_log2_TMM_TPM.rds   log2(TMM-normalized TPM + 1); gene symbols x the 444 samples
panther_clinical = readRDS(file.path(data_dir, "PANTHER", "PANTHER_clinical.rds"))
tmm_normalized_tpm_symbols = readRDS(file.path(data_dir, "PANTHER", "PANTHER_expression_log2_TMM_TPM.rds"))

panther_clinical$TREAT = factor(panther_clinical$TREAT, levels = c("FEC-D", "tddEC-D"))   # FEC-D = reference

# All randomized patients (Table 1) and the RNA-seq cohort, in the sample order of the expression matrix
all_clinical_original = panther_clinical
all_clinical = panther_clinical[match(colnames(tmm_normalized_tpm_symbols), panther_clinical$sample_id), , drop = FALSE]
rownames(all_clinical) = NULL

stopifnot(
  !anyDuplicated(panther_clinical$patno),
  !anyNA(all_clinical$patno),
  identical(as.character(all_clinical$sample_id), colnames(tmm_normalized_tpm_symbols)),
  sum(!is.na(panther_clinical$sample_id)) == ncol(tmm_normalized_tpm_symbols))

###### ----- SET2,3 signature -----
###### Compute SET2,3 signature ######

# Calculate SET2,3 and SET ER/PR scores using the external pathological-scoring function in ER-positive/HER2-negative patients.

# ER+ / HER2-
all_clinical_ERpos_HERneg = all_clinical[all_clinical$Clinical_Subtype == 'ER+/HER2-',]
tmm_normalized_tpm_symbols_ERpos_HERneg = tmm_normalized_tpm_symbols[, colnames(tmm_normalized_tpm_symbols) %in% all_clinical_ERpos_HERneg$sample_id]

all_clinical_ERpos_HERneg = set23_signature_pathological(gene_exp_data = tmm_normalized_tpm_symbols_ERpos_HERneg, clinical_data = all_clinical_ERpos_HERneg)$clinical
all_clinical_ERpos_HERneg = as.data.frame(all_clinical_ERpos_HERneg)

###### Plot HR vs.Signatures (Continuous plot) - SET signatures ######

# Plot treatment HRs across continuous SET2,3 and SET ER/PR scores in ER-positive/HER2-negative population using marker-specific adjustment.

all_clinical_adj_ERpos_HERneg = all_clinical_ERpos_HERneg %>%
  dplyr::mutate(
    Clinical_Subtype = factor(Clinical_Subtype),
    tumdiam_z = as.numeric(scale(tumdiam)),
    posnodes_z = as.numeric(scale(posnodes)),
    histgrade = factor(histgrade))

Supp_Figure_2A = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_2_3", biomarker_label = "SET 2,3", endpoint = "BCRFS", subtype_label = "ER+/HER2-",
                                              adjust_vars = c("histgrade"), adjust_terms = c("histgrade"), y_limits = c(-0.10, 4))
Supp_Figure_2B = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_2_3", biomarker_label = "SET 2,3", endpoint = "EFS", subtype_label = "ER+/HER2-",
                                              adjust_vars = c("histgrade"), adjust_terms = c("histgrade"), y_limits = c(-0.10, 4))
Supp_Figure_2C = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_2_3", biomarker_label = "SET 2,3", endpoint = "DDFS", subtype_label = "ER+/HER2-",
                                              adjust_vars = c("histgrade"), adjust_terms = c("histgrade"), y_limits = c(-0.10, 4))
Supp_Figure_2D = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_2_3", biomarker_label = "SET 2,3", endpoint = "OS", subtype_label = "ER+/HER2-",
                                              adjust_vars = c("histgrade"), adjust_terms = c("histgrade"), y_limits = c(-0.10, 4))

Supp_Figure_2E = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_ER_PR", biomarker_label = "SET ER,PR", endpoint = "BCRFS", subtype_label = "ER+/HER2-",
                                              adjust_vars = adjust_vars_no_clin_sub, adjust_terms = adjust_terms_no_clin_sub, y_limits = c(-0.10, 4))
Supp_Figure_2F = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_ER_PR", biomarker_label = "SET ER,PR", endpoint = "EFS", subtype_label = "ER+/HER2-",
                                              adjust_vars = adjust_vars_no_clin_sub, adjust_terms = adjust_terms_no_clin_sub, y_limits = c(-0.10, 4))
Supp_Figure_2G = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_ER_PR", biomarker_label = "SET ER,PR", endpoint = "DDFS", subtype_label = "ER+/HER2-",
                                              adjust_vars = adjust_vars_no_clin_sub, adjust_terms = adjust_terms_no_clin_sub, y_limits = c(-0.10, 4))
Supp_Figure_2H = make_plot_hr_single_endpoint(data = all_clinical_adj_ERpos_HERneg, biomarker = "SET_ER_PR", biomarker_label = "SET ER,PR", endpoint = "OS", subtype_label = "ER+/HER2-",
                                              adjust_vars = adjust_vars_no_clin_sub, adjust_terms = adjust_terms_no_clin_sub, y_limits = c(-0.10, 4))

###### ----- APM-related transcriptional signatures -----
###### Initialize APM gene sets ######

# Initialize the registry containing the BC360 and custom APM gene panels.
apm_gene_sets = list()

###### Antigen processing and presentation machinery (APM) - BC360 ######

# Extract the 21-gene antigen-presentation module from the BC360 annotation.
BC_360_panel = read_excel(file.path(data_dir, "annotation", "LBL-10541-03-Breast-Cancer-360-Gene-List.xlsx"), sheet = 3, skip = 1)
BC_360_panel = BC_360_panel %>% filter(!grepl("Bruker|©", Gene))
BC_360_panel = BC_360_panel[!is.na(BC_360_panel$Gene),]
BC_360_panel = BC_360_panel %>% mutate(across(everything(), ~ ifelse(. == "+", 1, ifelse(. == "-", 0, .))))
BC360_apm_genes = BC_360_panel$Gene[BC_360_panel$`Antigen Presentation` == 1]

apm_gene_sets[["APM_BC360"]] = BC360_apm_genes

###### BC360 antigen-presentation module (21 genes), applied unchanged: treatment interaction ######

# GSVA with the APM-score settings: all RNA-matched PANTHER patients; background = all finite, non-constant genes (the same background as Part 1). 
# The score is standardized over these patients; continuous (per SD) treatment-interaction LRT for BCRFS. The HLA6 submodule test is in Part 1 (01_step1_BC360_gate.csv).
bc360_ids = intersect(as.character(all_clinical$sample_id), colnames(tmm_normalized_tpm_symbols))
bc360_expr = as.matrix(tmm_normalized_tpm_symbols[, bc360_ids, drop = FALSE])
bc360_expr = bc360_expr[rowSums(!is.finite(bc360_expr)) == 0, , drop = FALSE]
bc360_expr = bc360_expr[apply(bc360_expr, 1L, stats::sd) > 1e-8, , drop = FALSE]

bc360_genes_used = intersect(BC360_apm_genes, rownames(bc360_expr))

if (length(bc360_genes_used) < length(BC360_apm_genes))
  warning("BC360 genes missing or constant: ", paste(setdiff(BC360_apm_genes, bc360_genes_used), collapse = ", "))

bc360_gsva = score_signature_GSVA(bc360_expr, list(BC360_APM_gsva = bc360_genes_used), suffix = "", min_size = 2)

z = as.numeric(scale(bc360_gsva$BC360_APM_gsva[match(bc360_ids, bc360_gsva$sample_id)]))

all_clinical$BC360_APM_gsva = z[match(as.character(all_clinical$sample_id), bc360_ids)]

bc360_int = get_linear_interaction_pvalues(data = all_clinical, time_col = "timeBCRFS", status_col = "BCRFS.status",
                                           continuous_var = "BC360_APM_gsva", adjust_vars = adjust_vars, adjust_terms = adjust_terms)

bc360_continuous_results = data.frame(biomarker = "BC360_APM_gsva", n = bc360_int$n, events = bc360_int$events, p_interaction_continuous_adj = bc360_int$p_adjusted)
bc360_continuous_results

write_table(bc360_continuous_results, file.path(output_dir, "BC360_module_treatment_interaction.csv"))

###### BC360 APM: main groups and functional subgroups ######

# Validate the 21-gene BC360 APM panel before subgroup annotation.

stopifnot(
  is.character(BC360_apm_genes),
  !anyNA(BC360_apm_genes),
  length(BC360_apm_genes) == 21L,
  !anyDuplicated(BC360_apm_genes))

# Same functional divisions as custom APM80, using only BC360 genes.
bc360_apm_modules = list(

  `MHC-I` = list(
    `Surface MHC-I chains` =
      c("HLA-A", "HLA-B", "HLA-C"),

    `Peptide transport and loading` =
      c("TAP1", "TAP2", "TAPBP"),

    `Immunoproteasome catalytic subunits` =
      c("PSMB9", "PSMB10"),

    `Constitutive proteasome catalytic subunits` =
      c("PSMB7")),

  `MHC-II` = list(
    `HLA-DR, -DQ and -DP chains` =
      c("HLA-DRA", "HLA-DRB1",
        "HLA-DQA1", "HLA-DQB1",
        "HLA-DPA1", "HLA-DPB1"),

    `HLA-DM / HLA-DO peptide editing` =
      c("HLA-DMA", "HLA-DMB", "HLA-DOB")))

# These three BC360 genes are outside the custom APM80 hierarchy.
# Assign these three genes to separate main groups and single-gene subgroups.
bc360_apm_modules[["Lipid antigen processing"]] = list(`CD1E-associated lipid antigen processing` = c("CD1E"))
bc360_apm_modules[["T-cell co-receptor"]] = list(`CD8 alpha co-receptor` = c("CD8A"))
bc360_apm_modules[["Extracellular immune regulation"]] = list(`Thrombospondin-1` = c("THBS1"))

# Derive the five BC360 main-group gene lists.
BC360_APM_MAIN_GROUPS = lapply(bc360_apm_modules, unlist, use.names = FALSE)

# Stable, panel-specific IDs for all nine individual subgroups.
bc360_subgroup_prefix = c(
  `MHC-I` = "bc360_core_mhc_i_genes",
  `MHC-II` = "bc360_core_mhc_ii_genes",
  `Lipid antigen processing` = "bc360_lipid_antigen_processing_genes",
  `T-cell co-receptor` = "bc360_t_cell_coreceptor_genes",
  `Extracellular immune regulation` = "bc360_extracellular_regulation_genes")

stopifnot(identical(names(bc360_apm_modules), names(bc360_subgroup_prefix)))

BC360_APM_SUBGROUPS = do.call(c, lapply(names(bc360_apm_modules), function(module) {
  groups = bc360_apm_modules[[module]]
  names(groups) = paste0(bc360_subgroup_prefix[[module]], "_", seq_along(groups))
  groups}))

# Labels and gene lists are derived from the same definitions.
BC360_APM_SUBGROUP_CATALOG = do.call(rbind,
                                     lapply(names(bc360_apm_modules), function(module) {
                                       groups = bc360_apm_modules[[module]]
                                       data.frame(
                                         module = module,
                                         module_label = module,
                                         subgroup_number = seq_along(groups),
                                         subgroup_id = paste0(bc360_subgroup_prefix[[module]], "_", seq_along(groups)),
                                         subgroup_label = names(groups),
                                         n_genes = lengths(groups),
                                         genes = vapply(groups, paste, collapse = ";", FUN.VALUE = character(1)),
                                         stringsAsFactors = FALSE)}))
rownames(BC360_APM_SUBGROUP_CATALOG) = NULL

# Every original BC360 APM gene must occur exactly once.
bc360_grouped_genes = unlist(BC360_APM_MAIN_GROUPS, use.names = FALSE)
bc360_subgroup_genes = unlist(BC360_APM_SUBGROUPS, use.names = FALSE)

stopifnot(
  length(BC360_APM_MAIN_GROUPS) == 5L,
  identical(unname(lengths(BC360_APM_MAIN_GROUPS)), c(9L, 9L, 1L, 1L, 1L)),
  identical(unname(lengths(bc360_apm_modules)), c(4L, 2L, 1L, 1L, 1L)),
  length(BC360_APM_SUBGROUPS) == 9L,
  nrow(BC360_APM_SUBGROUP_CATALOG) == 9L,
  !anyDuplicated(names(BC360_APM_SUBGROUPS)),
  identical(names(BC360_APM_SUBGROUPS), BC360_APM_SUBGROUP_CATALOG$subgroup_id),
  length(BC360_apm_genes) == 21L,
  !anyDuplicated(BC360_apm_genes),
  length(bc360_grouped_genes) == 21L,
  !anyDuplicated(bc360_grouped_genes),
  identical(bc360_grouped_genes, bc360_subgroup_genes),
  setequal(bc360_grouped_genes, BC360_apm_genes))

# Use BC360_APM_MAIN_GROUPS for the five broad groups;
# use BC360_APM_SUBGROUPS for the nine individual functional subgroups.
# Display names and IDs:
# BC360_APM_SUBGROUP_CATALOG[, c("module_label", "subgroup_label", "subgroup_id", "n_genes")]

###### Antigen processing and presentation machinery (APM) - 80-gene Comprehensive Gene Set ######

# The APM80 panel comprises 80 genes organized into seven modules and 28 submodules.
# The complete gene set, module/submodule definitions and all gene memberships were predefined independently of any clinical outcome.

# Descriptive subgroup names are the source of truth for figure labels.
custom_apm_modules = list(

  `MHC-I` = list(
    `Surface MHC-I chains and B2M` =
      c("HLA-A", "HLA-B", "HLA-C", "B2M"),
    # TAPBPL is a TAPBP paralogue; functions are not interchangeable.
    `Peptide transport and loading` =
      c("TAP1", "TAP2", "TAPBP", "TAPBPL"),
    `Immunoproteasome catalytic subunits` =
      c("PSMB8", "PSMB9", "PSMB10"),
    `Constitutive proteasome catalytic subunits` =
      c("PSMB5", "PSMB6", "PSMB7"),
    `Proteasome activators` =
      c("PSME1", "PSME2"),
    `ER peptide trimming` =
      c("ERAP1", "ERAP2"),
    `ER folding and loading-associated proteins` =
      c("CALR", "CANX", "PDIA3")),

  `MHC-II` = list(
    `HLA-DR, -DQ and -DP chains` =
      c("HLA-DRA", "HLA-DRB1", "HLA-DRB5",
        "HLA-DQA1", "HLA-DQA2", "HLA-DQB1", "HLA-DQB2",
        "HLA-DPA1", "HLA-DPB1"),
    # HLA-DM promotes CLIP exchange; HLA-DO modulates HLA-DM activity.
    `HLA-DM / HLA-DO peptide editing` =
      c("HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB"),
    `Invariant chain: assembly and trafficking` =
      c("CD74"),
    # Proteases can generate or destroy epitopes; no score direction is assumed.
    # LGMN (asparagine endopeptidase) is included in this functional subgroup.
    `Endolysosomal antigen-processing proteases` =
      c("CTSS", "CTSL", "CTSB", "CTSD", "LGMN"),
    # IFI30 also participates in MHC-I cross-presentation. The single primary
    # assignment avoids double-counting its membership in the APM80 panel.
    `Lysosomal thiol reductase` =
      c("IFI30")),

  # Internal key: Surface control; display label: Post-transcriptional MHC-I regulation.
  `Surface control` = list(
    # Includes ER-retained MHC-I capture; this is not exclusively surface turnover.
    `Autophagy-associated MHC-I degradation` =
      c("NBR1"),
    `Endocytic MHC-I turnover` =
      c("MAL2"),
    # Accessibility is distinct from HLA-I abundance.
    `Glycosylation-related HLA-I accessibility` =
      c("SPPL3")),

  `Cross presentation` = list(
    `Candidate trafficking and redox factors` =
      c("LNPEP", "RAB7A", "SEC22B", "CYBB"),
    `Candidate antigen export / dislocation factors` =
      c("DERL1", "DERL2", "VCP", "SEC61A1")),

  `Regulators` = list(
    `MHC-I transcriptional coactivator` =
      c("NLRC5"),
    `MHC-II transcriptional coactivator` =
      c("CIITA"),
    `RFX DNA-binding complex` =
      c("RFX5", "RFXANK", "RFXAP"),
    # Retained regardless of abundance. Variation/estimability is assessed separately; the ligand transcript is not an IFN activity signature.
    `Interferon-gamma ligand transcript` =
      c("IFNG"),
    # Type-I-associated set; JAK1 and STAT1 also participate in IFN-gamma signalling, and TYK2 has other cytokine roles.
    `Focused type-I interferon signalling machinery` =
      c("JAK1", "TYK2", "STAT1", "STAT2", "IRF9"),
    # Receptor/JAK2 components only: this is not the complete IFN-gamma pathway.
    `Interferon-gamma receptor/JAK2 components` =
      c("IFNGR1", "IFNGR2", "JAK2"),
    # IRF1 and IRF2 have context-dependent regulatory effects; no fixed activating/repressive direction or pathway activation is inferred from RNA.
    `IRF1/IRF2 transcriptional regulators` =
      c("IRF1", "IRF2")),

  `Nonclassical evasion` = list(
    `Nonclassical HLA-I molecules` =
      c("HLA-E", "HLA-G", "HLA-F")),

  `Repressors` = list(
    `PRC2 chromatin regulators` =
      c("EZH2", "SUZ12", "EED"),
    `DNA methyltransferases` =
      c("DNMT1", "DNMT3A", "DNMT3B"),
    `Histone deacetylases` =
      c("HDAC1", "HDAC2", "HDAC3")))

# Stable internal module keys; use these labels in tables and figures.
APM80_MODULE_LABELS = setNames(names(custom_apm_modules), names(custom_apm_modules))
APM80_MODULE_LABELS[["Surface control"]] = "Post-transcriptional MHC-I regulation"

# Derived views retain the object names and legacy prefixes used by the analyses.
FIXED_APM80 = lapply(custom_apm_modules, unlist, use.names = FALSE)

legacy_prefix = c(
  `MHC-I` = "apm_core_mhc_i_genes",
  `MHC-II` = "apm_core_mhc_ii_genes",
  `Surface control` = "apm_surface_control_genes",
  `Cross presentation` = "apm_cross_presentation_genes",
  `Regulators` = "apm_regulators_genes",
  `Nonclassical evasion` = "apm_nonclassical_evasion_genes",
  `Repressors` = "apm_repressors_genes")
stopifnot(identical(names(custom_apm_modules), names(legacy_prefix)))

FIXED_APM80_SUBGROUPS = do.call(c, lapply(names(custom_apm_modules), function(module) {
  groups = custom_apm_modules[[module]]
  prefix = legacy_prefix[[module]]
  names(groups) = if (module == "Nonclassical evasion") {
    prefix
  } else {
    paste0(prefix, "_", seq_along(groups))
  }
  groups
}))

panel_genes = unlist(FIXED_APM80, use.names = FALSE)
subgroup_genes = unlist(FIXED_APM80_SUBGROUPS, use.names = FALSE)

# Display/annotation catalogue derived from the same gene definitions.
APM80_SUBGROUP_CATALOG = do.call(rbind, lapply(names(custom_apm_modules), function(module) {
  groups = custom_apm_modules[[module]]
  prefix = legacy_prefix[[module]]
  data.frame(
    module = module,
    module_label = unname(APM80_MODULE_LABELS[[module]]),
    subgroup_number = seq_along(groups),
    subgroup_id = if (module == "Nonclassical evasion") prefix else
      paste0(prefix, "_", seq_along(groups)),
    subgroup_label = names(groups),
    n_genes = lengths(groups),
    genes = vapply(groups, paste, collapse = ";", FUN.VALUE = character(1)),
    stringsAsFactors = FALSE)}))
rownames(APM80_SUBGROUP_CATALOG) = NULL

stopifnot(
  length(panel_genes) == 80L,
  !anyDuplicated(panel_genes),
  identical(panel_genes, subgroup_genes),
  length(FIXED_APM80) == 7L,
  identical(unname(lengths(FIXED_APM80)), c(21L, 20L, 3L, 8L, 16L, 3L, 9L)),
  identical(unname(lengths(custom_apm_modules)), c(7L, 5L, 3L, 2L, 7L, 1L, 3L)),
  length(FIXED_APM80_SUBGROUPS) == 28L,
  nrow(APM80_SUBGROUP_CATALOG) == 28L,
  sum(lengths(FIXED_APM80_SUBGROUPS) == 1L) == 8L,
  "LGMN" %in% FIXED_APM80_SUBGROUPS$apm_core_mhc_ii_genes_4,
  identical(FIXED_APM80_SUBGROUPS$apm_core_mhc_ii_genes_5, "IFI30"),
  !"apm_core_mhc_ii_genes_6" %in% names(FIXED_APM80_SUBGROUPS),
  all(paste0("apm_regulators_genes_", 5:7) %in% names(FIXED_APM80_SUBGROUPS)),
  !"apm_regulators_genes_8" %in% names(FIXED_APM80_SUBGROUPS))

# Add the custom APM80 panel to the gene-set registry.
if (!exists("apm_gene_sets", inherits = FALSE)) apm_gene_sets = list()
stopifnot(is.list(apm_gene_sets))
apm_gene_sets[["APM_CUSTOM_SET"]] = panel_genes

# Verify consistency between module definitions, subgroup lists and the complete panel.
stopifnot(
  is.list(custom_apm_modules), is.list(FIXED_APM80), is.list(FIXED_APM80_SUBGROUPS),
  identical(names(custom_apm_modules), names(FIXED_APM80)),
  identical(FIXED_APM80, lapply(custom_apm_modules, unlist, use.names = FALSE)),
  identical(panel_genes, unlist(FIXED_APM80, use.names = FALSE)),
  identical(panel_genes, unlist(FIXED_APM80_SUBGROUPS, use.names = FALSE)),
  length(FIXED_APM80) == 7L, length(FIXED_APM80_SUBGROUPS) == 28L,
  length(panel_genes) == 80L, !anyDuplicated(panel_genes))

###### APM gene-set coverage ######

# Report the number and proportion of genes measured in PANTHER for each APM panel.
apm_coverage = data.frame(
  gene_set = names(apm_gene_sets),
  n_defined = lengths(apm_gene_sets),
  n_present_in_PANTHER = vapply(apm_gene_sets, function(g) length(intersect(unique(g), rownames(tmm_normalized_tpm_symbols))), integer(1)))

apm_coverage$coverage = with(apm_coverage, round(n_present_in_PANTHER / n_defined, 3))
print(apm_coverage, row.names = FALSE)

###### APM signature: definition (HLA6 - IFN5), PANTHER score and APM groups ######

#   Part 1  Select the final signature using BCRFS treatment-interaction tests of BC360 anchors and eligible predefined APM80 submodules.
#   Part 2  Scoring definition (CONFIG), used again for SCAN-B, MATADOR and TEX.
#   Part 3  PANTHER score on all RNA-matched patients; adds APM_score and APM_group to all_clinical.
#
# Needs: BC360_APM_SUBGROUPS / _CATALOG, FIXED_APM80_SUBGROUPS, APM80_SUBGROUP_CATALOG, panel_genes, adjust_terms, output_dir, all_clinical and tmm_normalized_tpm_symbols. 
# The scoring functions (score_signature_GSVA(), resolve_signature_pair(), score_signature(), ...) are in the Functions section.

###### Part 1. Signature selection using predefined gene sets ######

# The APM80 panel and all module/submodule gene memberships were predefined independently of outcomes and remain fixed during signature selection.
# The final signature components and score direction are selected using adjusted BCRFS treatment-by-score interaction tests in PANTHER.

# BCRFS; Cox model with adjust_terms; continuous treatment x score interaction tested
# with a 1-df likelihood-ratio test; the same complete cases throughout.
#
#   Step 1  Anchor gate. GSVA-score the multi-gene BC360 APM subgroups and test each one. Subgroups with interaction p < 0.05 become anchors. 
#           If none passes, no signature is defined.
#   Step 2  Partners. Each anchor is MHC-I or MHC-II (its BC360 main group). 
#           Its candidate partners are the APM80 subgroups with >= 2 genes, sharing no gene with the anchor, from the main groups that act directly on that class:
#              MHC-I  <- Regulators, Repressors, Post-transcriptional MHC-I regulation
#              MHC-II <- Regulators, Repressors
#   Step 3  Combination. score = z(GSVA(anchor) + s * GSVA(partner)), s = +1 or -1.
#           The combination with the largest interaction chi2 over all anchors, partners and both directions is the best combination.
#           Partner gate (SIGDEF_PARTNER_GATE = TRUE): the anchors alone also compete, so the partner is kept only if the best combination has a smaller interaction p than every passing anchor alone (ties: no partner).
#
# The anchor becomes the positive gene set and the partner, if kept, the negative gene set;
# The same three steps, re-run on treatment labels shuffled within ER strata, give the whole-pipeline permutation p in workflows/permutation_test/scripts/PANTHER_APM_score_Discovery_Permutation.R;
# Part 1 reads its summary.

SIGDEF_ENDPOINT = "BCRFS"
SIGDEF_TREAT_LEVELS = c("FEC-D", "tddEC-D")   # first = reference
SIGDEF_GATE_P = 0.05
SIGDEF_MIN_GENES = 2L
SIGDEF_PARTNER_GATE = TRUE                    # keep a partner only if it lowers the anchors' interaction p
SIGDEF_DIRECT_ACTORS = list(
  `MHC-I`  = c("Regulators", "Repressors", "Post-transcriptional MHC-I regulation"),
  `MHC-II` = c("Regulators", "Repressors"))
SIGDEF_OUTPUT_DIR = file.path(output_dir, "APM_signature_definition")

# Short labels (display only); identical to those in the permutation script, so the selected combination can be matched to the permutation output.
SIGDEF_SHORT = c(
  bc360_core_mhc_i_genes_1 = "MHC-I chains", bc360_core_mhc_i_genes_2 = "TAP",
  bc360_core_mhc_i_genes_3 = "PSMB9/10",     bc360_core_mhc_ii_genes_1 = "HLA6",
  bc360_core_mhc_ii_genes_2 = "DM/DO",
  apm_regulators_genes_3 = "RFX",  apm_regulators_genes_5 = "IFN5",
  apm_regulators_genes_6 = "IFNGR/JAK2", apm_regulators_genes_7 = "IRF1/IRF2",
  apm_repressors_genes_1 = "PRC2", apm_repressors_genes_2 = "DNMT", apm_repressors_genes_3 = "HDAC")

# Whole-pipeline permutation (permutation script; ER-stratified, 5,000 permutations).
# TRUE (default) runs the permutation script here (5,000 permutations); FALSE reads a saved summary from SIGDEF_PERMUTATION_DIR if there is one.
SIGDEF_RUN_PERMUTATION = TRUE
SIGDEF_PERMUTATION_SCRIPT = file.path(code_dir, "workflows", "permutation_test", "scripts", "PANTHER_APM_score_Discovery_Permutation.R")
SIGDEF_PERMUTATION_DIR = file.path(output_dir, "APM_discovery_permutation")

APM_SIGNATURE_DEFINITION = local({

  for (obj in c("BC360_APM_SUBGROUPS", "BC360_APM_SUBGROUP_CATALOG", "FIXED_APM80_SUBGROUPS", "APM80_SUBGROUP_CATALOG", "adjust_terms", "all_clinical", "tmm_normalized_tpm_symbols"))
    if (!exists(obj)) stop(obj, " not found: run the gene-set definitions above first.")
  
  if (!exists("score_signature_GSVA", mode = "function")) stop("score_signature_GSVA() not found.")
  
  dir.create(SIGDEF_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  short = function(id) unname(ifelse(id %in% names(SIGDEF_SHORT), SIGDEF_SHORT[id], id))

  ## Anchors: the multi-gene BC360 APM subgroups
  bc = BC360_APM_SUBGROUP_CATALOG[BC360_APM_SUBGROUP_CATALOG$n_genes >= SIGDEF_MIN_GENES, , drop = FALSE]
  
  stopifnot(nrow(bc) == 5L, all(bc$module %in% names(SIGDEF_DIRECT_ACTORS)), all(bc$subgroup_id %in% names(BC360_APM_SUBGROUPS)))
  
  anchors = data.frame(id = bc$subgroup_id, short = short(bc$subgroup_id), class = bc$module,
                       label = bc$subgroup_label, n_genes = bc$n_genes, genes = bc$genes, stringsAsFactors = FALSE)

  ## Step 2 pools (outcome-blind; fixed before any model is fitted)
  a80 = APM80_SUBGROUP_CATALOG
  
  stopifnot(all(unlist(SIGDEF_DIRECT_ACTORS) %in% a80$module_label))
  
  candidates = do.call(rbind, lapply(seq_len(nrow(anchors)), function(k) {
    pool = a80[a80$module_label %in% SIGDEF_DIRECT_ACTORS[[anchors$class[k]]], , drop = FALSE]
    
    shares = vapply(pool$subgroup_id, function(id)
      length(intersect(FIXED_APM80_SUBGROUPS[[id]], BC360_APM_SUBGROUPS[[anchors$id[k]]])) > 0L, logical(1))
    
    data.frame(anchor_id = anchors$id[k], anchor = anchors$short[k], anchor_class = anchors$class[k],
               partner_id = pool$subgroup_id, partner = short(pool$subgroup_id),
               main_group = pool$module_label, subgroup_label = pool$subgroup_label,
               n_genes = pool$n_genes, genes = pool$genes,
               candidate = pool$n_genes >= SIGDEF_MIN_GENES & !shares,
               reason = ifelse(pool$n_genes < SIGDEF_MIN_GENES, "single gene", ifelse(shares, "shares genes with anchor", "candidate")),
               stringsAsFactors = FALSE, row.names = NULL)
  }))

  ## GSVA scores: all RNA-matched patients, cohort background (finite, non-constant
  ## genes), expression only and computed once
  ids = intersect(as.character(all_clinical$sample_id), colnames(tmm_normalized_tpm_symbols))
  X = as.matrix(tmm_normalized_tpm_symbols[, ids, drop = FALSE])
  X = X[rowSums(!is.finite(X)) == 0, , drop = FALSE]
  X = X[apply(X, 1L, stats::sd) > 1e-8, , drop = FALSE]
  partner_ids = unique(candidates$partner_id[candidates$candidate])
  sets = c(BC360_APM_SUBGROUPS[anchors$id], FIXED_APM80_SUBGROUPS[partner_ids])
  absent = lapply(sets, setdiff, rownames(X))
  
  if (any(lengths(absent)))
    warning("Genes missing or constant in the GSVA background: ", paste(names(absent)[lengths(absent) > 0], collapse = ", "))
  
  if (any(lengths(sets) - lengths(absent) < 2L)) stop("A set has fewer than 2 scorable genes.")
  gs = score_signature_GSVA(X, sets, suffix = "", min_size = 1)
  S = as.matrix(gs[match(ids, gs$sample_id), names(sets), drop = FALSE])
  rownames(S) = ids
  stopifnot(!anyNA(S))

  ## Analysis set: complete cases for the adjusted model
  cl = as.data.frame(all_clinical)
  cl$sample_id = as.character(cl$sample_id)
  cl = cl[match(ids, cl$sample_id), , drop = FALSE]
  time_col = paste0("time", SIGDEF_ENDPOINT)
  status_col = paste0(SIGDEF_ENDPOINT, ".status")
  adj_vars = all.vars(stats::as.formula(paste("~", adjust_terms)))
  need = unique(c("sample_id", time_col, status_col, "TREAT", adj_vars))
  if (length(setdiff(need, names(cl)))) stop("Not in all_clinical: ", paste(setdiff(need, names(cl)), collapse = ", "))
  D = cl[stats::complete.cases(cl[, need]), need, drop = FALSE]
  D$TREAT = factor(as.character(D$TREAT), levels = SIGDEF_TREAT_LEVELS)
  D$.time = D[[time_col]]
  D$.status = D[[status_col]]
  stopifnot(!anyNA(D$TREAT), all(D$.status %in% c(0, 1)))
  rows = match(D$sample_id, ids)   # rows of S for the analysis patients

  F0 = stats::as.formula(paste("Surv(.time, .status) ~ TREAT + SCORE +", adjust_terms))
  F1 = stats::as.formula(paste("Surv(.time, .status) ~ TREAT * SCORE +", adjust_terms))
  INT = paste0("TREAT", SIGDEF_TREAT_LEVELS[2], ":SCORE")
  zs = function(x) (x - mean(x)) / stats::sd(x)   # over all RNA-matched patients

  int_test = function(score) {
    D$SCORE = score[rows]
    f0 = coxph(F0, data = D)
    f1 = coxph(F1, data = D)
    stopifnot(f0$n == f1$n, f0$nevent == f1$nevent)
    b = stats::coef(f1)
    V = stats::vcov(f1)
    chi2 = max(0, 2 * (f1$loglik[2] - f0$loglik[2]))
    slope = function(w) {
      est = sum(w * b[names(w)])
      se = sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w))
      exp(est + c(0, -1.96, 1.96) * se)
    }
    i = slope(stats::setNames(1, INT))
    c(chi2 = chi2, p = stats::pchisq(chi2, 1, lower.tail = FALSE),
      int_HR_per_SD = i[1], int_low95 = i[2], int_high95 = i[3],
      HR_per_SD_reference_arm = slope(c(SCORE = 1))[1],
      HR_per_SD_comparison_arm = slope(stats::setNames(c(1, 1), c("SCORE", INT)))[1])
  }

  ## Step 1: gate
  gate = cbind(anchors[, c("short", "id", "class", "label", "n_genes", "genes")],
               t(vapply(anchors$id, function(id) int_test(zs(S[, id])), numeric(7))), row.names = NULL)
  gate$passes_gate = gate$p < SIGDEF_GATE_P
  if (!any(gate$passes_gate))
    stop("No BC360 subgroup has interaction p < ", SIGDEF_GATE_P, ": the pipeline defines no signature.")

  ## Steps 2-3: candidates of the passing anchors, both directions
  grid = do.call(rbind, lapply(gate$id[gate$passes_gate], function(a) {
    p = candidates[candidates$anchor_id == a & candidates$candidate, , drop = FALSE]
    do.call(rbind, lapply(c(1, -1), function(s)
      data.frame(p[, c("anchor_id", "anchor", "partner_id", "partner", "main_group", "subgroup_label")], sign = s)))
  }))
  grid = cbind(grid, t(vapply(seq_len(nrow(grid)), function(j)
    int_test(zs(S[, grid$anchor_id[j]] + grid$sign[j] * S[, grid$partner_id[j]])), numeric(7))))
  grid$combination = paste(grid$anchor, ifelse(grid$sign > 0, "+", "-"), grid$partner)
  grid = grid[order(-grid$chi2), , drop = FALSE]   # stable: ties keep the script-3 order
  rownames(grid) = NULL
  best_combo = grid[1L, ]

  ## Partner gate: the best anchor alone competes with the best combination
  passed = gate[gate$passes_gate, , drop = FALSE]
  best_anchor = passed[which.max(passed$chi2), ]
  partner_kept = !SIGDEF_PARTNER_GATE || best_combo$chi2 > best_anchor$chi2
  grid$selected = partner_kept & seq_len(nrow(grid)) == 1L
  sel = if (partner_kept) best_combo else
    data.frame(best_anchor[, c("chi2", "p", "int_HR_per_SD", "int_low95", "int_high95",
                               "HR_per_SD_reference_arm", "HR_per_SD_comparison_arm")],
               anchor_id = best_anchor$id, anchor = best_anchor$short, partner_id = NA_character_, partner = "",
               sign = NA_real_, combination = best_anchor$short, stringsAsFactors = FALSE)

  if (partner_kept && sel$sign > 0)
    stop("Selected: ", sel$combination, ". The signature helpers score positive minus negative only; ",
         "a summed score needs its own method before this selection can be carried forward.")

  ## Whole-pipeline permutation (permutation script)
  if (SIGDEF_RUN_PERMUTATION) {
    # Fresh run: move the saved replicates aside, otherwise the permutation script resumes from them
    old = file.path(SIGDEF_PERMUTATION_DIR, c("05_permutation_replicates.csv", "06_permutation_summary.csv"))
    old = old[file.exists(old)]
    file.rename(old, sub("\\.csv$", "_previous.csv", old))
    perm_env = new.env(parent = globalenv())
    perm_env$OUT_DIR_OVERRIDE = SIGDEF_PERMUTATION_DIR   # the permutation script writes its output here
    source(SIGDEF_PERMUTATION_SCRIPT, local = perm_env)
  }
  perm_file = file.path(SIGDEF_PERMUTATION_DIR, "06_permutation_summary.csv")
  perm = if (file.exists(perm_file)) utils::read.csv(perm_file, stringsAsFactors = FALSE) else NULL
  if (is.null(perm)) {
    message("No permutation summary at ", perm_file, ". Run the permutation script or set SIGDEF_RUN_PERMUTATION = TRUE.")
  } else if (!identical(perm$observed_selected, sel$combination) || abs(perm$observed_chi2 - sel$chi2) > 1e-3) {
    warning("The permutation script selected ", perm$observed_selected, " (chi2 ", round(perm$observed_chi2, 2),
            "); this section selected ", sel$combination, " (chi2 ", round(sel$chi2, 2),
            "). Check that both used the same data and settings.")
  }

  selected = data.frame(
    signature = sel$combination, endpoint = SIGDEF_ENDPOINT, n = nrow(D), events = sum(D$.status),
    positive_set = sel$anchor, positive_genes = paste(BC360_APM_SUBGROUPS[[sel$anchor_id]], collapse = ";"),
    negative_set = sel$partner,
    negative_genes = if (partner_kept) paste(FIXED_APM80_SUBGROUPS[[sel$partner_id]], collapse = ";") else "",
    partner_gate = SIGDEF_PARTNER_GATE, partner_kept = partner_kept,
    best_combination = best_combo$combination, best_combination_chi2 = best_combo$chi2,
    best_anchor_alone = best_anchor$short, best_anchor_alone_chi2 = best_anchor$chi2,
    anchor_alone_chi2 = gate$chi2[gate$id == sel$anchor_id], anchor_alone_p = gate$p[gate$id == sel$anchor_id],
    chi2 = sel$chi2, p_nominal = sel$p,
    int_HR_per_SD = sel$int_HR_per_SD, int_low95 = sel$int_low95, int_high95 = sel$int_high95,
    HR_per_SD_FEC_D = sel$HR_per_SD_reference_arm, HR_per_SD_tddEC_D = sel$HR_per_SD_comparison_arm,
    n_combinations = nrow(grid),
    n_permutations = if (is.null(perm)) NA else perm$n_permutations,
    perm_strata = if (is.null(perm)) NA else perm$perm_strata,
    whole_pipeline_p = if (is.null(perm)) NA else perm$whole_pipeline_p,
    whole_pipeline_p_95ci = if (is.null(perm)) NA else
      sprintf("%.3f-%.3f", perm$mc_95ci_low, perm$mc_95ci_high),
    stringsAsFactors = FALSE)

  utils::write.csv(gate, file.path(SIGDEF_OUTPUT_DIR, "01_step1_BC360_gate.csv"), row.names = FALSE)
  utils::write.csv(candidates, file.path(SIGDEF_OUTPUT_DIR, "02_step2_partner_candidates.csv"), row.names = FALSE)
  utils::write.csv(grid, file.path(SIGDEF_OUTPUT_DIR, "03_step3_combinations.csv"), row.names = FALSE)
  utils::write.csv(selected, file.path(SIGDEF_OUTPUT_DIR, "04_selected_signature.csv"), row.names = FALSE)

  message(sprintf(paste0(
    "Signature definition (%s, n = %d, %d events)\n",
    "  Anchor gate passed: %s\n",
    "  Partner gate: best of %d combinations %s (chi2 %.2f) vs best anchor alone %s (chi2 %.2f): %s\n",
    "  Selected: %s; chi2 %.2f, nominal p %.2g; interaction HR per SD %.2f (%.2f-%.2f)\n",
    "  Whole-pipeline p: %s"),
    SIGDEF_ENDPOINT, nrow(D), sum(D$.status),
    paste(sprintf("%s (p %s)", gate$short[gate$passes_gate], formatC(gate$p[gate$passes_gate], digits = 2, format = "g")), collapse = ", "),
    nrow(grid), best_combo$combination, best_combo$chi2, best_anchor$short, best_anchor$chi2,
    if (!SIGDEF_PARTNER_GATE) "gate off" else if (partner_kept) "partner kept" else "no partner",
    sel$combination, sel$chi2, sel$p, sel$int_HR_per_SD, sel$int_low95, sel$int_high95,
    if (is.null(perm)) "not available" else
      sprintf("%.3f (%d permutations, strata %s)", perm$whole_pipeline_p, perm$n_permutations, perm$perm_strata)))

  list(positive_name = make.names(sel$anchor), positive_genes = unname(BC360_APM_SUBGROUPS[[sel$anchor_id]]),
       negative_name = if (partner_kept) make.names(sel$partner) else NULL,
       negative_genes = if (partner_kept) unname(FIXED_APM80_SUBGROUPS[[sel$partner_id]]) else character(0),
       positive_id = sel$anchor_id, negative_id = if (partner_kept) sel$partner_id else NA_character_,
       signature = sel$combination, partner_kept = partner_kept,
       score = stats::setNames(if (partner_kept) zs(S[, sel$anchor_id] + sel$sign * S[, sel$partner_id])
                               else zs(S[, sel$anchor_id]), ids),   # cross-check for Part 3
       gate = gate, candidates = candidates, combinations = grid, selected = selected,
       permutation = perm, output_dir = SIGDEF_OUTPUT_DIR)
})

# The selected gene sets, by name, for CONFIG (here: HLA6 and IFN5). Without a partner,
# negative_name is NULL and Part 3 scores the anchor alone (positive-only mode).
assign(APM_SIGNATURE_DEFINITION$positive_name, APM_SIGNATURE_DEFINITION$positive_genes)

if (APM_SIGNATURE_DEFINITION$partner_kept) assign(APM_SIGNATURE_DEFINITION$negative_name, APM_SIGNATURE_DEFINITION$negative_genes)

###### Part 2. Scoring definition (shared with SCAN-B, MATADOR and TEX) ######

# APM_score = z[GSVA(positive) - GSVA(negative)] with the gene sets selected in Part 1 (score_signature(), Functions section).
# Each ID is an APM80 subgroup label/ID or the name of a gene vector in the session (here the HLA6 and IFN5 vectors assigned above);
# negative NULL = positive-only score.
CONFIG = list(
  positive_subgroup_id = APM_SIGNATURE_DEFINITION$positive_name,   # "HLA6"
  negative_subgroup_id = APM_SIGNATURE_DEFINITION$negative_name)   # "IFN5"; NULL if no partner

###### Part 3. PANTHER score and APM groups -> APM_score, APM_group ######

# All RNA-matched PANTHER patients (the same patients as Part 1), before any endpoint
# or covariate exclusion. GSVA uses the whole expression matrix as background; the
# score is standardized and split at its median (High >= median) over these patients.

PANTHER_APM_SCORES = local({
  pair = resolve_signature_pair(CONFIG)
  
  ids = intersect(as.character(all_clinical$sample_id), colnames(tmm_normalized_tpm_symbols))
  e = as.matrix(tmm_normalized_tpm_symbols[, ids, drop = FALSE])
  arm = as.character(all_clinical$TREAT[match(ids, as.character(all_clinical$sample_id))])
  
  score = score_signature(e, pair$plus, pair$minus, "PANTHER APM score")
  out = data.frame(sample_id = ids, TREAT = arm, APM_score = unname(score), APM_group = apm_group(score), stringsAsFactors = FALSE)
  
  # Same score as the one Part 1 selected (same patients, GSVA background and standardization)
  diff = max(abs(out$APM_score - APM_SIGNATURE_DEFINITION$score[out$sample_id]))
  if (!is.finite(diff) || diff > 1e-6)
    warning("APM_score differs from the Part 1 score (max |difference| ", signif(diff, 3), ").")
  
  utils::write.csv(out, file.path(SIGDEF_OUTPUT_DIR, "05_PANTHER_APM_score_and_group.csv"), row.names = FALSE)
  
  message(sprintf("APM_score: %s, signed GSVA; %d RNA-matched patients (APM-high %d, APM-low %d)",
                  pair$signature_label, nrow(out), sum(out$APM_group == "APM-high"), sum(out$APM_group == "APM-low")))
  
  attr(out, "signature_label") = pair$signature_label
  attr(out, "method") = "signed GSVA"
  out
})

# Add to all_clinical by sample ID; every clinical row is kept, patients without
# RNA get NA. Rerunning replaces both columns.
apm_rows = match(as.character(all_clinical$sample_id), PANTHER_APM_SCORES$sample_id)
all_clinical$APM_score = PANTHER_APM_SCORES$APM_score[apm_rows]
all_clinical$APM_group = PANTHER_APM_SCORES$APM_group[apm_rows]
stopifnot(!anyNA(all_clinical$APM_group[!is.na(all_clinical$APM_score)]), sum(!is.na(all_clinical$APM_score)) == nrow(PANTHER_APM_SCORES))
rm(apm_rows)

# Also add the APM columns to the ER+/HER2- subset created in the SET2,3 section
m = match(all_clinical_ERpos_HERneg$sample_id, all_clinical$sample_id)
all_clinical_ERpos_HERneg$APM_score = all_clinical$APM_score[m]
all_clinical_ERpos_HERneg$APM_group = all_clinical$APM_group[m]

###### ----- Table 1 and Supplementary Table 2 -----
### Tables with Clinicopathologic Characteristics ###

biomarker_selected = 'APM_group'

if (nzchar(chrome_path)) Sys.setenv(CHROMOTE_CHROME = chrome_path)
# On Windows chromote looks for Chrome only in the machine-wide registry; also look for a per-user Chrome and for Microsoft Edge
if (.Platform$OS.type == "windows" && !nzchar(Sys.getenv("CHROMOTE_CHROME"))) {
  app_path = function(exe, hive) tryCatch(gsub('"', "", utils::readRegistry(paste0("SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\App Paths\\", exe, "\\"), hive = hive)[["(Default)"]]), error = function(e) NA_character_)
  browsers = c(app_path("chrome.exe", "HLM"), app_path("chrome.exe", "HCU"), app_path("msedge.exe", "HLM"), app_path("msedge.exe", "HCU"))
  browsers = browsers[!is.na(browsers) & file.exists(browsers)]
  if (length(browsers)) Sys.setenv(CHROMOTE_CHROME = browsers[1])
}

## Table 1 by Treatment arm (ALL) ##
# Present treatment-stratified characteristics for the full, Swedish and translational PANTHER cohorts in adjacent tables.

# input data frames
df1 = all_clinical_original
df2 = all_clinical_original[all_clinical_original$patno < 12000,]
df3 = all_clinical

# create one table per data frame
tbl1 = make_table1(df1)
tbl2 = make_table1(df2)
tbl3 = make_table1(df3)

# Merge the three cohort tables side by side.
table1_combined = tbl_merge(tbls = list(tbl1, tbl2, tbl3), tab_spanner = c("All PANTHER cohort", "Swedish PANTHER cohort", "Translational PANTHER cohort"))

# convert to gt and save as one png
gt_table1_combined = as_gt(table1_combined) |>
  cols_width( label ~ px(270), everything() ~ px(90)) |>
  tab_options(table.font.size = px(18),             # body
              column_labels.font.size = px(20),     # column headers
              heading.title.font.size = px(24))  |> # title
  tab_header(title = "Table 1: Patient Characteristics by Treatment Group in PANTHER trial")

gtsave(gt_table1_combined, vwidth = 2200, zoom = 2, expand = 10, filename = file.path(output_dir, "Table1_all.png"))

## Supplementary Table 2: by APM group ##
# Summarize and compare baseline characteristics between APM-low and APM-high groups.

table2 = tbl_summary(
  all_clinical,
  include = c(age, mensstat, optype, tumdiam_cat, posnodes_cat, histgrade, Clinical_Subtype, subtypecd_ssp.subtype, all_of(biomarker_selected)),
  by = all_of(biomarker_selected),  # split table by group
  missing = "ifany",                # Show missing-data rows when missing values are present.
  label = list(
    age = "Age (years)",
    mensstat = "Menopausal Status",
    optype = "Operation Type",
    tumdiam_cat = "Tumor Diameter (cm)",
    posnodes_cat = "Positive Nodes",
    histgrade = "Histological Grade",
    subtypecd_ssp.subtype = "Intrinsic Subtype",
    Clinical_Subtype = "Clinical Subtype")) |>
  add_n() |> # add column with total number of non-missing observations
  add_p() |> # test for a difference between groups
  modify_header(label = "**Variable**") |> # update the column header
  bold_labels()

# Convert to gt object
gt_table1 = as_gt(table2) |> tab_header(title = "Supplementary Table 2: Patient Characteristics by APM group in PANTHER")

gtsave(gt_table1, filename = file.path(output_dir, "Supp_Table_2_APM_group.png"))

###### ----- Check PH assumption ----- ######

# 1) ER as a covariate: does ER violate PH? (the reason for stratifying)
# 2) Primary model, ER-stratified: TREAT * APM_score + strata(ERpos) + newHER2 + scale(tumdiam) + scale(posnodes) + histgrade
ph_endpoints = c("BCRFS", "EFS", "DDFS", "OS")
ph_results_er_covariate = lapply(ph_endpoints, check_ph_er_covariate, data = all_clinical)
ph_results = lapply(ph_endpoints, check_ph_er_stratified, data = all_clinical)

ph_table_all_clean = bind_rows(lapply(c(ph_results_er_covariate, ph_results), `[[`, "table")) %>%
  mutate(
    model = factor(model, levels = c("ER as covariate", "ER-stratified")),
    endpoint = factor(endpoint, levels = ph_endpoints),
    term = case_when(
      term == "TREAT" ~ "Treatment arm",
      term == "APM_score" ~ "APM score",
      term == "TREAT:APM_score" ~ "Treatment x APM score",
      term == "ERpos" ~ "ER status",
      term == "newHER2" ~ "HER2 status",
      term == "scale(tumdiam)" ~ "Tumor size",
      term == "scale(posnodes)" ~ "Positive lymph nodes",
      term == "histgrade" ~ "Histological grade",
      term == "GLOBAL" ~ "Global test",
      TRUE ~ term),
    chisq = round(chisq, 3),
    p_print = format_p_regular(p),
    p = signif(p, 3)) %>%
  arrange(model, endpoint)

ph_table_all_clean

# Summarize ER and global PH tests, plus significant terms in ER-stratified models.
ph_table_all_clean %>%
  filter((model == "ER as covariate" & term %in% c("ER status", "Global test")) |
           (model == "ER-stratified" & (term == "Global test" | p < 0.05))) %>%
  select(model, endpoint, term, chisq, df, p_print)

write_table(ph_table_all_clean, file.path(output_dir, "PH_checks_ER_covariate_vs_stratified.csv"))

###### ----- Association Analysis for the continuous APM score -----
###### Continuous APM score - HR (Forest plot) - Per arm ######

# Estimate continuous-score HRs overall and in separate treatment-arm models.
# Obtain treatment-interaction p-values from pooled adjusted Cox models.
sig_cat = 'APM_score'
hr_column = paste(rep(' ', 40), collapse = ' ')

for (endpoint in c('BCRFS', 'EFS', 'DDFS', 'OS')) {
  
  # Title row
  row_cols = c('white', 'white', 'slategray1')
  fonts_3 = c('bold')
  
  results_mv_all_signatures = data.frame('  ' = character(), 'N' = numeric(),
                                         'est' =  numeric(), 'low' = numeric(), 'hi' = numeric(), 'se' = numeric(), ` ` = character(), `HR (95% CI)` = numeric(),
                                         'p value' = numeric(), check.names = FALSE)
  
  results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('  ' = 'Treatment arm', 'N' = '',
                                                                          'est' = NA, 'low' = NA, 'hi' = NA, 'se' = NA, ` ` = hr_column, `HR (95% CI)` = '',
                                                                          'p value' = '', check.names = FALSE))
  
  for (arm in c('FEC-D', 'tddEC-D')) {
    
    all_clinical_arm = all_clinical[all_clinical$TREAT == arm, ]
    
    fonts_3 = c(fonts_3, 'plain')
    
    data_cox = all_clinical_arm[!is.na(all_clinical_arm$ERpos) & !is.na(all_clinical_arm$newHER2) & !is.na(all_clinical_arm$tumdiam) & !is.na(all_clinical_arm$posnodes) & !is.na(all_clinical_arm$histgrade),]
    res.cox = coxph(as.formula(paste0('Surv(time', endpoint, '/365.25, ', endpoint, '.status) ~ ', sig_cat, ' + ', adjust_terms)), data = data_cox)
    cox_model_results = tidy(res.cox, exponentiate = TRUE, conf.int = TRUE)
    
    for (i in seq_along(cox_model_results$term)) {
      term = cox_model_results$term[i]
      
      results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('  ' = term,
                                                                              'N' = if (i == 1) res.cox$n else '',
                                                                              'est' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'estimate']),2),
                                                                              'low' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'conf.low']),2),
                                                                              'hi' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'conf.high']),2),
                                                                              'se' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'std.error']),2),
                                                                              ` ` = hr_column,
                                                                              `HR (95% CI)` = format_hr_ci(as.numeric(cox_model_results[cox_model_results$term == term, 'estimate']), as.numeric(cox_model_results[cox_model_results$term == term, 'conf.low']), as.numeric(cox_model_results[cox_model_results$term == term, 'conf.high'])),
                                                                              'p value' = format_p_regular(as.numeric(cox_model_results[cox_model_results$term == term, 'p.value']), digits = 3), check.names = FALSE))
    }
  }
  results_mv_all_signatures[['  ']][results_mv_all_signatures[['  ']] == "APM_score"][1] = 'FEC-D arm'
  results_mv_all_signatures[['  ']][results_mv_all_signatures[['  ']] == "APM_score"] = 'tddEC-D arm'
  results_mv_all_signatures[['  ']][results_mv_all_signatures[['  ']] == 'scale(tumdiam)'] = 'Tumor size'
  results_mv_all_signatures[['  ']][results_mv_all_signatures[['  ']] == 'scale(posnodes)'] = 'Lymph nodes'
  
  results_mv_all_signatures_hide = results_mv_all_signatures[results_mv_all_signatures[['  ']] %in% c("Treatment arm", "FEC-D arm", "tddEC-D arm"), , drop = FALSE]
  
  # Pooled adjusted model without the interaction
  pint_formula = as.formula(paste0('Surv(time', endpoint, '/365.25, ', endpoint, '.status) ~ ', 'TREAT + ', sig_cat, ' + ', adjust_terms))
  
  pint_keep = complete.cases(all_clinical[, all.vars(pint_formula), drop = FALSE]) & all_clinical$TREAT %in% c('FEC-D', 'tddEC-D')
  
  pint_data = droplevels(all_clinical[pint_keep, , drop = FALSE])
  
  pint_base = survival::coxph(pint_formula, data = pint_data, na.action = stats::na.fail)
  
  # Add the continuous score × treatment interaction
  pint_full = survival::coxph(update(pint_formula, as.formula(paste0('. ~ . + TREAT:', sig_cat))), data = pint_data, na.action = stats::na.fail)
  
  pint_panther = lmtest::lrtest(pint_base, pint_full)$`Pr(>Chisq)`[2L]
  
  results_mv_all_signatures_hide[['p interaction']] = ''
  results_mv_all_signatures_hide[['p interaction']][results_mv_all_signatures_hide[['  ']] == 'FEC-D arm'] = format_p_regular(pint_panther, digits = 3)
  
  tm3 = forest_theme(core = list(fg_params = list(fontface = fonts_3),
                                 bg_params = list(fill = row_cols),
                                 padding = unit(c(12, 10), 'mm')),
                     # Reference line width/type/color
                     refline_gp = gpar(lwd = 1, lty = 'dashed', col = 'black'),
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = 'dashed',
                     vertline_col = 'grey40')
  
  fp3 = forest(results_mv_all_signatures_hide[ , c('  ', 'N', ' ', 'HR (95% CI)', 'p value', 'p interaction')],
               est = results_mv_all_signatures_hide$est,
               lower = results_mv_all_signatures_hide$low,
               upper = results_mv_all_signatures_hide$hi,
               sizes = 0.6,
               ci_column = 3,
               ref_line = 1,
               arrow_lab = c('Favors higher APM score', 'Favors lower APM score'),
               xlim = c(0, 3),
               ticks_at = c(0, 1, 2, 3),
               theme = tm3)
  
  switch(
    endpoint,
    "BCRFS" = { Fig_2B = fp3 },
    "EFS"   = { Supp_Fig_9A = fp3 },
    "DDFS"  = { Supp_Fig_9B = fp3 },
    "OS"    = { Supp_Fig_9C = fp3 })
}

###### Heterogeneity of the APM score x treatment interaction by ER status and clinical subtype (three-way interaction) ######

# Test whether the continuous-score treatment interaction differs by ER status or clinical subtype using adjusted three-way interaction LRTs.
endpoint = 'BCRFS'

APM_THREE_WAY = rbind(
  data.frame(subgroup = 'ER status',
             three_way_interaction_lrt(all_clinical, endpoint, 'ERpos', adjust_terms, adjust_vars, group_main = FALSE)),
  data.frame(subgroup = 'Clinical subtype',
             three_way_interaction_lrt(all_clinical, endpoint, 'Clinical_Subtype', adjust_terms_no_her2_sub, adjust_vars_no_her2_sub, group_main = TRUE)))

APM_THREE_WAY

write_table(APM_THREE_WAY, file.path(output_dir, paste0('APM_score_treatment_three_way_interaction_', endpoint, '.csv')))

###### Continuous APM score - HR (Forest plot) - Per arm by ER status / Clinical subtype ######

# Display separate-arm score HRs within ER and clinical-subtype groups, with pooled interaction and three-way heterogeneity p-values.
endpoint = 'BCRFS'
hr_column = paste(rep(' ', 30), collapse = ' ')

supp6_panels = list(
  A = list(label_col = 'ER status', by = 'ER status',
           groups = list('ER+' = all_clinical$ERpos == 'Yes',
                         'ER-' = all_clinical$ERpos == 'No')),
  B = list(label_col = 'Clinical subtype', by = 'clinical subtype',
           groups = list('ER+/HER2-' = all_clinical$Clinical_Subtype == 'ER+/HER2-',
                         'HER2+'     = all_clinical$Clinical_Subtype == 'HER2+',
                         'TNBC'      = all_clinical$Clinical_Subtype == 'TNBC')))

for (panel in names(supp6_panels)) {
  spec = supp6_panels[[panel]]
  
  # Subgroup-specific adjustment: no ER term within ER subgroups; no ER/HER2 terms within ER+/HER2- and TNBC; no HER2 term within HER2+
  tab = NULL
  
  for (lab in names(spec$groups)) {
    terms = if (panel == 'A') adjust_terms_no_er_sub else if (lab == 'HER2+') adjust_terms_no_her2_sub else adjust_terms_no_clin_sub
    vars  = if (panel == 'A') adjust_vars_no_er_sub  else if (lab == 'HER2+') adjust_vars_no_her2_sub  else adjust_vars_no_clin_sub
    tab = rbind(tab, apm_score_by_arm_rows(all_clinical[which(spec$groups[[lab]]), , drop = FALSE], endpoint, terms, vars,
                                           paste0(lab, ' patients'), spec$label_col))
  }
  tab[[spec$label_col]] = paste0(tab[[spec$label_col]], '   ')   # room for the bold labels (column width is measured in plain font)
  
  n_blocks = nrow(tab) / 3
  
  # Heterogeneity row: three-way interaction p (all patients) in the p interaction column
  het_row = tab[1, , drop = FALSE]   # a subgroup header row: no estimate
  het_row[[spec$label_col]] = paste0('Heterogeneity across ', if (panel == 'A') 'ER status' else 'subtypes')
  het_row[['p interaction']] = format_p_regular(APM_THREE_WAY$p[APM_THREE_WAY$subgroup == spec$label_col], digits = 3)
  tab = rbind(tab, het_row)
  
  tm_6 = forest_theme(core = list(fg_params = list(fontface = c(rep(c('bold', 'plain', 'plain'), n_blocks), 'italic')),
                                  bg_params = list(fill = c(rep(c('white', 'white', 'slategray1'), n_blocks), 'white')),
                                  padding = unit(c(1, 8), 'mm')),
                      refline_gp = gpar(lwd = 1, lty = 'dashed', col = 'black'),
                      vertline_lwd = 1,
                      vertline_lty = 'dashed',
                      vertline_col = 'grey40')
  
  fp_6 = forest(tab[, c(spec$label_col, 'N', ' ', 'HR (95% CI)', 'p value', 'p interaction')],
                est = tab$est,
                lower = tab$low,
                upper = tab$hi,
                sizes = 0.6,
                ci_column = 3,
                ref_line = 1,
                arrow_lab = c('Favors higher APM score', 'Favors lower APM score'),
                xlim = c(0, 3),
                ticks_at = c(0, 1, 2, 3),
                theme = tm_6)
  
  switch(panel,
         "A" = { Supp_Fig_4A = fp_6 },
         "B" = { Supp_Fig_4B = fp_6 })
}

###### Sensitivity analysis: age added to the Cox models (continuous score and APM group) ######

# Compare models with and without age using the same complete cases per endpoint, for continuous-score and median-group analyses.
age_endpoints = c("BCRFS", "EFS", "DDFS", "OS")
age_model_terms = c("Without age" = adjust_terms,
                    "With age"    = paste(adjust_terms, "scale(age)", sep = " + "))

# Calculate the likelihood-ratio p-value for two nested Cox models.
lrt_p_age = function(reduced, full) {
  ll_r = logLik(reduced); ll_f = logLik(full)
  pchisq(2 * (as.numeric(ll_f) - as.numeric(ll_r)), df = attr(ll_f, "df") - attr(ll_r, "df"), lower.tail = FALSE)
}

# Format a selected Cox coefficient as an HR, 95% confidence interval and Wald p-value.
hr_row = function(fit, term) {
  ci = summary(fit)$conf.int[term, ]
  p  = summary(fit)$coefficients[term, "Pr(>|z|)"]
  sprintf("%.2f (%.2f-%.2f), p=%.3f", ci["exp(coef)"], ci["lower .95"], ci["upper .95"], p)
}

age_sens_rows = list()

for (ep in age_endpoints) {
  
  time_var   = paste0("time", ep)
  status_var = paste0(ep, ".status")
  
  # One complete-case population per endpoint (age required), used by every model below
  d = all_clinical %>%
    dplyr::filter(TREAT %in% c("FEC-D", "tddEC-D")) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(c(time_var, status_var, "APM_score", "APM_group", "age", adjust_vars)), ~ !is.na(.x))) %>%
    dplyr::mutate(time_years = .data[[time_var]] / 365.25,
                  TREAT      = factor(TREAT, levels = c("FEC-D", "tddEC-D")),
                  APM_group  = factor(APM_group, levels = c("APM-low", "APM-high")))
  
  surv_term = paste0("survival::Surv(time_years, ", status_var, ")")
  
  for (model_name in names(age_model_terms)) {
    
    adj = age_model_terms[[model_name]]
    fit = function(rhs, dat = d) survival::coxph(as.formula(paste(surv_term, "~", rhs, "+", adj)), data = dat)
    
    # Continuous score: interaction LRT and within-arm HR per SD
    p_int_cont = lrt_p_age(fit("TREAT + APM_score"), fit("TREAT * APM_score"))
    hr_fec = hr_row(fit("APM_score", dplyr::filter(d, TREAT == "FEC-D")),   "APM_score")
    hr_tdd = hr_row(fit("APM_score", dplyr::filter(d, TREAT == "tddEC-D")), "APM_score")
    
    # APM group: interaction LRT and treatment HR within each group
    p_int_grp = lrt_p_age(fit("TREAT + APM_group"), fit("TREAT * APM_group"))
    hr_high = hr_row(fit("TREAT", dplyr::filter(d, APM_group == "APM-high")), "TREATtddEC-D")
    hr_low  = hr_row(fit("TREAT", dplyr::filter(d, APM_group == "APM-low")),  "TREATtddEC-D")
    
    age_sens_rows[[length(age_sens_rows) + 1]] = data.frame(
      endpoint = ep, model = model_name, n = nrow(d), events = sum(d[[status_var]] == 1),
      p_int_continuous = round(p_int_cont, 3), HR_per_SD_FEC_D = hr_fec, HR_per_SD_tddEC_D = hr_tdd,
      p_int_APM_group = round(p_int_grp, 3), HR_tdd_vs_FEC_APM_high = hr_high, HR_tdd_vs_FEC_APM_low = hr_low)
  }
}

age_sens = dplyr::bind_rows(age_sens_rows)
print(age_sens)
write.csv(age_sens, file.path(output_dir, "Sensitivity_age_APM_interactions.csv"), row.names = FALSE)

###### ----- Continuous scores: APM score, and head-to-head with SET2,3 and SET ER/PR ----- ######

# Endpoints
endpoints = c("BCRFS", "EFS", "DDFS", "OS")

# A. APM_score in all patients
apm_all_sets = list(
  
  APM_score_all = list(
    marker = "APM_score",
    cohort = "All patients",
    data = all_clinical,
    
    # No forced matching with other markers
    matched_vars = character(0),
    
    adjustment_terms = adjust_terms,
    
    adjustment_vars = adjust_vars))

# Head-to-head marker comparison in ER-positive/HER2-negative population
#
# All markers are analysed in identical complete-case patients within each cohort, endpoint and model type.
#
# Clinical_Subtype is excluded from the ER-positive/HER2-negative population.
#
# Marker-specific adjustment:
#   APM_score, all patients: strata(ERpos) + newHER2 + tumdiam + posnodes + histgrade
#
#   APM_score, ER-positive/HER2-negative: tumdiam + posnodes + histgrade
#
#   SET_2_3, ER-positive/HER2-negative: histgrade only
#
#   SET_ER_PR, ER-positive/HER2-negative: tumdiam + posnodes + histgrade

comparison_sets = list(
  
  ## ER-positive/HER2-negative population
  
  APM_score_ERpos_HER2neg = list(marker = "APM_score", cohort = "ER-positive/HER2-negative", data = all_clinical_ERpos_HERneg,
                                 
                                 matched_vars = c("APM_score", "SET_2_3", "SET_ER_PR", "tumdiam", "posnodes", "histgrade"),
                                 
                                 adjustment_terms = c("scale(tumdiam)", "scale(posnodes)", "histgrade"),
                                 adjustment_vars = c("tumdiam", "posnodes", "histgrade")),
  
  SET_2_3_ERpos_HER2neg = list(marker = "SET_2_3", cohort = "ER-positive/HER2-negative", data = all_clinical_ERpos_HERneg,
                               
                               matched_vars = c("APM_score", "SET_2_3", "SET_ER_PR", "tumdiam", "posnodes", "histgrade"),
                               
                               # SET2,3 already incorporates tumour size and nodal burden
                               adjustment_terms = c("histgrade"),
                               adjustment_vars = c("histgrade")),
  
  SET_ER_PR_ERpos_HER2neg = list(
    marker = "SET_ER_PR",
    cohort = "ER-positive/HER2-negative",
    data = all_clinical_ERpos_HERneg,
    
    matched_vars = c("APM_score", "SET_2_3", "SET_ER_PR", "tumdiam", "posnodes", "histgrade"),
    
    adjustment_terms = c("scale(tumdiam)", "scale(posnodes)", "histgrade"),
    adjustment_vars = c( "tumdiam", "posnodes", "histgrade")))

comparison_sets = lapply(comparison_sets, function(s) c(s, standardize = TRUE))

# APM score in the complete PANTHER cohort
apm_all_prognostic_results = run_biomarker_models(analysis_sets = apm_all_sets, endpoints = endpoints)
apm_all_interaction_results = run_interaction_models(analysis_sets = apm_all_sets, endpoints = endpoints)

# Head-to-head comparison of APM score, SET2,3 and SET ER/PR
comparison_prognostic_results = run_biomarker_models(analysis_sets = comparison_sets, endpoints = endpoints)
comparison_interaction_results = run_interaction_models(analysis_sets = comparison_sets, endpoints = endpoints)

# Create wide-format output tables
apm_all_prognostic_wide = make_biomarker_wide(results = apm_all_prognostic_results, endpoints = endpoints)
apm_all_interaction_wide = make_interaction_wide(results = apm_all_interaction_results, endpoints = endpoints)

comparison_prognostic_wide = make_biomarker_wide(results = comparison_prognostic_results, endpoints = endpoints)
comparison_interaction_wide = make_interaction_wide(results = comparison_interaction_results, endpoints = endpoints)

results_xlsx = list(
  "APM_score_all_prognostic" = apm_all_prognostic_wide,
  "APM_score_all_interaction" = apm_all_interaction_wide,
  "APM_SET_comparison_prognostic" = comparison_prognostic_wide,
  "APM_SET_comparison_interaction" = comparison_interaction_wide)

openxlsx::write.xlsx(x = results_xlsx, file = file.path(output_dir, "APM_Cox_model_prognostic_and_interaction_results.xlsx"), overwrite = TRUE)

###### ----- Association Analysis for the binary APM group -----
###### SubGroup Analysis (median cut-offs) - HR (Forest plot) - APM score ######

# Estimate unadjusted and adjusted treatment HRs overall and within cohort-median APM groups; test treatment-by-group interactions for each endpoint.
sig_cat = 'APM_group'
hr_column = paste(rep(' ', 30), collapse = ' ')
for (endpoint in c('BCRFS', 'EFS', 'DDFS', 'OS')) {

  fonts = c('bold')
  row_cols = c('white')

  results_mv_all_signatures = data.frame('Biomarker' = character(), `FEC-D` = character(), `tddEC-D` = character(),
                                         'est' =  numeric(), 'low' = numeric(), 'hi' = numeric(), 'se' = numeric(), ` ` = character(), `HR (95% CI)` = numeric(),
                                         'p value' = numeric(), 'p interaction' = numeric(), check.names = FALSE)

  # All patients (MV)
  model_all_mv = coxph(as.formula(paste0("Surv(time", endpoint, "/365.25, ", endpoint, ".status) ~ TREAT + ", adjust_terms)), data = all_clinical, model = TRUE)
  model_all_mv_results = tidy(model_all_mv, exponentiate = TRUE, conf.int = TRUE)
  n_mv = table(factor(as.character(model_all_mv$model$TREAT), levels = c("FEC-D", "tddEC-D")))

  fonts = c('plain')

  results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('Biomarker' = 'All patients',
                                                                          `FEC-D` = as.integer(n_mv[["FEC-D"]]),
                                                                          `tddEC-D` = as.integer(n_mv[["tddEC-D"]]),
                                                                          'est' = round(as.numeric(model_all_mv_results[1, 'estimate']),2),
                                                                          'low' = round(as.numeric(model_all_mv_results[1, 'conf.low']),2),
                                                                          'hi' = round(as.numeric(model_all_mv_results[1, 'conf.high']),2),
                                                                          'se' = round(as.numeric(model_all_mv_results[1, 'std.error']),2),
                                                                          ` ` = hr_column,
                                                                          `HR (95% CI)` = format_hr_ci(as.numeric(model_all_mv_results[1, 'estimate']), as.numeric(model_all_mv_results[1, 'conf.low']), as.numeric(model_all_mv_results[1, 'conf.high'])),
                                                                          'p value' = format_p_regular(as.numeric(model_all_mv_results[1, 'p.value']), digits = 3),
                                                                          'p interaction' = '', check.names = FALSE))

  fonts = c(fonts, 'bold')
  new_col = if (tail(row_cols, 1) == "white") { rep("slategray1", length(unique(all_clinical[, sig_cat])) + 1) } else { rep("white", length(unique(all_clinical[, sig_cat])) + 1) }
  row_cols = c(row_cols, new_col)

  results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('Biomarker' = sig_cat, `FEC-D` = '', `tddEC-D` = '',
                                                                          'est' = NA, 'low' = NA, 'hi' = NA, 'se' = NA, ` ` = hr_column, `HR (95% CI)` = '',
                                                                          'p value' = '', 'p interaction' = '', check.names = FALSE))

  # interaction test for Adjusted Cox model
  res.cox_base = coxph(as.formula(paste0('Surv(time', endpoint, '/365.25, ', endpoint, '.status) ~ TREAT + ', sig_cat, ' + ', adjust_terms)), data = all_clinical)
  res.cox_int = coxph(as.formula(paste0('Surv(time', endpoint, '/365.25, ', endpoint, '.status) ~ TREAT * ', sig_cat, ' + ', adjust_terms)), data = all_clinical)
  mv_interaction_lr = lrtest(res.cox_base, res.cox_int)
  mv_interaction_p_value = mv_interaction_lr$Pr[2]

  for (bt in sort(unique(all_clinical[, sig_cat]))) {

    fonts = c(fonts, 'plain')

    # Multivariable model
    data_cox = all_clinical[all_clinical[, sig_cat] == bt & !is.na(all_clinical$TREAT) & !is.na(all_clinical$ERpos) & !is.na(all_clinical$newHER2) & !is.na(all_clinical$tumdiam) & !is.na(all_clinical$posnodes) & !is.na(all_clinical$histgrade),]
    res.cox = coxph(as.formula(paste0('Surv(time', endpoint, '/365.25, ', endpoint, '.status) ~ TREAT + ', adjust_terms)), data = data_cox, model = TRUE)

    cox_model_results = tidy(res.cox, exponentiate = TRUE, conf.int = TRUE)

    results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('Biomarker' = paste0(bt),
                                                                            `FEC-D` = sum(res.cox$model$TREAT == 'FEC-D'),
                                                                            `tddEC-D` = sum(res.cox$model$TREAT == 'tddEC-D'),
                                                                            'est' = NA,
                                                                            'low' = NA,
                                                                            'hi' = NA,
                                                                            'se' = NA,
                                                                            ` ` = hr_column,
                                                                            `HR (95% CI)` = '',
                                                                            'p value' = '',
                                                                            'p interaction' = '', check.names = FALSE))

    for (term in cox_model_results$term) {

      results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('Biomarker' = term,
                                                                              `FEC-D` = '',
                                                                              `tddEC-D` = '',
                                                                              'est' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'estimate']),2),
                                                                              'low' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'conf.low']),2),
                                                                              'hi' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'conf.high']),2),
                                                                              'se' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'std.error']),2),
                                                                              ` ` = hr_column,
                                                                              `HR (95% CI)` = format_hr_ci(as.numeric(cox_model_results[cox_model_results$term == term, 'estimate']), as.numeric(cox_model_results[cox_model_results$term == term, 'conf.low']), as.numeric(cox_model_results[cox_model_results$term == term, 'conf.high'])),
                                                                              'p value' = format_p_regular(as.numeric(cox_model_results[cox_model_results$term == term, 'p.value']), digits = 3),
                                                                              'p interaction' = if (which(bt == sort(unique(all_clinical[,sig_cat]))) == 1 & which(term == cox_model_results$term) == 1) {format_p_regular(mv_interaction_p_value, digits = 3)} else {''}, check.names = FALSE))
    }
  }
  results_mv_all_signatures$Biomarker[results_mv_all_signatures$Biomarker == 'APM_group'] = 'APM score'
  results_mv_all_signatures$Biomarker[results_mv_all_signatures$Biomarker == 'scale(tumdiam)'] = 'Tumor size'
  results_mv_all_signatures$Biomarker[results_mv_all_signatures$Biomarker == 'scale(posnodes)'] = 'Lymph nodes'

  results_mv_all_signatures_hide = results_mv_all_signatures[results_mv_all_signatures$Biomarker %in% c("All patients", "APM score", "APM-low", "APM-high", "TREATtddEC-D"), , drop = FALSE]

  df = results_mv_all_signatures_hide
  # rows containing the treatment results
  treat_rows = which(df$Biomarker == "TREATtddEC-D")
  # rows where the values should be copied (one row above)
  target_rows = treat_rows - 1
  # columns to copy
  cols = c("est","low","hi","se","HR (95% CI)","p value","p interaction")
  df[target_rows, cols] = df[treat_rows, cols]
  # remove the treatment rows
  df = df[-treat_rows, ]
  results_mv_all_signatures_hide = df

  results_mv_all_signatures_hide = results_mv_all_signatures_hide %>% dplyr::rename('  ' = Biomarker)

  tm3 = forest_theme(core = list(fg_params = list(fontface = fonts),
                                 bg_params = list(fill = row_cols),
                                 padding = unit(c(1, 10), 'mm')),
                     # Reference line width/type/color
                     refline_gp = gpar(lwd = 1, lty = 'dashed', col = 'black'),
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = 'dashed',
                     vertline_col = 'grey40')

  fp3 = forest(results_mv_all_signatures_hide[, c('  ', 'FEC-D', 'tddEC-D', ' ', 'HR (95% CI)', 'p value', 'p interaction')],
               est = results_mv_all_signatures_hide$est,
               lower = results_mv_all_signatures_hide$low,
               upper = results_mv_all_signatures_hide$hi,
               sizes = 0.6,
               ci_column = 4,
               ref_line = 1,
               vert_line = round(as.numeric(model_all_mv_results[1, 'estimate']),2),
               arrow_lab = c('Favors tddEC-D', 'Favors FEC-D'),
               xlim = c(0, 4),
               ticks_at = c(0, 1, 2, 4),
               theme = tm3)

  if (sig_cat == "APM_group") {
    switch(
      endpoint,
      "BCRFS" = { Fig_2A = fp3 },
      "EFS"   = { Supp_Fig_7A = fp3 },
      "DDFS"  = { Supp_Fig_7B = fp3 },
      "OS"    = { Supp_Fig_7C = fp3 })
  }
}

###### Survival analysis - KM curves (by APM score) ######

# Generate KM curves by APM group, by treatment within each APM group
sig_cat = 'APM_group'
for (endpoint in c('BCRFS', 'EFS', 'DDFS', 'OS')) {

  # All (by APM score)
  km_all_apm = plot_km(dataset = all_clinical, predictor = sig_cat, predictor_name = 'APM score', cohort = 'All', surv_endpoint = endpoint, colors = group_col_apm)

  if (endpoint == "BCRFS") {
    Supp_Fig_5A = km_all_apm
  } else if (endpoint == "EFS") {
    Supp_Fig_5B = km_all_apm
  } else if (endpoint == "DDFS") {
    Supp_Fig_5C = km_all_apm
  } else if (endpoint == "OS") {
    Supp_Fig_5D = km_all_apm
  }

  # APM-high Subgroup (by arm)
  km_apm_high = plot_km(dataset = all_clinical[all_clinical[[sig_cat]] == 'APM-high', ], predictor = 'TREAT', predictor_name = 'Treatment Arm', cohort = 'APM-high', surv_endpoint = endpoint, colors = group_col_treat)

  # APM-low Subgroup (by arm)
  km_apm_low = plot_km(dataset = all_clinical[all_clinical[[sig_cat]] == 'APM-low', ], predictor = 'TREAT', predictor_name = 'Treatment Arm', cohort = 'APM-low', surv_endpoint = endpoint, colors = group_col_treat)

  switch(
    endpoint,
    "BCRFS" = { Fig_3A = km_apm_high ; Fig_3B = km_apm_low },
    "EFS"   = { Supp_Fig_8A = km_apm_high ; Supp_Fig_8B = km_apm_low },
    "DDFS"  = { Supp_Fig_8C = km_apm_high ; Supp_Fig_8D = km_apm_low },
    "OS"    = { Supp_Fig_8E = km_apm_high ; Supp_Fig_8F = km_apm_low })
}

###### ----- APM score vs. Clinical Subtypes / PAM50 -----
###### APM score (continuous) vs Clinical_Subtype ######

# Compare continuous APM scores across clinical subtypes using a Kruskal-Wallis test.

plot_df1 = all_clinical %>% filter(!is.na(APM_score), !is.na(Clinical_Subtype))

# Plot
p1 = ggplot(
  plot_df1,
  aes(x = Clinical_Subtype, y = APM_score, fill = Clinical_Subtype)) +

  geom_boxplot(outlier.shape = NA, alpha = 0.75) +

  geom_jitter(width = 0.15, alpha = 0.45, size = 1.5) +

  stat_compare_means(method = "kruskal.test") +

  labs(
    title = "APM score across clinical subtypes",
    x = "Clinical subtype",
    y = "APM score") +

  theme_bw(base_size = 14) +

  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))

# Save plot
ggsave(filename = file.path(output_dir, "APM_vs_ClinicalSubtype.png"), plot = p1, width = 8, height = 6,  dpi = 300, bg = "white")

###### APM score (continuous) vs PAM50 subtype ######

# Compare continuous APM scores across intrinsic subtypes using a Kruskal-Wallis test.

plot_df2 = all_clinical %>% filter(!is.na(APM_score), !is.na(subtypecd_ssp.subtype))

# Plot
p2 = ggplot(
  plot_df2,
  aes(x = subtypecd_ssp.subtype, y = APM_score, fill = subtypecd_ssp.subtype)) +

  geom_boxplot(outlier.shape = NA, alpha = 0.75) +

  geom_jitter(width = 0.15, alpha = 0.45, size = 1.8) +

  stat_compare_means(method = "kruskal.test") +

  labs(
    title = "APM score across intrinsic subtypes",
    x = "Intrinsic subtypes",
    y = "APM score") +

  theme_bw(base_size = 14) +

  theme(legend.position = "none")

# Save plot
ggsave(filename = file.path(output_dir, "APM_vs_PAM50Subtype.png"), plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

###### APM group (median, Low vs. High) vs PAM50 subtype ######

# Display APM-group proportions within intrinsic subtypes and test their association using Pearson's chi-square test.

plot_df4 = all_clinical %>% filter(!is.na(APM_group), !is.na(subtypecd_ssp.subtype))

n_panther_plot = nrow(plot_df4)

# Contingency table
table_pam50_cat = table(
  plot_df4$subtypecd_ssp.subtype,
  plot_df4$APM_group)

# Chi-square test
chisq_pam50_cat = chisq.test(table_pam50_cat)

# Plot dataframe
plot_counts4 = plot_df4 %>%
  dplyr::count(subtypecd_ssp.subtype, APM_group) %>%
  dplyr::group_by(subtypecd_ssp.subtype) %>%
  dplyr::mutate(Percentage = 100 * n / sum(n))

# Plot
p4 = ggplot(
  plot_counts4,
  aes(
    x = subtypecd_ssp.subtype,
    y = Percentage,
    fill = APM_group)) +

  geom_col(width = 0.75) +

  geom_text(
    aes(label = paste0(formatC(Percentage, format = "f", digits = 1), "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4) +

  ggplot2::scale_fill_manual(
    values = group_col_apm,
    breaks = c("APM-low", "APM-high"),
    labels = c("APM-low", "APM-high"),
    drop = FALSE) +

  labs(
    title = "APM groups across intrinsic subtypes in the PANTHER cohort",
    subtitle = paste0("N = ", n_panther_plot, "; Chi-square ", p_lab(chisq_pam50_cat$p.value)),
    x = "Intrinsic subtype",
    y = "Percentage (%)",
    fill = "APM group") +

  theme_bw(base_size = 14)

# Save plot
ggsave(filename = file.path(output_dir, "APM_group_vs_PAM50Subtype.png"), plot = p4, width = 8, height = 6, dpi = 300, bg = "white")

Supp_Figure_10A = p4

###### APM score (continuous, per SD) x treatment, adjusted for intrinsic subtype (PAM50) ######

# Test the continuous APM-score treatment interaction after adjustment for intrinsic subtype alone and together with clinical covariates.

pam50_adj_apm = do.call(rbind, lapply(c("BCRFS", "EFS", "DDFS", "OS"), function(ep) {
  s = paste0("Surv(time", ep, "/365.25, ", ep, ".status) ~ ")
  lrt = function(adj) {
    b = coxph(as.formula(paste0(s, "TREAT + APM_score + ", adj)), data = all_clinical)
    f = coxph(as.formula(paste0(s, "TREAT * APM_score + ", adj)), data = all_clinical)
    stopifnot(b$n == f$n)   # same patients in both models
    c(n = f$n, p = lmtest::lrtest(b, f)$`Pr(>Chisq)`[2])
  }
  a = lrt("subtypecd_ssp.subtype")
  m = lrt(paste("subtypecd_ssp.subtype +", adjust_terms))
  data.frame(endpoint = ep,
             n_subtype_only = a[["n"]], p_interaction_subtype_only = format_p_regular(a[["p"]]),
             n_subtype_clinical = m[["n"]], p_interaction_subtype_clinical = format_p_regular(m[["p"]]))
}))

pam50_adj_apm

write_table(pam50_adj_apm, file.path(output_dir, "APM_score_interaction_adjusted_for_PAM50.csv"))

###### Intrinsic subtype x treatment interaction (LRT), all endpoints ######

# Unadjusted, and adjusted for tumour size, nodes and grade (ER/HER2 left out because intrinsic subtype largely captures them)

pam50_interaction = do.call(rbind, lapply(c("BCRFS", "EFS", "DDFS", "OS"), function(endpoint) {
  surv = paste0("Surv(time", endpoint, "/365.25, ", endpoint, ".status) ~ ")
  lrt = function(adj = "") {
    base = coxph(as.formula(paste0(surv, "TREAT + subtypecd_ssp.subtype", adj)), data = all_clinical)
    full = coxph(as.formula(paste0(surv, "TREAT * subtypecd_ssp.subtype", adj)), data = all_clinical)
    stopifnot(base$n == full$n)   # same patients in both models
    c(n = full$n, events = full$nevent, p = lmtest::lrtest(base, full)$`Pr(>Chisq)`[2])
  }
  u = lrt()
  a = lrt(paste(" +", adjust_terms_no_clin_sub))
  data.frame(endpoint,
             n_unadj = u[["n"]], p_interaction_unadjusted = format_p_regular(u[["p"]]),
             n_adj = a[["n"]], events_adj = a[["events"]], p_interaction_adjusted = format_p_regular(a[["p"]]))
}))

pam50_interaction

write_table(pam50_interaction, file.path(output_dir, "PAM50_treatment_interaction.csv"))

###### ----- APM group vs. 21-gene RS group ----- ######

# Calculate RNA-seq-derived 21-gene recurrence scores in ER+/HER2- tumours and test the association between recurrence-risk groups and APM groups.

## 21-gene RS (for ER+/HER2- patients only)

data(sig.oncotypedx)
sig.oncotypedx['CTSL2', 'symbol'] = 'CTSV'
rownames(sig.oncotypedx)[rownames(sig.oncotypedx) == 'CTSL2'] = 'CTSV'

oncotype_rs = OncotypeDX(data = t(tmm_normalized_tpm_symbols_ERpos_HERneg[sig.oncotypedx$symbol,]), annot = sig.oncotypedx, do.mapping = TRUE)
OncotypeDx = data.frame('sample_id' = names(oncotype_rs$score), 'OncotypeDx_score' = oncotype_rs$score, 'OncotypeDx_risk_group' = oncotype_rs$risk, check.names = FALSE)

all_clinical_ERpos_HERneg = left_join(all_clinical_ERpos_HERneg, OncotypeDx, by = 'sample_id')

## 21-gene RS risk groups vs APM groups (ER+/HER2- patients)
rs_apm = all_clinical_ERpos_HERneg[!is.na(all_clinical_ERpos_HERneg$OncotypeDx_risk_group) & !is.na(all_clinical_ERpos_HERneg$APM_group), ]
rs_apm$RS_group = factor(rs_apm$OncotypeDx_risk_group, levels = c(0, 0.5, 1), labels = c("Low (RS < 18)", "Intermediate (RS 18 to < 31)", "High (RS >= 31)"))

rs_tab = table(rs_apm$RS_group, rs_apm$APM_group)
rs_pct = round(100 * prop.table(rs_tab, margin = 2), 1)                                            # % within each APM group
rs_p = chisq.test(rs_tab)$p.value                                                                  # 3 RS groups x 2 APM groups
rs_p_collapsed = chisq.test(table(rs_apm$RS_group != "Low (RS < 18)", rs_apm$APM_group))$p.value   # low vs intermediate/high

summary_rs = data.frame(
  section = "RNA-seq-derived 21-gene RS risk group vs APM group (ER+/HER2-)",
  measure = c(rownames(rs_tab), "Chi-square, 3 RS groups x APM group", "Chi-square, RS low vs intermediate/high x APM group"),
  n = c(as.integer(rowSums(rs_tab)), sum(rs_tab), sum(rs_tab)),
  events = NA_integer_,
  result = c(sprintf("APM-low %d (%.1f%%); APM-high %d (%.1f%%)", rs_tab[, "APM-low"], rs_pct[, "APM-low"], rs_tab[, "APM-high"], rs_pct[, "APM-high"]), "", ""),
  p_value = c(rep("", nrow(rs_tab)), format_p_regular(rs_p), format_p_regular(rs_p_collapsed)),
  row.names = NULL)

summary_rs

write_table(summary_rs, file.path(output_dir, "RS_group_vs_APM_group.csv"))

###### ----- Deconvolution methods -----
###### MCPcounter  ######

# Estimate immune and stromal population scores from PANTHER expression and merge them with clinical records.
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1070-5

MCPcounter_gene = read.table(curl::curl('https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt'), sep='\t', stringsAsFactors = FALSE, header = TRUE, colClasses = 'character', check.names = FALSE)
MCPcounter_gene$`HUGO symbols`[MCPcounter_gene$`HUGO symbols` == 'C'] = 'SDPR'
MCPcounter_gene$`HUGO symbols`[MCPcounter_gene$`HUGO symbols` == 'DIPK2B'] = 'CXorf36'

MCPcounter = as.data.frame(t(MCPcounter.estimate(tmm_normalized_tpm_symbols, featuresType = 'HUGO_symbols', genes = MCPcounter_gene)))
colnames(MCPcounter) = paste0('MCPcounter_', colnames(MCPcounter))
colnames(MCPcounter) = gsub(' ', '_', colnames(MCPcounter))
MCPcounter$sample_id = rownames(MCPcounter)

all_clinical = left_join(all_clinical, MCPcounter, by = 'sample_id')

###### Danaher ######

# Calculate mean expression of available cell-type marker genes and merge the resulting population scores with clinical records.
# https://jitc.bmj.com/content/5/1/18

# load gene annotation
gannot = read.csv(file.path(data_dir, "annotation", "gannot-cip-junecelltype.csv"), row.names = 1)

# Genes and cell types included in the TCGA2/Danaher annotation
celltypegenes = rownames(gannot)[!is.na(gannot$Cell.Type.TCGA2)]
celltypes = unique(gannot$Cell.Type.TCGA2)
celltypes = celltypes[!is.na(celltypes)]

# Construct gene-by-cell-type membership matrix
celltypematrix = matrix(0, nrow = length(celltypegenes), ncol = length(celltypes), dimnames = list(celltypegenes, celltypes))

for (i in seq_along(celltypes)) {
  tempgenes = rownames(gannot)[!is.na(gannot$Cell.Type.TCGA2) & gannot$Cell.Type.TCGA2 == celltypes[i]]
  celltypematrix[tempgenes, celltypes[i]] = 1
}

# Confirm marker availability in RNA-seq expression matrix
available_celltypegenes = intersect(rownames(tmm_normalized_tpm_symbols), rownames(celltypematrix))

missing_celltypegenes = setdiff(rownames(celltypematrix), rownames(tmm_normalized_tpm_symbols))

if (length(missing_celltypegenes) > 0) { warning("Some Danaher marker genes are absent from the expression matrix: ", paste(missing_celltypegenes, collapse = ", ")) }

# Restrict annotation matrix to genes available in expression data
celltypematrix_available = celltypematrix[available_celltypegenes, , drop = FALSE]

# Remove a cell type only if none of its marker genes is available
celltypematrix_available = celltypematrix_available[ , colSums(celltypematrix_available) > 0, drop = FALSE]

# Number of genes contributing to each immune population score
danaher_gene_counts = colSums(celltypematrix_available)

# Calculate Danaher immune population scores
# Expression matrix: samples × genes
df = t(tmm_normalized_tpm_symbols[available_celltypegenes, , drop = FALSE])

# Weight matrix: genes × cell types
# Each marker receives equal weight within its immune population
celltype_weights = sweep(celltypematrix_available, 2, danaher_gene_counts, FUN = "/")

# Danaher score = mean expression of marker genes for each cell type
Danaher = as.data.frame(df %*% celltype_weights)

# Standardise output column names
colnames(Danaher) = paste0("Danaher_", colnames(Danaher))
colnames(Danaher) = gsub("-", "_", colnames(Danaher), fixed = TRUE)
colnames(Danaher) = gsub(" ", "_", colnames(Danaher), fixed = TRUE)
# Add sample identifier
Danaher$sample_id = rownames(Danaher)

all_clinical = left_join(all_clinical, Danaher, by = 'sample_id')

###### Boxplots for deconvolution methods ######

# Compare population scores between APM groups using Wilcoxon tests, with BH correction across populations within each method.

mcpcounter_order = c("B lineage", "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "NK cells", "Neutrophils", "Monocytic lineage", "Myeloid dendritic cells", "Endothelial cells", "Fibroblasts")

danaher_order = c("CD45", "B cells", "T cells", "CD8 T cells", "Exhausted CD8", "Cytotoxic cells", "Th1 cells", "Treg", "NK cells", "NK CD56dim cells", "Mast cells", "DC", "Macrophages", "Neutrophils" )

mcpcounter_results = plot_apm_immune_scores(
  data = all_clinical,
  score_prefix = "MCPcounter_",
  group_var = "APM_group",
  y_label = "MCP-counter score",
  title = NULL,
  feature_order = mcpcounter_order,
  display_p = "p.adj",
  p_adjust_method = "BH",
  y_multiplier = 1.30)

mcpcounter_results$plot

Supp_Figure_11A = mcpcounter_results$plot
ggsave(filename = file.path(output_dir, "ImmDeconv_MCPcounter_APM_groups_boxplots.png"),
       plot = mcpcounter_results$plot, width = 12, height = 7, dpi = 300, bg = "white")

danaher_results = plot_apm_immune_scores(
  data = all_clinical,
  score_prefix = "Danaher_",
  group_var = "APM_group",
  y_label = "Danaher immune cell score",
  title = NULL,
  feature_order = danaher_order,
  display_p = "p.adj",
  p_adjust_method = "BH",
  y_multiplier = 1.30)

danaher_results$plot

Supp_Figure_12A = danaher_results$plot
ggsave(filename = file.path(output_dir, "ImmDeconv_Danaher_APM_groups_boxplots.png"),
       plot = danaher_results$plot, width = 12, height = 7, dpi = 300, bg = "white")

###### ----- 14-gene immune IgG signature -----
# https://bmcmedgenomics.biomedcentral.com/articles/10.1186/1755-8794-4-3

# Score on all RNA-matched patients: mean log2 TMM-TPM of the available signature genes.
# The score is standardized over these patients (HR per SD) and split at the cohort median (High >= median), fixed before any endpoint or covariate exclusion.

immune_igg_genes = c('CD27', 'CD79A', 'CXCL8', 'HLA-C', 'IGKC', 'IGLC6', 'IGLV3-25', 'IL2RG', 'JCHAIN', 'LAX1', 'NTN3', 'PIM2', 'POU2AF1', 'TNFRSF17')

igg_ids = intersect(as.character(all_clinical$sample_id), colnames(tmm_normalized_tpm_symbols))
igg_expr = as.matrix(tmm_normalized_tpm_symbols[, igg_ids, drop = FALSE])
igg_expr = igg_expr[rowSums(!is.finite(igg_expr)) == 0, , drop = FALSE]
igg_expr = igg_expr[apply(igg_expr, 1L, stats::sd) > 1e-8, , drop = FALSE]

igg_genes_used = intersect(immune_igg_genes, rownames(igg_expr))
if (length(igg_genes_used) < length(immune_igg_genes))
  warning("IgG genes missing or constant: ", paste(setdiff(immune_igg_genes, igg_genes_used), collapse = ", "))

igg_scores = data.frame(
  sample_id = igg_ids,
  immune_igg = colMeans(igg_expr[igg_genes_used, , drop = FALSE]))
stopifnot(!anyNA(igg_scores))

# Standardize and add to all_clinical by sample ID (rerunning replaces the column)
all_clinical$immune_igg = as.numeric(scale(igg_scores$immune_igg))[match(as.character(all_clinical$sample_id), igg_scores$sample_id)]

# Median split (cohort median): treatment HR within Low / High and interaction LRT
igg_median_results = dplyr::bind_rows(lapply("immune_igg", function(b)
  get_median_cutoff_pint(data = all_clinical, biomarker = b, endpoint = "BCRFS", adjust_vars = adjust_vars, adjust_terms = adjust_terms)))

# Continuous score (per SD): interaction LRT, as for the APM score
igg_continuous_results = data.frame(
  biomarker = c("immune_igg"),
  p_interaction_continuous_adj = vapply("immune_igg", function(b)
    get_linear_interaction_pvalues(data = all_clinical, time_col = "timeBCRFS", status_col = "BCRFS.status", continuous_var = b,
                                   adjust_vars = adjust_vars, adjust_terms = adjust_terms)$p_adjusted, numeric(1)), row.names = NULL)

igg_median_results[, c("biomarker", "group", "n", "events", "HR_CI", "p_value_print", "p_interaction_print")]
igg_continuous_results

## Summary table (BCRFS, adjusted LRT): continuous score (per SD) and cohort-median split
labels_igg = c(immune_igg = "IgG signature (mean)")
igg_median_df = as.data.frame(igg_median_results)

summary_igg = do.call(rbind, lapply(names(labels_igg), function(b) {
  m = igg_median_df[igg_median_df$biomarker == b, ]
  lo = m[m$group == "Low", ]
  hi = m[m$group == "High", ]
  data.frame(measure = paste0(labels_igg[[b]], c(", continuous score (per SD)", ", median split")),
             n = m$N_total[1], events = sum(m$events),
             result = c("", sprintf("tddEC-D vs FEC-D: low %s, p %s; high %s, p %s", lo$HR_CI, lo$p_value_print, hi$HR_CI, hi$p_value_print)),
             p_interaction = c(format_p_regular(igg_continuous_results$p_interaction_continuous_adj[igg_continuous_results$biomarker == b]),
                               m$p_interaction_print[1]),
             row.names = NULL)
}))
summary_igg

write_table(summary_igg, file.path(output_dir, "IgG_signature_treatment_interaction.csv"))

###### ----- Immune signals and differential treatment benefit (exploratory) -----

# Q1  Which immune scores show a treatment interaction, as the APM score does?
# Q2  Does the APM interaction remain when an immune score is given its own treatment interaction?
# Q3  Cytotoxic vs humoral: which interaction remains when both are in one model?
# BCRFS; adjusted Cox models as for the APM score (adjust_terms); all scores per SD; same patients in every model.

imm_endpoint = "BCRFS"
imm_time = paste0("time", imm_endpoint)
imm_status = paste0(imm_endpoint, ".status")
surv_imm = paste0("Surv(", imm_time, "/365.25, ", imm_status, ") ~ ")

immune_vars = c(grep("^(MCPcounter|Danaher)_", names(all_clinical), value = TRUE), "immune_igg")
immune_vars = immune_vars[vapply(immune_vars, function(v) isTRUE(sd(all_clinical[[v]], na.rm = TRUE) > 0), logical(1))]

# Q3 pairs (Define cytotoxic/CD8 and humoral score pairs for joint interaction models.)
cyto = c("MCPcounter_Cytotoxic_lymphocytes", "MCPcounter_CD8_T_cells", "Danaher_Cytotoxic_cells", "Danaher_CD8_T_cells")
humo = c("immune_igg", "MCPcounter_B_lineage", "Danaher_B_cells")
missing_q3 = setdiff(c(cyto, humo), immune_vars)
if (length(missing_q3)) stop("Not found: ", paste(missing_q3, collapse = ", "), "\nAvailable: ", paste(immune_vars, collapse = ", "))

# Standardize over all RNA-matched patients (as for the APM score), then one complete-case population
d_imm = all_clinical
z_vars = paste0(c("APM_score", immune_vars), "_z")
for (v in c("APM_score", immune_vars)) d_imm[[paste0(v, "_z")]] = as.numeric(scale(d_imm[[v]]))
need = c(imm_time, imm_status, "TREAT", adjust_vars, z_vars)
d_imm = droplevels(d_imm[complete.cases(d_imm[, need]), , drop = FALSE])
cat("Immune interaction analyses:", nrow(d_imm), "patients,", sum(d_imm[[imm_status]]), "events\n")

# Fit an immune-score Cox model with the common clinical adjustment specification.
fit_imm = function(rhs, data = d_imm) coxph(as.formula(paste0(surv_imm, rhs, " + ", adjust_terms)), data = data)

# Calculate the likelihood-ratio p-value for two nested immune-score models.
lrt_imm = function(small, big) {
  chisq = 2 * (as.numeric(logLik(big)) - as.numeric(logLik(small)))
  pchisq(chisq, df = attr(logLik(big), "df") - attr(logLik(small), "df"), lower.tail = FALSE)
}

# Format a selected immune-score coefficient as an HR and 95% confidence interval.
hr_imm = function(fit, term) {
  r = summary(fit)$conf.int[term, ]
  trimws(format_hr_ci(r[["exp(coef)"]], r[["lower .95"]], r[["upper .95"]]))
}

# Convert immune-score variable names into readable table labels.
imm_label = function(z) {
  if (z == "immune_igg_z") return("14-gene IgG")
  gsub("_", " ", sub("^(MCPcounter|Danaher)_", "\\1: ", sub("_z$", "", z)))
}

# Format columns whose names begin with p_ for display and export.
fmt_p = function(df) { for (v in grep("^p_", names(df), value = TRUE)) df[[v]] = format_p_regular(df[[v]]); df }

## Q1: treatment x score interaction, one score at a time (APM score first, as reference)
q1 = do.call(rbind, lapply(z_vars, function(z) {
  main = fit_imm(paste("TREAT +", z))
  int  = fit_imm(paste("TREAT *", z))
  b = coef(int)[[paste0("TREATtddEC-D:", z)]]
  data.frame(score = imm_label(z),
             rho_with_APM = round(cor(d_imm[[z]], d_imm$APM_score_z, method = "spearman"), 2),
             HR_per_SD_FEC_D   = hr_imm(fit_imm(z, d_imm[d_imm$TREAT == "FEC-D", ]), z),
             HR_per_SD_tddEC_D = hr_imm(fit_imm(z, d_imm[d_imm$TREAT == "tddEC-D", ]), z),
             direction = ifelse(b < 0, "higher favours tddEC-D", "higher favours FEC-D"),
             p_interaction = lrt_imm(main, int))
}))
is_imm = q1$score != "APM score"
q1$p_interaction_BH = NA
q1$p_interaction_BH[is_imm] = p.adjust(q1$p_interaction[is_imm], method = "BH")
q1 = q1[order(is_imm, q1$p_interaction), ]

## Q2: APM interaction with and without each immune score's own treatment interaction
q2 = do.call(rbind, lapply(setdiff(z_vars, "APM_score_z"), function(z) {
  full     = fit_imm(paste("TREAT * APM_score_z + TREAT *", z))
  no_apm_x = fit_imm(paste("TREAT *", z, "+ APM_score_z"))
  no_imm_x = fit_imm(paste("TREAT * APM_score_z +", z))
  data.frame(score = imm_label(z),
             p_APM_interaction_given_immune = lrt_imm(no_apm_x, full),
             p_immune_interaction_given_APM = lrt_imm(no_imm_x, full),
             APM_interaction_retained_pct = round(100 * coef(full)[["TREATtddEC-D:APM_score_z"]] / coef(no_imm_x)[["TREATtddEC-D:APM_score_z"]]))
}))
q2 = q2[order(q2$p_APM_interaction_given_immune, decreasing = TRUE), ]

## Q3: cytotoxic vs humoral, both interactions in one model
pairs = expand.grid(cytotoxic = cyto, humoral = humo, stringsAsFactors = FALSE)
q3 = do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  c_z = paste0(pairs$cytotoxic[i], "_z"); h_z = paste0(pairs$humoral[i], "_z")
  full   = fit_imm(paste("TREAT *", c_z, "+ TREAT *", h_z))
  no_c_x = fit_imm(paste("TREAT *", h_z, "+", c_z))
  no_h_x = fit_imm(paste("TREAT *", c_z, "+", h_z))
  data.frame(cytotoxic = imm_label(c_z), humoral = imm_label(h_z),
             rho = round(cor(d_imm[[c_z]], d_imm[[h_z]], method = "spearman"), 2),
             p_cytotoxic_interaction_given_humoral = lrt_imm(no_c_x, full),
             p_humoral_interaction_given_cytotoxic = lrt_imm(no_h_x, full))
}))

fmt_p(q1); fmt_p(q2); fmt_p(q3)

write_table(fmt_p(q1), file.path(output_dir, paste0("Immune_interactions_Q1_single_", imm_endpoint, ".csv")))
write_table(fmt_p(q2), file.path(output_dir, paste0("Immune_interactions_Q2_APM_adjusted_", imm_endpoint, ".csv")))
write_table(fmt_p(q3), file.path(output_dir, paste0("Immune_interactions_Q3_cytotoxic_vs_humoral_", imm_endpoint, ".csv")))


###### ----- Validation cohorts -----
###### ----- SCAN-B -----
###### SCAN-B: load data ######

# Load clinical, follow-up, expression and gene-annotation data.
# Retain operation samples with follow-up and a defined clinical subtype.

clinical_data_scanb = read_excel(file.path(data_dir, "SCAN-B", "Supplementary_Data_Table_1_2023-01-13.xlsx"), col_types = c(rep('guess', 39) , 'text', rep('guess', 47)))

# Keep patients with Follow-up data
clinical_data_scanb_selected = clinical_data_scanb %>% filter(Follow.up.cohort == TRUE)

# Check if there are any duplicated patients
sum(duplicated(clinical_data_scanb_selected$Patient))

# make ER and HER2 columns factors
clinical_data_scanb_selected = clinical_data_scanb_selected %>%
  dplyr::mutate(
    ER = factor(ER, levels = c("Negative", "Positive")),
    HER2 = factor(HER2, levels = c("Negative", "Positive")))

# Define clinical subtypes (using ER/HER2)
clinical_data_scanb_selected = clinical_data_scanb_selected %>%
  mutate(Clinical_Subtype = case_when(
    ER == 'Positive' & HER2 == 'Negative' ~ 'ER+/HER2-',
    HER2 == 'Positive' ~ 'HER2+',
    ER == 'Negative' & HER2 == 'Negative' ~ 'TNBC'))
clinical_data_scanb_selected$Clinical_Subtype = fct_relevel(clinical_data_scanb_selected$Clinical_Subtype, 'ER+/HER2-', 'HER2+', 'TNBC')

# Keep only non-NA values for clinical subtype
clinical_data_scanb_selected = clinical_data_scanb_selected %>% filter(!is.na(Clinical_Subtype))

# Keep only operation biopsy type
clinical_data_scanb_selected = clinical_data_scanb_selected %>% filter(BiopsyType == 'OP')

## read expression data
load(file.path(data_dir, "SCAN-B", "SCANB.9206.genematrix_noNeg.RData"))
gene_exp_scanb_selected = SCANB.9206.genematrix_noNeg[,clinical_data_scanb_selected$GEX.assay]

## read annotation data
load(file.path(data_dir, "SCAN-B", "Gene.ID.ann.RData"))

stopifnot(identical(rownames(SCANB.9206.genematrix_noNeg), Gene.ID.ann$Gene.ID))

# put gene symbols on the expression data
rownames(gene_exp_scanb_selected) = Gene.ID.ann$Gene.Name

### Data Pre-processing

# Derive nodal and grade categories, harmonize intrinsic subtype labels and rename overall-survival variables.

colnames(clinical_data_scanb_selected)[colnames(clinical_data_scanb_selected) == "GEX.assay"] = 'sample_id'

# Nodal category and grade, used in the adjusted Cox models of the KM plots
clinical_data_scanb_selected$LN.spec = factor(clinical_data_scanb_selected$LN.spec,
                                              levels = c("N0", "SubMicroMet", "1to3", "4toX"))

clinical_data_scanb_selected$Grade_cat = factor(case_when(clinical_data_scanb_selected$NHG == 1 ~ 'G1',
                                                          clinical_data_scanb_selected$NHG == 2 ~ 'G2',
                                                          clinical_data_scanb_selected$NHG == 3 ~ 'G3'),
                                                levels = c('G1', 'G2', 'G3'))

clinical_data_scanb_selected$SSP.Subtype = case_when(
  clinical_data_scanb_selected$SSP.Subtype == 'LumA' ~ 'Luminal A',
  clinical_data_scanb_selected$SSP.Subtype == 'LumB' ~ 'Luminal B',
  clinical_data_scanb_selected$SSP.Subtype == 'Her2' ~ 'HER2-enriched',
  clinical_data_scanb_selected$SSP.Subtype == 'Basal' ~ 'Basal-like',
  TRUE ~ NA)
clinical_data_scanb_selected$SSP.Subtype = factor(clinical_data_scanb_selected$SSP.Subtype, levels = c('Luminal A', 'Luminal B', 'HER2-enriched', 'Basal-like'))

names(clinical_data_scanb_selected)[names(clinical_data_scanb_selected) == "OS_days"] = "timeOS"
names(clinical_data_scanb_selected)[names(clinical_data_scanb_selected) == "OS_event"] = "OS.status"

### Criteria for PANTHER patients
# Approximate PANTHER eligibility using the available clinical variables.
# The inclusive age definition includes the ambiguous 31-35 age group.

age_var = "Age (5-year range, e.g., 35(31-35), 40(36-40), 45(41-45) etc.)"

scanb_eligibility = clinical_data_scanb_selected %>%
  mutate(
    age_group = as.numeric(.data[[age_var]]),
    
    # Available grouped age variable:
    # 25 = 21–25, ..., 65 = 61–65
    age_18_65 = age_group >= 25 & age_group <= 65,
    
    # Criterion 1: node-positive
    criterion_node_positive = LN == 1,
    
    # Criterion 2: node-negative, ER-negative, size >=20 mm, grade 3
    criterion_node_negative_ERneg = LN == 0 & ER == "Negative" & Size.mm >= 20 & NHG == 3,
    
    # Criterion 3: younger than 35 years
    # Inclusive approach: includes the ambiguous 31–35 age group
    criterion_age_under35 = age_group <= 35,
    
    eligible_inclusive = age_18_65 & (criterion_node_positive | criterion_node_negative_ERneg | criterion_age_under35))

### Normalization of gene expression data
# Restrict expression data to eligible patients and apply log2(expression + 1).
clinical_data_scanb_selected_panther_criteria = scanb_eligibility %>% filter(eligible_inclusive)
gene_exp_scanb_selected = gene_exp_scanb_selected[,clinical_data_scanb_selected_panther_criteria$sample_id]
gene_exp_scanb_selected = log2(gene_exp_scanb_selected + 1)

### Remove duplicated genes
# Retain the row with the highest mean log expression for each gene symbol.
n_rows_before = nrow(gene_exp_scanb_selected)
gene_exp_scanb_selected = collapse_max_mean(gene_exp_scanb_selected, rownames(gene_exp_scanb_selected))
message("SCAN-B: ", n_rows_before, " rows -> ", nrow(gene_exp_scanb_selected), " unique gene symbols")

###### APM score (HLA6 minus IFN5, frozen from PANTHER) -> APM_score ######

# Apply the PANTHER-selected signature, standardize scores within eligible SCAN-B patients and define cohort-median APM groups.

# The positive and negative gene sets fixed in PANTHER (CONFIG).
scanb_pair = resolve_signature_pair(CONFIG)

clinical_data_scanb_selected_panther_criteria = as.data.frame(clinical_data_scanb_selected_panther_criteria)

clinical_data_scanb_selected_panther_criteria$sample_id = as.character(clinical_data_scanb_selected_panther_criteria$sample_id)

scanb_ids = clinical_data_scanb_selected_panther_criteria$sample_id

stopifnot(!anyNA(scanb_ids), !anyDuplicated(scanb_ids), !anyDuplicated(colnames(gene_exp_scanb_selected)))

# Check that every signature gene is present.
signature_genes = c(scanb_pair$plus, scanb_pair$minus)

gene_matches = vapply(signature_genes, function(g) sum(rownames(gene_exp_scanb_selected) == g, na.rm = TRUE), integer(1))

if (any(gene_matches != 1L)) {
  stop("Missing or duplicated signature genes: ", paste(signature_genes[gene_matches != 1L], collapse = ", "))
}

sample_idx = match(scanb_ids, colnames(gene_exp_scanb_selected))

if (anyNA(sample_idx)) stop("Clinical samples missing from expression data.")

# Expression has already been log2-transformed above.
# GSVA uses the full gene background (every gene of the SCAN-B matrix).
expr_apm = as.matrix(gene_exp_scanb_selected[, sample_idx, drop = FALSE])

stopifnot(!anyDuplicated(rownames(expr_apm)), all(signature_genes %in% rownames(expr_apm)))

stopifnot(ncol(expr_apm) >= 3L, all(is.finite(expr_apm)))

scanb_score = score_signature(expr_apm, scanb_pair$plus, scanb_pair$minus, "SCAN-B APM score")

scanb_apm = data.frame(sample_id = names(scanb_score), APM_score = unname(scanb_score))

# Replace existing score/group columns when rerunning.
clinical_data_scanb_selected_panther_criteria =
  clinical_data_scanb_selected_panther_criteria %>%
  dplyr::select(-dplyr::any_of(c("APM_score", "APM_group"))) %>%
  dplyr::left_join(scanb_apm, by = "sample_id") %>%
  dplyr::mutate(APM_group = apm_group(APM_score))

###### APM group analysis ######
### Survival analysis
# Compare survival between APM groups for BCFi, RFi, DRFi and OS, displaying log-rank tests and adjusted Cox estimates.

for (endpoint in c('BCFi', 'RFi', 'DRFi', 'OS')) {
  # All (by APM score)
  km_all_apm = plot_km(dataset = clinical_data_scanb_selected_panther_criteria, predictor = "APM_group", predictor_name = 'APM score', cohort = 'SCAN-B', surv_endpoint = endpoint, colors = group_col_apm)
  
  if (endpoint == "BCFi") {
    Supp_Fig_6A = km_all_apm
  } else if (endpoint == "RFi") {
    Supp_Fig_6B = km_all_apm
  } else if (endpoint == "DRFi") {
    Supp_Fig_6C = km_all_apm
  } else if (endpoint == "OS") {
    Supp_Fig_6D = km_all_apm
  }
}

### APM score (Continuous) vs Clinical_Subtype

# Compare continuous APM scores across clinical subtypes using a Kruskal-Wallis test.

plot_df1 = clinical_data_scanb_selected_panther_criteria %>% filter(!is.na(APM_score), !is.na(Clinical_Subtype))

# Plot
p1 = ggplot(
  plot_df1,
  aes(
    x = Clinical_Subtype,
    y = APM_score,
    fill = Clinical_Subtype)) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.8) +
  
  stat_compare_means(method = "kruskal.test") +
  
  labs(
    title = "APM score across clinical subtypes",
    x = "Clinical subtype",
    y = "APM score") +
  
  theme_bw(base_size = 14) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1))

# Save plot
ggsave(filename = file.path(output_dir, "APM_vs_ClinicalSubtype_SCAN_B.png"), plot = p1, width = 8, height = 6, dpi = 300, bg = "white")

### APM score (Continuous) vs PAM50 subtype

# Compare continuous APM scores across SSP intrinsic subtypes using a Kruskal-Wallis test.

plot_df2 = clinical_data_scanb_selected_panther_criteria %>% filter(!is.na(APM_score), !is.na(SSP.Subtype))

# Plot
p2 = ggplot(
  plot_df2,
  aes(
    x = SSP.Subtype,
    y = APM_score,
    fill = SSP.Subtype)) +
  
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75) +
  
  geom_jitter(
    width = 0.15,
    alpha = 0.45,
    size = 1.8) +
  
  stat_compare_means(method = "kruskal.test") +
  
  labs(
    title = "APM score across intrinsic subtypes",
    x = "intrinsic subtype",
    y = "APM score") +
  
  theme_bw(base_size = 14) +
  
  theme(legend.position = "none")

# Save plot
ggsave(filename = file.path(output_dir, "APM_vs_PAM50Subtype_SCAN_B.png"), plot = p2, width = 8, height = 6, dpi = 300, bg = "white")

### APM group (median, Low vs. High) vs PAM50 subtype

# Display APM-group proportions within SSP intrinsic subtypes and test their association using Pearson's chi-square test.

plot_df4 = clinical_data_scanb_selected_panther_criteria %>%
  filter(
    !is.na(APM_group),
    !is.na(SSP.Subtype))

n_scanb_plot = nrow(plot_df4)

# Contingency table
table_pam50_cat = table(
  plot_df4$SSP.Subtype,
  plot_df4$APM_group)

# Chi-square test
chisq_pam50_cat = chisq.test(table_pam50_cat)

# Plot dataframe
plot_counts4 = plot_df4 %>%
  dplyr::count(SSP.Subtype, APM_group) %>%
  dplyr::group_by(SSP.Subtype) %>%
  dplyr::mutate(Percentage = 100 * n / sum(n))

# Plot
p4 = ggplot(
  plot_counts4,
  aes(
    x = SSP.Subtype,
    y = Percentage,
    fill = APM_group)) +
  
  geom_col(width = 0.75) +
  
  geom_text(
    aes(label = paste0(formatC(Percentage, format = "f", digits = 1), "%")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 4) +
  
  ggplot2::scale_fill_manual(
    values = group_col_apm,
    breaks = c("APM-low", "APM-high"),
    labels = c("APM-low", "APM-high"),
    drop = FALSE) +
  
  labs(
    title = "APM groups across intrinsic subtypes in the SCAN-B cohort",
    subtitle = paste0("N = ", n_scanb_plot, "; Chi-square ", p_lab(chisq_pam50_cat$p.value)),
    x = "Intrinsic subtype",
    y = "Percentage (%)",
    fill = "APM group") +
  
  theme_bw(base_size = 14)

# Save plot
ggsave(filename = file.path(output_dir, "APM_group_vs_PAM50Subtype_SCAN_B.png"), plot = p4, width = 8, height = 6, dpi = 300, bg = "white")
Supp_Figure_10B = p4

###### SCAN-B: MCPcounter ######

# Estimate immune and stromal population scores in eligible SCAN-B patients and merge them with clinical records.

MCPcounter = as.data.frame(t(MCPcounter.estimate(gene_exp_scanb_selected, featuresType = 'HUGO_symbols', genes = MCPcounter_gene)))
colnames(MCPcounter) = paste0('MCPcounter_', colnames(MCPcounter))
colnames(MCPcounter) = gsub(' ', '_', colnames(MCPcounter))
MCPcounter$sample_id = rownames(MCPcounter)

clinical_data_scanb_selected_panther_criteria = left_join(clinical_data_scanb_selected_panther_criteria, MCPcounter, by = 'sample_id')

###### SCAN-B: Danaher ######

# Calculate mean expression of available cell-type marker genes and merge the resulting scores with SCAN-B clinical records.

# Confirm marker availability in RNA-seq expression matrix
available_celltypegenes = intersect(rownames(gene_exp_scanb_selected), rownames(celltypematrix))
missing_celltypegenes = setdiff(rownames(celltypematrix), rownames(gene_exp_scanb_selected))

if (length(missing_celltypegenes) > 0) { warning("Some Danaher marker genes are absent from the expression matrix: ", paste(missing_celltypegenes, collapse = ", ")) }

# Restrict annotation matrix to genes available in expression data
celltypematrix_available = celltypematrix[available_celltypegenes, , drop = FALSE]

# Remove a cell type only if none of its marker genes is available
celltypematrix_available = celltypematrix_available[ , colSums(celltypematrix_available) > 0, drop = FALSE]

# Number of genes contributing to each immune population score
danaher_gene_counts = colSums(celltypematrix_available)

# Calculate Danaher immune population scores
# Expression matrix: samples × genes
df = t(gene_exp_scanb_selected[available_celltypegenes, , drop = FALSE])

# Weight matrix: genes × cell types
# Each marker receives equal weight within its immune population
celltype_weights = sweep(celltypematrix_available, 2, danaher_gene_counts, FUN = "/")

# Danaher score = mean expression of marker genes for each cell type
Danaher = as.data.frame(df %*% celltype_weights)

# Standardise output column names
colnames(Danaher) = paste0("Danaher_", colnames(Danaher))
colnames(Danaher) = gsub("-", "_", colnames(Danaher), fixed = TRUE)
colnames(Danaher) = gsub(" ", "_", colnames(Danaher), fixed = TRUE)
# Add sample identifier
Danaher$sample_id = rownames(Danaher)

clinical_data_scanb_selected_panther_criteria = left_join(clinical_data_scanb_selected_panther_criteria, Danaher, by = 'sample_id')

###### SCAN-B: boxplots for deconvolution methods ######

# Compare population scores between APM groups using Wilcoxon tests, with BH correction across populations within each method.

mcpcounter_results = plot_apm_immune_scores(
  data = clinical_data_scanb_selected_panther_criteria,
  score_prefix = "MCPcounter_",
  group_var = "APM_group",
  y_label = "MCP-counter score",
  title = NULL,
  feature_order = mcpcounter_order,
  display_p = "p.adj",
  p_adjust_method = "BH",
  y_multiplier = 1.30)

mcpcounter_results$plot

Supp_Figure_11B = mcpcounter_results$plot
ggsave(filename = file.path(output_dir, "ImmDeconv_MCPcounter_APM_groups_SCAN_B_boxplots.png"),
       plot = mcpcounter_results$plot, width = 12, height = 7, dpi = 300, bg = "white")

danaher_results = plot_apm_immune_scores(
  data = clinical_data_scanb_selected_panther_criteria,
  score_prefix = "Danaher_",
  group_var = "APM_group",
  y_label = "Danaher immune cell score",
  title = NULL,
  feature_order = danaher_order,
  display_p = "p.adj",
  p_adjust_method = "BH",
  y_multiplier = 1.30)

danaher_results$plot

Supp_Figure_12B = danaher_results$plot
ggsave(filename = file.path(output_dir, "ImmDeconv_Danaher_APM_groups_SCAN_B_boxplots.png"),
       plot = danaher_results$plot, width = 12, height = 7, dpi = 300, bg = "white")

###### ----- MATADOR trial -----
###### MATADOR: load data ######

# Harmonize clinical and TPM expression data, recover grade from older records, and identify patients with matched RNA-seq samples.

## Clinical data
MATADOR_clinical = read.csv(file.path(data_dir, "MATADOR", "20250401 CF32 MATADOR IHC manuscript Opdam 2025 clean survival data file.csv"), sep = ';', check.names = F)
MATADOR_clinical = MATADOR_clinical[,-1]
MATADOR_clinical[] = lapply(MATADOR_clinical, function(x) { if (is.character(x)) gsub(",", ".", x) else x })
MATADOR_clinical$sample_id = paste0("X", MATADOR_clinical$PATIENT_NO)

# Make columns numeric
cols = c("Age.at.randomization", "TILS.Mark", "RFSyear", "OSyear", "RFSevt", "OSevt", names(MATADOR_clinical)[which(names(MATADOR_clinical) == "ERmeanperc"):which(names(MATADOR_clinical) == "ThioRedmeancyto")])
MATADOR_clinical[cols] = lapply(MATADOR_clinical[cols], as.numeric)

# In the original MATADOR file, OSevt is coded as 1 = alive/censored and 2 = death/event.
# For consistency with the analysis pipeline, overwrite OSevt as: 0 = censored/alive, 1 = death/event.
MATADOR_clinical = MATADOR_clinical %>%
  dplyr::mutate(
    OSevt_original = OSevt,
    OSevt = dplyr::case_when(
      OSevt == 2 ~ 1,  # death/event
      OSevt == 1 ~ 0,  # alive/censored
      TRUE ~ NA_real_))

# Make columns factor
cols = c("Randomized.treatment", "pathsubtype", names(MATADOR_clinical)[which(names(MATADOR_clinical) == "Agebin"):which(names(MATADOR_clinical) == "Histological.tumor.type")])
MATADOR_clinical[cols] = lapply(MATADOR_clinical[cols], as.factor)
MATADOR_clinical$pathsubtype = factor(MATADOR_clinical$pathsubtype, levels = c("HR+HER2-", "HER2 positive", "triple negative"))

# Define ER and HER2 status for subgroup and adjusted Cox analyses.
# ER status is based on ERmeanperc using a 10% cutoff when available.
# When ERmeanperc is missing, infer ER status from HR+/HER2- or triple-negative subtype.
# HER2 status is extracted from the final pathological subtype variable.
MATADOR_clinical = MATADOR_clinical %>%
  dplyr::mutate(

    ER_status = dplyr::case_when(
      !is.na(ERmeanperc) & ERmeanperc >= 10 ~ "Positive",
      !is.na(ERmeanperc) & ERmeanperc < 10  ~ "Negative",

      # fallback when ERmeanperc is missing
      is.na(ERmeanperc) & pathsubtype == "HR+HER2-" ~ "Positive",
      is.na(ERmeanperc) & pathsubtype == "triple negative" ~ "Negative",
      TRUE ~ NA_character_),

    HER2_status = dplyr::case_when(
      pathsubtype == "HER2 positive" ~ "Positive",
      pathsubtype %in% c("HR+HER2-", "triple negative") ~ "Negative",
      TRUE ~ NA_character_),

    ER_status = factor(ER_status, levels = c("Negative", "Positive")),
    HER2_status = factor(HER2_status, levels = c("Negative", "Positive")))

# Load  the older clinical data file to get the full grade information (iScience, 2024)
MATADOR_clinical_old = read.csv(file.path(data_dir, "MATADOR", "20240507 CF32 MATADOR iscience mendeley data repository.csv"), check.names = F)
MATADOR_clinical_old = MATADOR_clinical_old[,-1]
colnames(MATADOR_clinical_old)[colnames(MATADOR_clinical_old) == "X"] = "PATIENT_NO"
MATADOR_clinical_old$grade = factor(MATADOR_clinical_old$grade, levels = c("good", "intermediate", "poor"))

# Add full grade info to the MATADOR clinical data
MATADOR_clinical = MATADOR_clinical %>% left_join(MATADOR_clinical_old %>% select(PATIENT_NO, grade), by = "PATIENT_NO")

MATADOR_clinical = MATADOR_clinical %>%
  dplyr::mutate(Randomized.treatment = dplyr::recode(Randomized.treatment, "6 x AC dd" = "ddAC", "6 x TAC" = "TAC"),
                Randomized.treatment = factor(Randomized.treatment, levels = c("TAC", "ddAC")))

## Expression data
# The file contains gene-level RNA-seq expression values reported as TPM (Transcripts Per Million), and the values are not log-transformed.
MATADOR_counts_TPM = read.delim(file.path(data_dir, "MATADOR", "GSE167977_20210223_matadorCounts.csv.gz"), row.names = 1, check.names = FALSE)

## Annotation from MATADOR Expression data
gene_annot_MATADOR = data.frame(
  GENEID = sub("_.*", "", rownames(MATADOR_counts_TPM)),
  SYMBOL = sub(".*_", "", rownames(MATADOR_counts_TPM)),
  row.names = rownames(MATADOR_counts_TPM),
  stringsAsFactors = FALSE)

rownames(MATADOR_counts_TPM) = gene_annot_MATADOR$GENEID

## Clinical data with RNA-seq samples
MATADOR_clinical_with_RNA = MATADOR_clinical %>% filter(sample_id %in% colnames(MATADOR_counts_TPM))

## TPM -> log2(TPM + 1)
tpm = log2(MATADOR_counts_TPM + 1)
MATADOR_tmm_normalized_tpm = as.data.frame(tpm)

## log2(TPM + 1) with gene symbols (no TMM step; object names follow the PANTHER pipeline)
MATADOR_tmm_normalized_tpm_symbols = MATADOR_tmm_normalized_tpm
MATADOR_tmm_normalized_tpm_symbols = MATADOR_tmm_normalized_tpm_symbols[order(rowSums(MATADOR_tmm_normalized_tpm_symbols), decreasing = TRUE),]
MATADOR_tmm_normalized_tpm_symbols$symbol = gene_annot_MATADOR[match(rownames(MATADOR_tmm_normalized_tpm_symbols), gene_annot_MATADOR$GENEID),]$SYMBOL
MATADOR_tmm_normalized_tpm_symbols = MATADOR_tmm_normalized_tpm_symbols[!duplicated(MATADOR_tmm_normalized_tpm_symbols$symbol),]
rownames(MATADOR_tmm_normalized_tpm_symbols) = MATADOR_tmm_normalized_tpm_symbols$symbol
MATADOR_tmm_normalized_tpm_symbols$symbol = NULL

## Reorder expression matrix columns to match clinical sample order
MATADOR_tmm_normalized_tpm_symbols = MATADOR_tmm_normalized_tpm_symbols[, match(MATADOR_clinical_with_RNA$sample_id[MATADOR_clinical_with_RNA$sample_id %in% colnames(MATADOR_tmm_normalized_tpm_symbols)], colnames(MATADOR_tmm_normalized_tpm_symbols)), drop = FALSE]
stopifnot(identical(colnames(MATADOR_tmm_normalized_tpm_symbols), MATADOR_clinical_with_RNA$sample_id[MATADOR_clinical_with_RNA$sample_id %in% colnames(MATADOR_tmm_normalized_tpm_symbols)]))

###### APM score (HLA6 minus IFN5, frozen from PANTHER) -> APM_score ######

# Uses CONFIG and the HLA6 / IFN5 gene vectors from the PANTHER APM signature section.
# Nothing is selected or tuned in MATADOR. The genes, the sign and the scoring method are
# the ones fixed in PANTHER (CONFIG; signed GSVA):
#   APM_score = z(GSVA(HLA6) - GSVA(IFN5))
#   GSVA: Gaussian kernel; background = the MATADOR expression matrix (log2(TPM + 1), all genes; GSVA drops constant genes).
#   z-standardization uses all RNA-matched MATADOR patients, before any endpoint or covariate exclusion.

matador_pair = resolve_signature_pair(CONFIG)   # positive = HLA6, negative = IFN5

# All RNA-matched MATADOR patients, in clinical order
matador_ids = as.character(MATADOR_clinical_with_RNA$sample_id)

matador_expr = as.matrix(MATADOR_tmm_normalized_tpm_symbols)

if (anyNA(matador_ids) || anyDuplicated(matador_ids) || !all(matador_ids %in% colnames(matador_expr)))
  stop("MATADOR sample IDs are missing, duplicated or not in the expression matrix.")

matador_expr = matador_expr[, matador_ids, drop = FALSE]

if (any(!is.finite(matador_expr))) stop("Non-finite MATADOR expression values.")

matador_missing_genes = setdiff(c(matador_pair$plus, matador_pair$minus), rownames(matador_expr))

if (length(matador_missing_genes))
  stop("Signature genes missing from MATADOR: ", paste(matador_missing_genes, collapse = ", "))

MATADOR_APM_SCORES = local({
  score = score_signature(matador_expr, matador_pair$plus, matador_pair$minus, "MATADOR APM score")
  out = data.frame(sample_id = names(score), Randomized.treatment = as.character(MATADOR_clinical_with_RNA$Randomized.treatment),
                   APM_score = unname(score), stringsAsFactors = FALSE)
  attr(out, "signature_label") = matador_pair$signature_label
  attr(out, "method") = "signed GSVA"
  out
})

# Add to the clinical table by sample ID (rerunning replaces the column)
matador_rows = match(matador_ids, MATADOR_APM_SCORES$sample_id)
stopifnot(!anyNA(matador_rows))
MATADOR_clinical_with_RNA$APM_score = MATADOR_APM_SCORES$APM_score[matador_rows]

matador_output_dir = file.path(output_dir, "MATADOR_APM")
dir.create(matador_output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(MATADOR_APM_SCORES, file.path(matador_output_dir, "MATADOR_APM_score.csv"), row.names = FALSE)

message(sprintf("MATADOR APM_score: %s, %s; %d RNA-matched patients",
                attr(MATADOR_APM_SCORES, "signature_label"), attr(MATADOR_APM_SCORES, "method"), nrow(MATADOR_APM_SCORES)))

###### Continuous APM score - HR (Forest plot) - APM score - Per arm ######

# APM_score per SD (z over all RNA-matched MATADOR patients) within each arm; the score x arm
# interaction (p interaction) is a 1-df LRT in the pooled adjusted model. Figure 2D (RFS), 2E (OS).
sig_cat = 'APM_score'

hr_column = paste(rep(' ', 40), collapse = ' ')

for (endpoint in c('RFS', 'OS')) {

  # Title row
  fonts = c('bold', 'plain', 'plain')
  row_cols = c('white', 'white', 'slategray1')

  results_mv_all_signatures = data.frame('  ' = character(), 'N' = numeric(),
                                         'est' =  numeric(), 'low' = numeric(), 'hi' = numeric(), 'se' = numeric(), ` ` = character(), `HR (95% CI)` = numeric(),
                                         'p value' = numeric(), check.names = FALSE)

  results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('  ' = 'Treatment arm', 'N' = '',
                                                                          'est' = NA, 'low' = NA, 'hi' = NA, 'se' = NA, ` ` = hr_column, `HR (95% CI)` = '',
                                                                          'p value' = '', check.names = FALSE))

  for (arm in c('TAC', 'ddAC')) {

    all_clinical_arm = MATADOR_clinical_with_RNA[MATADOR_clinical_with_RNA$Randomized.treatment == arm, ]

    data_cox = droplevels(all_clinical_arm)
    res.cox = coxph(as.formula(paste0('Surv(', endpoint, 'year, ', endpoint, 'evt) ~ ', sig_cat, ' + strata(ER_status) + HER2_status + TNM.after.surgery..T + TNM.after.surgery..N + grade')), data = data_cox)
    cox_model_results = tidy(res.cox, exponentiate = TRUE, conf.int = TRUE)

    for (i in seq_along(cox_model_results$term)) {
      term = cox_model_results$term[i]

      biomarker_label = if (term == sig_cat) { paste0(arm, " arm") } else { term }

      results_mv_all_signatures = rbind(results_mv_all_signatures, data.frame('  ' = biomarker_label,
                                                                              'N' = if (i == 1) res.cox$n else '',
                                                                              'est' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'estimate']),2),
                                                                              'low' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'conf.low']),2),
                                                                              'hi' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'conf.high']),2),
                                                                              'se' = round(as.numeric(cox_model_results[cox_model_results$term == term, 'std.error']),2),
                                                                              ` ` = hr_column,
                                                                              `HR (95% CI)` = format_hr_ci(as.numeric(cox_model_results[cox_model_results$term == term, 'estimate']), as.numeric(cox_model_results[cox_model_results$term == term, 'conf.low']), as.numeric(cox_model_results[cox_model_results$term == term, 'conf.high'])),
                                                                              'p value' = format_p_regular(as.numeric(cox_model_results[cox_model_results$term == term, 'p.value']), digits = 3), check.names = FALSE))
    }
  }

  results_mv_all_signatures_hide = results_mv_all_signatures[results_mv_all_signatures[['  ']] %in% c("Treatment arm", "TAC arm", "ddAC arm"), , drop = FALSE]

  # Pooled adjusted model without the interaction
  pint_formula = as.formula(paste0('Surv(', endpoint, 'year, ', endpoint, 'evt) ~ Randomized.treatment + ', sig_cat, ' + strata(ER_status) + HER2_status', ' + TNM.after.surgery..T + TNM.after.surgery..N + grade'))

  pint_keep = complete.cases(MATADOR_clinical_with_RNA[, all.vars(pint_formula), drop = FALSE]) & MATADOR_clinical_with_RNA$Randomized.treatment %in% c('TAC', 'ddAC')

  pint_data = droplevels(MATADOR_clinical_with_RNA[pint_keep, , drop = FALSE])

  pint_base = survival::coxph(pint_formula, data = pint_data, na.action = stats::na.fail)

  pint_full = survival::coxph(update(pint_formula, as.formula(paste0('. ~ . + Randomized.treatment:', sig_cat))), data = pint_data, na.action = stats::na.fail)

  pint_matador = lmtest::lrtest(pint_base, pint_full)$`Pr(>Chisq)`[2L]

  results_mv_all_signatures_hide[['p interaction']] = ''
  results_mv_all_signatures_hide[['p interaction']][results_mv_all_signatures_hide[['  ']] == 'TAC arm'] = format_p_regular(pint_matador, digits = 3)

  tm3 = forest_theme(core = list(fg_params = list(fontface = fonts),
                                 bg_params = list(fill = row_cols),
                                 padding = unit(c(12, 10), 'mm')),
                     # Reference line width/type/color
                     refline_gp = gpar(lwd = 1, lty = 'dashed', col = 'black'),
                     # Vertical line width/type/color
                     vertline_lwd = 1,
                     vertline_lty = 'dashed',
                     vertline_col = 'grey40')

  fp3 = forest(results_mv_all_signatures_hide[, c('  ', 'N', ' ', 'HR (95% CI)', 'p value', 'p interaction')],
               est = results_mv_all_signatures_hide$est,
               lower = results_mv_all_signatures_hide$low,
               upper = results_mv_all_signatures_hide$hi,
               sizes = 0.6,
               ci_column = 3,
               ref_line = 1,
               arrow_lab = c('Favors higher APM score', 'Favors lower APM score'),
               xlim = c(0, 3),
               ticks_at = c(0, 1, 2, 3),
               theme = tm3)

  if (endpoint == "RFS") {
    Fig_2D = fp3
  } else if (endpoint == "OS") {
    Fig_2E = fp3
  }
}

###### ----- TEX trial -----
###### TEX: load data ######

# Load relapse-biopsy clinical data and deposited microarray expression; 
# match expression arrays to patients using CEL barcodes; map probe sets to genes with the GEO platform annotation.
load(file.path(data_dir, "TEX", "pd.unique.mets.RData"))

clinical_data_tex = pd.unique.mets
rm(pd.unique.mets)

## Expression: GEO GSE56493 (GPL10379), log2 RMA as deposited (aroma.affymetrix).
## Downloaded once to tex_geo_dir and read from there afterwards.
tex_geo_dir = file.path(data_dir, "TEX", "GEO")
dir.create(tex_geo_dir, recursive = TRUE, showWarnings = FALSE)
tex_geo = GEOquery::getGEO("GSE56493", GSEMatrix = TRUE, getGPL = FALSE, destdir = tex_geo_dir)
tex_geo = tex_geo[[which(vapply(tex_geo, Biobase::annotation, character(1)) == "GPL10379")]]
exp_data = Biobase::exprs(tex_geo)   # probe sets x GSM arrays, log2

## Link GEO arrays to SampleIDRelapse by CEL barcode: GEO title "<n>_<barcode>", workbook BARCODE "@<barcode>"
tex_xlsx = file.path(data_dir, "TEX", "TEX_clinical incl microarray SEB.xlsx")

for (s in readxl::excel_sheets(tex_xlsx)) {
  tex_barcodes = as.data.frame(readxl::read_excel(tex_xlsx, sheet = s, col_types = "text"))
  if (all(c("SampleIDRelapse", "BARCODE") %in% names(tex_barcodes))) break
}

tex_barcodes = tex_barcodes[!is.na(tex_barcodes$SampleIDRelapse) & !is.na(tex_barcodes$BARCODE), c("SampleIDRelapse", "BARCODE")]
tex_gsm = colnames(exp_data)[match(sub("^@", "", trimws(tex_barcodes$BARCODE)), sub("^[0-9]+_", "", Biobase::pData(tex_geo)$title))]
stopifnot(!anyNA(tex_gsm), !anyDuplicated(tex_gsm), !anyDuplicated(tex_barcodes$SampleIDRelapse))
exp_data = exp_data[, tex_gsm, drop = FALSE]
colnames(exp_data) = trimws(tex_barcodes$SampleIDRelapse)

# Between-array quantile normalization (the GEO values are already log2: no log2 step)
# Rows = probes; columns = TEX samples
exp_data_norm = preprocessCore::normalize.quantiles(
  x = as.matrix(exp_data),
  copy = TRUE)

# Preserve probe IDs and sample names
dimnames(exp_data_norm) = dimnames(exp_data)

## Map probe sets -> gene symbols: GEO platform annotation GPL10379, column GeneSymbol as deposited
## (no symbol: "" or "---"). Downloaded once to tex_geo_dir and read from there afterwards.
tex_gpl = GEOquery::Table(GEOquery::getGEO("GPL10379", destdir = tex_geo_dir))
tex_probe_symbol = stats::setNames(trimws(as.character(tex_gpl$GeneSymbol)), as.character(tex_gpl$ID))
tex_probe_symbol[tex_probe_symbol %in% c("", "---")] = NA
stopifnot(!anyDuplicated(names(tex_probe_symbol)), all(rownames(exp_data_norm) %in% names(tex_probe_symbol)))

# put gene symbols as rownames
rownames(exp_data_norm) = unname(tex_probe_symbol[rownames(exp_data_norm)])

## Drop probes without symbol
keep = !is.na(rownames(exp_data_norm)) & rownames(exp_data_norm) != ""
exp_data_norm = exp_data_norm[keep, , drop = FALSE]

## Keep the probe set with the highest mean expression per gene (Miller et al. 2011, BMC Bioinformatics 12:322)
# exp_data_norm_iqr contains the highest-mean probe set for each gene.
mean_vec = rowMeans(exp_data_norm, na.rm = TRUE)
df = data.frame(gene = rownames(exp_data_norm), mean = mean_vec, stringsAsFactors = FALSE)
keep_idx = unlist(tapply(seq_len(nrow(df)), df$gene, function(idx) { idx[which.max(df$mean[idx])] }))
exp_data_norm_iqr = exp_data_norm[keep_idx, , drop = FALSE]

### Relapse biopsy site: seven-category definitions

# RelapseSite_clean has seven anatomical categories:
# Lymph node, Skin/chest wall, Breast, Liver, Lung, Bone and Ascites.
# Laterality and spelling variants are merged.
# All lymph-node regions are combined into one Lymph node category.
# Skin lesions over the breast remain Skin/chest wall.

# Cox models below use RelapseSite_clean as a categorical covariate.
# Lymph node is the reference category.
# Unrecognized nonmissing source labels stop the script.

site_map = c(   # key: lower case, single spaces  ->  tidy label
  'ascites'                          = 'Ascites',
  'lung'                             = 'Lung',
  'liver'                            = 'Liver',
  'bone costa 1'                     = 'Bone',
  'bone spine'                       = 'Bone',
  'bone thoracal spine'              = 'Bone',
  'bone, pelvis'                     = 'Bone',
  'bone, pubic'                      = 'Bone',
  'breast'                           = 'Breast',
  'breast right'                     = 'Breast',
  'left breast'                      = 'Breast',
  'right breast'                     = 'Breast',
  'axill. lymph node'                = 'Lymph node',
  'axillary lymphnode'               = 'Lymph node',
  'left axillary lymphnode'          = 'Lymph node',
  'lymph node axilla'                = 'Lymph node',
  'lymph node axilla dx'             = 'Lymph node',
  'lymph node left axilla'           = 'Lymph node',
  'lymph node right axilla'          = 'Lymph node',
  'lymph node fossa supclav sin'     = 'Lymph node',
  'lymph node fossa supraclav'       = 'Lymph node',
  'lymph node fossa supraclav dx'    = 'Lymph node',   # also catches the double-space spelling
  'lymph node fossa supraclav sin'   = 'Lymph node',
  'lymph node left fossa supraclav'  = 'Lymph node',
  'lymph node right fossa supraclav' = 'Lymph node',
  'lymph node supraclav right'       = 'Lymph node',
  'lymph node left fossa infraclav'  = 'Lymph node',
  'lymph node neck'                  = 'Lymph node',
  'lymph node left neck'             = 'Lymph node',
  'lymph node right neck'            = 'Lymph node',
  'lymph node inguinal dx'           = 'Lymph node',
  'lymph node, retrosternal'         = 'Lymph node',
  'lymph node'                       = 'Lymph node',
  'skin'                             = 'Skin/chest wall',
  'scar, skin'                       = 'Skin/chest wall',
  'skin left clavicle'               = 'Skin/chest wall',
  'skin parasternal right'           = 'Skin/chest wall',
  'skin presternal'                  = 'Skin/chest wall',
  'skin right breast (scar)'         = 'Skin/chest wall',
  'skin right shoulder'              = 'Skin/chest wall',
  'skin scar thoracic wall right'    = 'Skin/chest wall',
  'skin sternum'                     = 'Skin/chest wall',
  'skin substernal'                  = 'Skin/chest wall',
  'skin thoracic wall'               = 'Skin/chest wall',
  'skin thoracic wall right'         = 'Skin/chest wall',
  'skin, breast'                     = 'Skin/chest wall',
  'skin, right breast'               = 'Skin/chest wall',
  'skin, sternal'                    = 'Skin/chest wall',
  'skin, thoracic wall right'        = 'Skin/chest wall')

key = gsub('\\s+', ' ', trimws(tolower(as.character(clinical_data_tex$RelapseSite))))
unmapped = setdiff(unique(key[!is.na(key)]), names(site_map))
if (length(unmapped)) stop('Relapse sites without a mapping: ', paste(unmapped, collapse = '; '))

clean_levels = c('Lymph node', 'Skin/chest wall', 'Breast', 'Liver', 'Lung', 'Bone', 'Ascites')
stopifnot(setequal(clean_levels, unique(site_map)))

clinical_data_tex$RelapseSite_clean = factor(unname(site_map[key]), levels = clean_levels)

###### APM score (HLA6 minus IFN5, frozen from PANTHER) -> APM_score ######

# Uses CONFIG and the HLA6 / IFN5 gene vectors from the PANTHER APM signature section.
# Nothing is selected or tuned in TEX. The genes, the sign and the scoring method are the ones fixed in PANTHER (CONFIG; signed GSVA):
#   APM_score = z(GSVA(HLA6) - GSVA(IFN5))
#   Gene sets: the genes of each set measured on the TEX array (HLA-DQA1 is absent, so HLA6 has 5 of its 6 genes in TEX); the genes used are written to the output.
#   GSVA: Gaussian kernel; background = the TEX expression matrix (exp_data_norm_iqr, all genes; GSVA drops constant genes).
#   Scoring population: all expression-matched ET/TEX patients, before clinical exclusions; APM_score is z-standardized over these patients.
#   Survival population: scoring population with ER_PT, RelapseSite_clean and BCSS time/event available. 

tex_required_cols = c("SampleIDRelapse", "ClinTrialTreat", "ER_PT", "RelapseSite_clean", "RelToDeath_cens", "BCSS_cens")
tex_missing_cols = setdiff(tex_required_cols, names(clinical_data_tex))
if (length(tex_missing_cols)) stop("Missing TEX clinical columns: ", paste(tex_missing_cols, collapse = ", "))

# Adjustment used by the Cox models below
tex_adjust = "strata(ER_PT) + RelapseSite_clean"

## Gene sets: the PANTHER sets restricted to the genes on the TEX array
tex_expr = as.matrix(exp_data_norm_iqr)
tex_pair = resolve_signature_pair(CONFIG)   # positive = HLA6, negative = IFN5
tex_genes_defined = list(positive = tex_pair$plus, negative = tex_pair$minus)
tex_pair$plus = intersect(tex_pair$plus, rownames(tex_expr))
tex_pair$minus = intersect(tex_pair$minus, rownames(tex_expr))
if (!length(tex_pair$plus) || (tex_pair$has_negative && !length(tex_pair$minus)))
  stop("A signature gene set has no genes on the TEX array.")

tex_gene_coverage = data.frame(
  role = c("positive", "negative"),
  gene_set = c(tex_pair$positive_subgroup_id, tex_pair$negative_subgroup_id),
  n_defined = lengths(tex_genes_defined),
  n_used = c(length(tex_pair$plus), length(tex_pair$minus)),
  genes_used = c(paste(tex_pair$plus, collapse = ";"), paste(tex_pair$minus, collapse = ";")),
  genes_missing = c(paste(setdiff(tex_genes_defined$positive, tex_pair$plus), collapse = ";"),
                    paste(setdiff(tex_genes_defined$negative, tex_pair$minus), collapse = ";")),
  row.names = NULL, stringsAsFactors = FALSE)
print(tex_gene_coverage)

## Scoring population: all expression-matched ET/TEX patients
tex_clinical_ids = as.character(clinical_data_tex$SampleIDRelapse)
if (anyDuplicated(tex_clinical_ids[!is.na(tex_clinical_ids)]))
  stop("Duplicated SampleIDRelapse in clinical_data_tex.")
tex_scoring_rows = which(clinical_data_tex$ClinTrialTreat %in% c("ET", "TEX") & !is.na(tex_clinical_ids) & tex_clinical_ids %in% colnames(tex_expr))
sample_ids_tex = tex_clinical_ids[tex_scoring_rows]
if (length(sample_ids_tex) < 2L) stop("Fewer than two expression-matched ET/TEX samples.")
tex_expr = tex_expr[, sample_ids_tex, drop = FALSE]
if (any(!is.finite(tex_expr))) stop("Non-finite TEX expression values.")

tex_score = score_signature(tex_expr, tex_pair$plus, tex_pair$minus, "TEX APM score")

TEX_APM_RESULT = list(
  signature_id = tex_pair$signature_id,
  signature_label = tex_pair$signature_label,
  method = "signed GSVA",
  gene_coverage = tex_gene_coverage,
  n_scoring_population = length(sample_ids_tex))

# APM_score for every scored patient (NA for the others); rerunning replaces the column
clinical_data_tex$APM_score = unname(tex_score[match(tex_clinical_ids, names(tex_score))])

## Survival population (ET/TEX, expression-matched, ER, relapse site and BCSS available)
clinical_data_tex_selected = clinical_data_tex[clinical_data_tex$ClinTrialTreat %in% c("ET", "TEX"), , drop = FALSE]
clinical_data_tex_selected$ClinTrialTreat = factor(clinical_data_tex_selected$ClinTrialTreat, levels = c("ET", "TEX"))

clinical_data_tex_selected_er = dplyr::filter(
  clinical_data_tex_selected,
  as.character(SampleIDRelapse) %in% sample_ids_tex,
  !is.na(ER_PT),
  !is.na(RelapseSite_clean),
  !is.na(RelToDeath_cens),
  !is.na(BCSS_cens),
  is.finite(APM_score))
stopifnot(nrow(clinical_data_tex_selected_er) > 0L)

# Relabel ER
tex_er = as.character(clinical_data_tex_selected_er$ER_PT)
if (any(!tex_er %in% c("0", "1"))) stop("ER_PT must be coded 0/1 before relabeling.")
clinical_data_tex_selected_er$ER_PT = factor(tex_er, levels = c("0", "1"), labels = c("ER−", "ER+"))

# Remove unused site levels while preserving the Lymph node reference
clinical_data_tex_selected_er$RelapseSite_clean = droplevels(clinical_data_tex_selected_er$RelapseSite_clean)

TEX_APM_RESULT$n_survival_population = nrow(clinical_data_tex_selected_er)

tex_output_dir = file.path(output_dir, "TEX_APM")
dir.create(tex_output_dir, recursive = TRUE, showWarnings = FALSE)
TEX_APM_RESULT$output_dir = tex_output_dir

utils::write.csv(tex_gene_coverage, file.path(tex_output_dir, "TEX_APM_signature_genes.csv"), row.names = FALSE)
utils::write.csv(data.frame(SampleIDRelapse = names(tex_score),
                            ClinTrialTreat = as.character(clinical_data_tex$ClinTrialTreat[tex_scoring_rows]),
                            APM_score = unname(tex_score),
                            survival_population = names(tex_score) %in% as.character(clinical_data_tex_selected_er$SampleIDRelapse)),
                 file.path(tex_output_dir, "TEX_APM_score.csv"), row.names = FALSE)

message(sprintf("TEX APM_score: %s, signed GSVA; %d scored ET/TEX patients, %d in the survival population",
                TEX_APM_RESULT$signature_label, length(sample_ids_tex), nrow(clinical_data_tex_selected_er)))

###### Continuous APM score - HR (Forest plot) - APM score - Per arm ######

# BCSS model: treatment * continuous APM_score + strata(ER_PT) + RelapseSite_clean (the same patients throughout).
# Score HRs from separate adjusted ET and TEX models; the interaction p is from the pooled model (LRT).
# Figure 2C: the arm rows in the PANTHER per-arm layout.

TEX_APM_CONTINUOUS_RESULT = local({

  ###### 1. Same data and adjustment; continuous APM score ######

  data_cont = as.data.frame(clinical_data_tex_selected_er)
  required = c("RelToDeath_cens", "BCSS_cens", "ClinTrialTreat", "ER_PT", "RelapseSite_clean", "APM_score")

  if (!all(required %in% names(data_cont)))
    stop("Missing columns: ", paste(setdiff(required, names(data_cont)), collapse = ", "))

  if (!is.numeric(data_cont$APM_score) || any(!is.na(data_cont$APM_score) & !is.finite(data_cont$APM_score)))
    stop("APM_score must contain finite numeric values or NA.")

  # APM_score already contains the score_signature()’s final cohort-standardized score.
  # Preserve its original RNA-matched ET/TEX SD in pooled and independent models.
  # The following reference statistics are an audit, not a second transformation.
  score_reference = tex_score
  score_reference = score_reference[is.finite(score_reference)]
  score_center = mean(score_reference)
  score_SD = stats::sd(score_reference)
  if (!is.finite(score_SD) || score_SD <= 0) stop("APM_score is invariant or has too few values.")
  if (abs(score_center) > 1e-6 || abs(score_SD - 1) > 1e-6)
    stop("Expected the shared helper's cohort-standardized APM_score; rerun TEX scoring before this block.")

  n_input = nrow(data_cont)
  complete = stats::complete.cases(data_cont[, required, drop = FALSE]) & as.character(data_cont$ClinTrialTreat) %in% c("ET", "TEX")
  data_cont = data_cont[complete, , drop = FALSE]
  data_cont$ClinTrialTreat = factor(as.character(data_cont$ClinTrialTreat), levels = c("ET", "TEX"))
  data_cont$ER_PT = droplevels(factor(data_cont$ER_PT))
  data_cont$RelapseSite_clean = droplevels(factor(data_cont$RelapseSite_clean))

  if (!nrow(data_cont) || any(table(data_cont$ClinTrialTreat) == 0L))
    stop("Both ET and TEX must be present in the complete-case cohort.")

  if (!is.numeric(data_cont$RelToDeath_cens) || any(!is.finite(data_cont$RelToDeath_cens) | data_cont$RelToDeath_cens < 0) || !all(data_cont$BCSS_cens %in% c(0, 1)))
    stop("Check supplied survival times and 0/1 BCSS_cens coding; no recoding was applied.")

  # Explicit ET reference, independent of the session's contrasts option.
  contrasts(data_cont$ClinTrialTreat) = stats::contr.treatment(c("ET", "TEX"), base = 1L)
  strata = survival::strata

  res.cox_cont_base = survival::coxph(survival::Surv(RelToDeath_cens, BCSS_cens) ~ ClinTrialTreat + APM_score + strata(ER_PT) + RelapseSite_clean,
                                      data = data_cont, model = TRUE, x = TRUE, na.action = stats::na.fail)

  res.cox_cont_int = survival::coxph(survival::Surv(RelToDeath_cens, BCSS_cens) ~ ClinTrialTreat * APM_score + strata(ER_PT) + RelapseSite_clean,
                                     data = data_cont, model = TRUE, x = TRUE, na.action = stats::na.fail)

  stopifnot(identical(rownames(res.cox_cont_base$model), rownames(res.cox_cont_int$model)))
  mv_cont_interaction_lr = lmtest::lrtest(res.cox_cont_base, res.cox_cont_int)
  mv_cont_interaction_p_value = mv_cont_interaction_lr$`Pr(>Chisq)`[2L]

  ###### 2. Independent adjusted models within ET and TEX ######

  # Use exactly the same complete-case patients and unchanged cohort-wide score.
  # Each arm has its own baseline hazards and relapse-site coefficients.
  # A reduced model is fitted only to provide the within-arm score LRT.
  independent_arm_models = stats::setNames(lapply(c("ET", "TEX"), function(arm_name) {
    d_arm = droplevels(data_cont[data_cont$ClinTrialTreat == arm_name, , drop = FALSE])

    if (nrow(d_arm) < 2L || sum(d_arm$BCSS_cens) == 0L ||
        !is.finite(stats::sd(d_arm$APM_score)) || stats::sd(d_arm$APM_score) <= 0)
      stop(arm_name, ": insufficient patients/events or invariant APM_score for an independent arm model.")

    if (nlevels(d_arm$RelapseSite_clean) < 2L)
      stop(arm_name, ": relapse site has fewer than two observed levels; inspect before fitting the specified adjusted model.")

    full = survival::coxph(survival::Surv(RelToDeath_cens, BCSS_cens) ~ APM_score + strata(ER_PT) + RelapseSite_clean, data = d_arm, model = TRUE, x = TRUE, na.action = stats::na.fail, singular.ok = FALSE)
    reduced = survival::coxph(survival::Surv(RelToDeath_cens, BCSS_cens) ~ strata(ER_PT) + RelapseSite_clean, data = d_arm, model = TRUE, x = TRUE, na.action = stats::na.fail, singular.ok = FALSE)
    stopifnot(identical(rownames(full$model), rownames(d_arm)), identical(rownames(full$model), rownames(reduced$model)))

    if (any(!is.finite(stats::coef(full))) || any(!is.finite(stats::vcov(full))) || any(!is.finite(stats::coef(reduced))) || any(!is.finite(stats::vcov(reduced))))
      stop(arm_name, ": non-estimable independent arm model; inspect the Cox warnings.")

    sm = summary(full)
    coef_score = sm$coefficients["APM_score", ]
    ci_score = sm$conf.int["APM_score", ]
    score_lr = lmtest::lrtest(reduced, full)
    result = data.frame(
      arm = arm_name, n = nrow(d_arm), events = sum(d_arm$BCSS_cens),
      beta = unname(coef_score["coef"]), se = unname(coef_score["se(coef)"]),
      est = unname(ci_score["exp(coef)"]),
      low = unname(ci_score["lower .95"]), hi = unname(ci_score["upper .95"]),
      p_Wald = unname(coef_score["Pr(>|z|)"]), p_LRT = score_lr$`Pr(>Chisq)`[2L],
      model_source = "Independent model within each arm", stringsAsFactors = FALSE)

    if (any(!is.finite(as.matrix(result[, c("beta", "se", "est", "low", "hi", "p_Wald", "p_LRT")]))))
      stop(arm_name, ": nonfinite independent-arm results; inspect the Cox warnings.")

    list(full = full, reduced = reduced, results = result)
  }), c("ET", "TEX"))
  per_arm_independent_results = do.call(rbind, lapply(independent_arm_models, `[[`, "results"))
  rownames(per_arm_independent_results) = NULL
  
  ###### 3. Figure 2C: independent-arm score HRs in the PANTHER per-arm layout ######
  
  hr_column = paste(rep(" ", 30), collapse = " ")

  # Score HRs from the independent ET and TEX models
  plot_cont_per_arm_independent = rbind(
    data.frame(
      '  ' = paste0(per_arm_independent_results$arm, ' arm'),
      'N' = per_arm_independent_results$n,
      'est' = per_arm_independent_results$est,
      'low' = per_arm_independent_results$low,
      'hi' = per_arm_independent_results$hi,
      'se' = per_arm_independent_results$se,
      ` ` = rep(hr_column, nrow(per_arm_independent_results)),
      `HR (95% CI)` = format_hr_ci(
        per_arm_independent_results$est,
        per_arm_independent_results$low,
        per_arm_independent_results$hi),
      'p value' = format_p_regular(per_arm_independent_results$p_Wald), check.names = FALSE))

  # PANTHER-style lower panel: heading followed by the two arm rows
  plot_cont_independent_arms = rbind(
    data.frame(
      '  ' = 'Treatment arm',
      'N' = '',
      'est' = NA_real_,
      'low' = NA_real_,
      'hi' = NA_real_,
      'se' = NA_real_,
      ` ` = hr_column,
      `HR (95% CI)` = '',
      'p value' = '',
      check.names = FALSE),
    plot_cont_per_arm_independent)

  plot_cont_independent_arms[['p interaction']] = ''

  plot_cont_independent_arms[['p interaction']][plot_cont_independent_arms[['  ']] == 'ET arm'] = format_p_regular(mv_cont_interaction_p_value, digits = 3)

  display_columns = c('  ', 'N', ' ', 'HR (95% CI)', 'p value', 'p interaction')

  # Lower panel: same layout as the PANTHER per-arm forest
  tm_cont3_2 = forestploter::forest_theme(
    core = list(
      fg_params = list(fontface = c('bold', 'plain', 'plain')),
      bg_params = list(fill = c('white', 'white', 'slategray1')),
      padding = unit(c(12, 10), 'mm')),
    refline_gp = grid::gpar(lwd = 1, lty = 'dashed', col = 'black'))

  fp_cont3_2 = forestploter::forest(
    plot_cont_independent_arms[, display_columns, drop = FALSE],
    est = plot_cont_independent_arms$est,
    lower = plot_cont_independent_arms$low,
    upper = plot_cont_independent_arms$hi,
    sizes = 0.6,
    ci_column = 3,
    ref_line = 1,
    arrow_lab = c('Favors higher APM score', 'Favors lower APM score'),
    xlim = c(0, 3),
    ticks_at = c(0, 1, 2, 3),
    theme = tm_cont3_2)

  # Figure 2C: the per-arm score panel
  Fig_2C = fp_cont3_2

  ###### 4. Save the figures and tables ######

  cont_output_dir = file.path(tex_output_dir, "continuous_APM_forest")
  dir.create(cont_output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(cont_output_dir)) stop("Could not create the output directory.")

  settings = data.frame(
    item = c("endpoint", "adjustment", "score_unit", "score_reference_population", "n_score_reference", "score_center", "score_SD", "n_input", "n_model", "events_model"),
    value = as.character(c("BCSS: supplied time/status", "strata(ER_PT) + RelapseSite_clean", "1 original RNA-matched ET/TEX cohort SD higher APM_score; no rescaling", "All RNA-matched ET/TEX patients used by the shared scoring helper before clinical exclusions",
                           length(score_reference), score_center, score_SD, n_input, nrow(data_cont), sum(data_cont$BCSS_cens))))

  utils::write.csv(per_arm_independent_results, file.path(cont_output_dir, "TEX_APM_continuous_per_arm_independent_results.csv"), row.names = FALSE)
  utils::write.csv(settings, file.path(cont_output_dir, "TEX_APM_continuous_settings.csv"), row.names = FALSE)

  # Save summary results only; model frames remain in the R workspace.
  saveRDS(list(per_arm_independent = per_arm_independent_results, settings = settings),
          file.path(cont_output_dir, "TEX_APM_continuous_results.rds"))

  print(per_arm_independent_results)

  list(per_arm_independent = per_arm_independent_results, settings = settings,
       models = list(additive = res.cox_cont_base, interaction = res.cox_cont_int, independent_arms = independent_arm_models),
       LRT = mv_cont_interaction_lr,
       display_per_arm_independent = plot_cont_per_arm_independent, output_dir = cont_output_dir,
       plots = list(Fig_2C = Fig_2C))
})

Fig_2C = TEX_APM_CONTINUOUS_RESULT$plots$Fig_2C

TEX_APM_RESULT$adjustment = tex_adjust
TEX_APM_RESULT$completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

saveRDS(TEX_APM_RESULT, file.path(tex_output_dir, "TEX_APM_RESULT.rds"))

capture.output(sessionInfo(), file = file.path(tex_output_dir, "TEX_sessionInfo.txt"))

message("TEX results saved in: ", tex_output_dir)

###### ----- Median Follow-up -----

# Estimate endpoint-specific median follow-up using reverse Kaplan-Meier for PANTHER, SCAN-B and MATADOR.

## PANTHER
endpoints = c("BCRFS", "EFS", "DDFS", "OS")
followup_summary = purrr::map_dfr(endpoints, ~ median_followup_reverse_km(all_clinical, .x, time_scale = "days"))

followup_summary

## SCAN-B
scanb_endpoints = tibble::tribble(
  ~endpoint, ~time_var,    ~status_var,
  "BCFi", "BCFi_days",  "BCFi_event",
  "RFi",  "RFi_days",   "RFi_event",
  "DRFi", "DRFi_days",  "DRFi_event",
  "OS",   "timeOS",     "OS.status")

scanb_followup_summary = purrr::pmap_dfr(
  scanb_endpoints,
  ~ median_followup_reverse_km(
    data = clinical_data_scanb_selected_panther_criteria,
    endpoint = ..1,
    time_var = ..2,
    status_var = ..3,
    time_scale = "days"))

scanb_followup_summary

## MATADOR
matador_endpoints = tibble::tribble(
  ~endpoint, ~time_var,  ~status_var,
  "RFS",     "RFSyear",  "RFSevt",
  "OS",      "OSyear",   "OSevt")

matador_followup_summary = purrr::pmap_dfr(
  matador_endpoints,
  ~ median_followup_reverse_km(
    data = MATADOR_clinical_with_RNA,
    endpoint = ..1,
    time_var = ..2,
    status_var = ..3,
    time_scale = "years"))

matador_followup_summary

###### ----- Main Figures -----
###### Figure 1 : CONSORT diagram ######

# Placeholder for the PANTHER CONSORT diagram; The diagram was created in Power Point.

###### Figure 2 : Forest plots (PANTHER, TEX, MATADOR) ######

# Assemble PANTHER median-group treatment effects (A), PANTHER and TEX continuous-score effects (B-C), and MATADOR RFS/OS score effects (D-E).

fig2_file = file.path(output_dir, "Figure_2.pdf")

# 1. Same column widths in B-E (same 6-column layout); A scaled to the same total width
w_BE = unit.pmax(Fig_2B$widths, Fig_2C$widths, Fig_2D$widths, Fig_2E$widths)
Fig_2B$widths = w_BE; Fig_2C$widths = w_BE; Fig_2D$widths = w_BE; Fig_2E$widths = w_BE
Fig_2A$widths = Fig_2A$widths * (in_cm(w_BE) / in_cm(Fig_2A$widths))

# 2. Each panel gets its letter at the top-left and a small gap below
panels = Map(label_panel, list(Fig_2A, Fig_2B, Fig_2C, Fig_2D, Fig_2E), LETTERS[1:5])

# 3. Stack the panels at their natural heights; the page size follows the content
heights = do.call(unit.c, lapply(panels, function(g) sum(g$heights)))
figure_2 = arrangeGrob(grobs = panels, ncol = 1, heights = heights)

page_w = convertWidth(sum(panels[[2]]$widths), "in", valueOnly = TRUE) + 0.4
page_h = convertHeight(sum(heights), "in", valueOnly = TRUE) + 0.4

pdf(fig2_file, width = page_w, height = page_h)
grid.draw(figure_2)
dev.off()

###### Figure 3 : KM curves for BCRFS (APM stratification per arm) ######

# Arrange treatment-arm BCRFS curves for APM-high and APM-low patients, including risk tables.

Fig_3A$plot = Fig_3A$plot +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Fig_3B$plot = Fig_3B$plot +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Figure_3 = survminer::arrange_ggsurvplots(
  x = list(Fig_3A, Fig_3B),
  nrow = 1,
  ncol = 2,
  print = FALSE,
  surv.plot.height = 0.78,
  risk.table.height = 0.22)

ggplot2::ggsave(filename = file.path(output_dir, "Figure_3.pdf"),
                plot = Figure_3, width = 15, height = 6, units = "in", device = grDevices::cairo_pdf, bg = "white")

###### Figure 4 : Complex Heatmap ######

# Plot gene-wise standardized expression for proliferation, ER and IgG modules, split by APM group and treatment, with clinical annotations.
# Shared genes appear once, assigned to the first listed module.

rownames(all_clinical) = all_clinical$sample_id

prolif_genes_nielsen = c('BIRC5', 'CCNB1', 'CDC20', 'NUF2', 'CEP55', 'NDC80', 'MKI67', 'PTTG1', 'RRM2', 'TYMS', 'UBE2C')
prolif_genes_wallden = c('ANLN', 'CCNE1', 'CDC20', 'CDC6', 'NUF2', 'CENPF', 'CEP55', 'EXO1', 'KIF2C', 'NDC80', 'MELK', 'MKI67', 'ORC6', 'PTTG1', 'RRM2', 'TYMS', 'UBE2C', 'UBE2T')
prolif_oncotype_genes = c('MKI67', 'AURKA', 'BIRC5', 'CCNB1', 'MYBL2')
er_oncotype_genes = c('ESR1', 'PGR', 'BCL2', 'SCUBE2')
er_set_genes = c('ABAT', 'ADCY1', 'AZGP1', 'CA12', 'CD2', 'CD3D', 'DNAJC12', 'ESR1', 'KCNE4', 'MAPT', 'MRPS30', 'NAT1', 'NPY1R', 'PDZK1', 'QDPR', 'SCUBE2', 'SLC39A6', 'STC2')
igg_genes = immune_igg_genes   # 14-gene IgG signature (IgG section)

heatmap_gene_sets = list(
  "11-gene Proliferation"    = prolif_genes_nielsen,
  "18-gene Proliferation"    = prolif_genes_wallden,
  "21-gene RS Proliferation" = prolif_oncotype_genes,
  "21-gene RS Estrogen"      = er_oncotype_genes,
  "18-gene SET ER/PR"        = er_set_genes,
  "14-gene IgG"              = igg_genes)

heatmap_gene_membership = tibble::tibble(gene = unlist(heatmap_gene_sets, use.names = FALSE), module = rep(names(heatmap_gene_sets), lengths(heatmap_gene_sets)))

overlapping_genes_heatmap = heatmap_gene_membership %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(n_modules = dplyr::n_distinct(module), modules = paste(module, collapse = " | "), .groups = "drop") %>%
  dplyr::filter(n_modules > 1) %>%
  dplyr::arrange(dplyr::desc(n_modules), gene)

overlapping_genes_heatmap

heatmap_gene_membership_unique = heatmap_gene_membership %>%
  dplyr::filter(!duplicated(gene)) %>%
  dplyr::filter(gene %in% rownames(tmm_normalized_tpm_symbols))

heatmap_gene_exp_matrix = as.matrix(tmm_normalized_tpm_symbols[heatmap_gene_membership_unique$gene, rownames(all_clinical), drop = FALSE])
heatmap_gene_exp_matrix = t(scale(t(heatmap_gene_exp_matrix)))

heatmap_genes_groups = heatmap_gene_membership_unique$module

names(heatmap_genes_groups) = heatmap_gene_membership_unique$gene

heatmap_genes_groups = factor(heatmap_genes_groups[rownames(heatmap_gene_exp_matrix)], levels = names(heatmap_gene_sets))

top_annot = data.frame('Group' = as.factor(all_clinical$APM_group),
                       'Treatment arm' = as.factor(all_clinical$TREAT),
                       'Size group' = all_clinical$tumdiam_cat,
                       'Nodal status' = as.factor(all_clinical$posnodes_cat),
                       'Histological grade' = all_clinical$histgrade,
                       'Ki-67' = all_clinical$prolifvalue,
                       'Clinical subtypes' = all_clinical$Clinical_Subtype,
                       'Intrinsic subtypes' = all_clinical$subtypecd_ssp.subtype,
                       check.names = F)
rownames(top_annot) = rownames(all_clinical)
stopifnot(identical(rownames(top_annot), colnames(heatmap_gene_exp_matrix)))

combined_split = factor(paste0(top_annot$Group, ' \n ', top_annot$`Treatment arm`),
                        levels = c('APM-high \n FEC-D', 'APM-high \n tddEC-D', 'APM-low \n FEC-D', 'APM-low \n tddEC-D'))

top_annot_left = HeatmapAnnotation(df = top_annot,
                                   col = list(Group = group_col_apm,
                                              'Treatment arm' = c('FEC-D' = 'dodgerblue3' , 'tddEC-D' = 'maroon3'),
                                              'Size group' = c('≤ 2 cm' = 'skyblue1', '> 2 cm' = 'skyblue4'),
                                              'Nodal status' = c('N0' = 'seagreen3', 'N1-3' = "mistyrose", 'N4-9' = "orangered3", 'N10+' = "darkred"),
                                              'Histological grade' = c('I' = 'tan1', 'II' = 'tan3', 'III' = 'tan4'),
                                              'Ki-67' = colorRamp2(c(0, 100), c('gold', 'gold4')),
                                              'Clinical subtypes' = c('ER+/HER2-' = '#3769B1', 'HER2+' = '#CC79A7', 'TNBC' = '#7F3B08'),
                                              'Intrinsic subtypes' = c('Luminal A' = '#2a3188', 'Luminal B' = '#419ad2', 'HER2-enriched' = '#d4279c', 'Basal-like' = '#97191e')),
                                   # Annotation row names
                                   annotation_name_gp = gpar(fontsize = 10),
                                   annotation_name_side = 'left',
                                   annotation_legend_param = list(
                                     # Legend titles
                                     title_gp = gpar(fontsize = 8, fontface = "bold"),
                                     # Legend scale labels
                                     labels_gp = gpar(fontsize = 8),
                                     grid_height = unit(2, "mm"),
                                     grid_width = unit(4, "mm")))

cairo_pdf(file.path(output_dir, "Figure_4.pdf"), width = 10, height = 12)
ComplexHeatmap::draw(
  Heatmap(
    heatmap_gene_exp_matrix,
    col = colorRamp2(c(-2, 0, 2), c('#3367A1', 'white', '#A2222D')),
    column_split = combined_split,
    row_split = heatmap_genes_groups,
    cluster_row_slices = FALSE,
    cluster_column_slices = FALSE,
    top_annotation = top_annot_left,
    row_labels = rownames(heatmap_gene_exp_matrix),
    column_labels = rep('', ncol(heatmap_gene_exp_matrix)),
    # gene symbols
    row_names_gp = gpar(fontsize = 8),
    # APM group names
    column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
    # gene signature groups
    row_title_gp = gpar(fontsize = 10, fontface = 'bold'),
    row_title_rot = 0,
    heatmap_legend_param = list(
      title = 'Z-score',
      # Z-score title
      title_gp = grid::gpar(fontsize = 8, fontface = "bold"),
      # Z-score scale labels
      labels_gp = grid::gpar(fontsize = 8),
      legend_height = unit(15, "mm"))),
  # Configure heatmap drawing and legend placement.
  merge_legends = TRUE,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  align_heatmap_legend = "heatmap_top",
  align_annotation_legend = "heatmap_top",
  legend_gap = unit(4, "mm"),
  padding = unit(c(2, 2, 8, 2), 'mm'))
dev.off()

###### Complementary to Figure 4 ######

# Compare mean module expression between APM groups in PANTHER and SCAN-B.
# Return nominal and BH-adjusted Wilcoxon p-values; plots display nominal p-values.

## Gene expression matrix
## rows = genes, columns = samples
expr_mat = tmm_normalized_tpm_symbols

## Phenotype / clinical data
pheno_df = all_clinical

## Define gene modules
gene_modules = list(
  Proliferation_Nielsen  = prolif_genes_nielsen,
  Proliferation_Wallden  = prolif_genes_wallden,
  Proliferation_Oncotype = prolif_oncotype_genes,
  ER_Oncotype            = er_oncotype_genes,
  ER_SET                 = er_set_genes,
  IgG_14gene             = igg_genes)

res_panther_modules = run_module_score_apm_comparison(
  expr_mat = tmm_normalized_tpm_symbols,
  pheno_df = all_clinical,
  gene_modules = gene_modules,
  cohort_name = "PANTHER",
  sample_id_col = "sample_id",
  apm_group_col = "APM_group")

res_panther_modules$wilcox_results
res_panther_modules$plot

res_scanb_modules = run_module_score_apm_comparison(
  expr_mat = gene_exp_scanb_selected,
  pheno_df = clinical_data_scanb_selected_panther_criteria,
  gene_modules = gene_modules,
  cohort_name = "SCAN-B",
  sample_id_col = "sample_id",
  apm_group_col = "APM_group")

res_scanb_modules$wilcox_results
res_scanb_modules$plot

ggsave(filename = file.path(output_dir, "PANTHER_module_scores_by_APM_group.png"),
       plot = res_panther_modules$plot, width = 10, height = 7, dpi = 300, bg = "white")

ggsave(filename = file.path(output_dir, "SCANB_module_scores_by_APM_group.png"),
       plot = res_scanb_modules$plot, width = 10, height = 7, dpi = 300, bg = "white")

###### ----- Supplementary Figures -----
###### Supp. Figure 1 : TEX CONSORT diagram ######

# Placeholder for the TEX CONSORT diagram; The diagram was created in Power Point.

###### Supp. Figure 2 : SET ######

# Assemble continuous treatment-effect plots for SET2,3 (A-D) and SET ER/PR (E-H) across four endpoints in ER+/HER2- patients.

Supp_Figure_2_1 = gridExtra::arrangeGrob(
  grobs = list(
    gridExtra::arrangeGrob(Supp_Figure_2A, top = grid::textGrob("A", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16))),
    gridExtra::arrangeGrob(Supp_Figure_2B, top = grid::textGrob("B", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16))),
    gridExtra::arrangeGrob(Supp_Figure_2C, top = grid::textGrob("C", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16))),
    gridExtra::arrangeGrob(Supp_Figure_2D, top = grid::textGrob("D", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16)))),
  nrow = 2, widths = c(1, 1))

Supp_Figure_2_2 = gridExtra::arrangeGrob(
  grobs = list(
    gridExtra::arrangeGrob(Supp_Figure_2E, top = grid::textGrob("E", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16))),
    gridExtra::arrangeGrob(Supp_Figure_2F, top = grid::textGrob("F", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16))),
    gridExtra::arrangeGrob(Supp_Figure_2G, top = grid::textGrob("G", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16))),
    gridExtra::arrangeGrob(Supp_Figure_2H, top = grid::textGrob("H", x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16)))),
  nrow = 2, widths = c(1, 1))

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_2_1.png"), plot = Supp_Figure_2_1,
                width = 12, height = 9, units = "in", dpi = 300, bg = "white")

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_2_2.png"), plot = Supp_Figure_2_2,
                width = 12, height = 9, units = "in", dpi = 300, bg = "white")

###### Supp. Figure 3 : Continuous APM score vs. HR, all endpoints ######

# Assemble spline-based treatment HR plots across the continuous APM score for all four PANTHER endpoints.

panel_labels = c("A", "B", "C", "D")

hr_plots = lapply(endpoints, function(ep) {
  make_plot_hr_single_endpoint(data = all_clinical, biomarker = "APM_score", biomarker_label = "APM score",
                               endpoint = ep, subtype_label = "All",
                               adjust_vars = adjust_vars, adjust_terms = adjust_terms,
                               y_limits = c(-0.10, 4))
})

Supp_Figure_3 = gridExtra::arrangeGrob(
  grobs = Map(function(p, lab) {
    gridExtra::arrangeGrob(p, top = grid::textGrob(lab, x = 0, hjust = 0, gp = grid::gpar(fontface = "bold", fontsize = 16)))
  }, hr_plots, panel_labels),
  nrow = 2, widths = c(1, 1))

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_3.png"), plot = Supp_Figure_3,
                width = 15, height = 10, units = "in", dpi = 300, bg = "white")

###### Supp. Figure 4 : Continuous APM score within treatment arms, by ER status (A) and clinical subtype (B) ######

# Assemble within-arm score HRs by ER status and clinical subtype, with treatment-interaction and heterogeneity p-values.

# Adjust forest-plot row heights and lower padding before panel assembly.
tighten_bottom = function(g, pad_mm = 1, xaxis_mm = NULL, arrow_mm = NULL, last_row_mm = NULL) {
  rows_of = function(prefix) unique(g$layout$t[startsWith(g$layout$name, prefix)])
  if (!is.null(last_row_mm)) g$heights[min(rows_of("xaxis")) - 1] = unit(last_row_mm, "mm")  # heterogeneity row
  if (!is.null(xaxis_mm))    g$heights[rows_of("xaxis")] = unit(xaxis_mm, "mm")
  if (!is.null(arrow_mm))    g$heights[rows_of("arrow")] = unit(arrow_mm, "mm")
  g$heights[nrow(g)] = unit(pad_mm, "mm")                                                      # empty padding row
  g
}

Supp_Fig_4A = tighten_bottom(Supp_Fig_4A, pad_mm = 0, arrow_mm = 2)
Supp_Fig_4B = tighten_bottom(Supp_Fig_4B, pad_mm = 1, arrow_mm = 2)

# Same column widths in A and B; each panel at its natural height; page size follows the content
w_4 = unit.pmax(Supp_Fig_4A$widths, Supp_Fig_4B$widths)
Supp_Fig_4A$widths = w_4; Supp_Fig_4B$widths = w_4
panels_4 = Map(label_panel, list(Supp_Fig_4A, Supp_Fig_4B), c("A", "B"))
heights_4 = do.call(unit.c, lapply(panels_4, function(g) sum(g$heights)))
Supp_Figure_4 = arrangeGrob(grobs = panels_4, ncol = 1, heights = heights_4)

png(file.path(output_dir, "Supp_Figure_4.png"),
    width = convertWidth(sum(w_4), "in", valueOnly = TRUE) + 0.4,
    height = convertHeight(sum(heights_4), "in", valueOnly = TRUE),
    units = "in", res = 300, bg = "white")
grid.draw(Supp_Figure_4)
dev.off()

###### Supp. Figure 5 : Prognostication by APM score (PANTHER) ######

# Combine PANTHER APM-group KM plots for BCRFS, EFS, DDFS and OS, including risk tables.

Supp_Fig_5A$plot = Supp_Fig_5A$plot +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_5B$plot = Supp_Fig_5B$plot +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_5C$plot = Supp_Fig_5C$plot +
  ggplot2::labs(tag = "C") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_5D$plot = Supp_Fig_5D$plot +
  ggplot2::labs(tag = "D") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Figure_5 = survminer::arrange_ggsurvplots(
  x = list(Supp_Fig_5A, Supp_Fig_5C, Supp_Fig_5B, Supp_Fig_5D),
  nrow = 2,
  ncol = 2,
  print = FALSE,
  surv.plot.height = 0.80,
  risk.table.height = 0.20)

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_5.png"),
                plot = Supp_Figure_5, width = 3500, height = 3000, units = "px", dpi = 300, bg = "white")

###### Supp. Figure 6 : Prognostication by APM score (SCAN-B) ######

# Combine SCAN-B APM-group KM plots for BCFi, RFi, DRFi and OS, including risk tables.

Supp_Fig_6A$plot = Supp_Fig_6A$plot +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_6B$plot = Supp_Fig_6B$plot +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_6C$plot = Supp_Fig_6C$plot +
  ggplot2::labs(tag = "C") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_6D$plot = Supp_Fig_6D$plot +
  ggplot2::labs(tag = "D") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Figure_6 = survminer::arrange_ggsurvplots(
  x = list(Supp_Fig_6A, Supp_Fig_6C, Supp_Fig_6B, Supp_Fig_6D),
  nrow = 2,
  ncol = 2,
  print = FALSE,
  surv.plot.height = 0.80,
  risk.table.height = 0.20)

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_6.png"),
                plot = Supp_Figure_6, width = 3500, height = 3000, units = "px", dpi = 300, bg = "white")

###### Supp. Figure 7 : Forest Plots for EFS, DDFS, OS by APM group ######

# Combine adjusted treatment HRs overall and within median-defined APM groups for EFS, DDFS and OS.

Supp_Fig_7A_p = add_left_space(Supp_Fig_7A, 1)
Supp_Fig_7B_p = add_left_space(Supp_Fig_7B, 1)
Supp_Fig_7C_p = add_left_space(Supp_Fig_7C, 1)

png(file.path(output_dir, "Supp_Figure_7.png"), width = 9.5, height = 9, units = "in", res = 300, bg = "white")
gridExtra::grid.arrange(Supp_Fig_7A_p, Supp_Fig_7B_p, Supp_Fig_7C_p, ncol = 1, heights = c(1, 1, 1))
grid::grid.text("A", x = unit(0.03, "npc"), y = unit(0.98, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
grid::grid.text("B", x = unit(0.03, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
grid::grid.text("C", x = unit(0.03, "npc"), y = unit(0.32, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()

###### Supp. Figure 8 : KM curves for EFS, DDFS, OS (APM stratification per arm) ######

# Arrange treatment-arm KM curves within APM-high and APM-low groups for EFS, DDFS and OS, including risk tables.

Supp_Fig_8A$plot = Supp_Fig_8A$plot +
  ggplot2::labs(tag = "A") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_8B$plot = Supp_Fig_8B$plot +
  ggplot2::labs(tag = "B") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_8C$plot = Supp_Fig_8C$plot +
  ggplot2::labs(tag = "C") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_8D$plot = Supp_Fig_8D$plot +
  ggplot2::labs(tag = "D") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_8E$plot = Supp_Fig_8E$plot +
  ggplot2::labs(tag = "E") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Fig_8F$plot = Supp_Fig_8F$plot +
  ggplot2::labs(tag = "F") +
  ggplot2::theme(plot.tag = ggplot2::element_text(size = 18, face = "bold"),
                 plot.tag.position = c(0.01, 1.02),
                 plot.margin = ggplot2::margin(t = 20, r = 5, b = 5, l = 5))

Supp_Figure_8 = survminer::arrange_ggsurvplots(
  x = list(Supp_Fig_8A, Supp_Fig_8C, Supp_Fig_8E, Supp_Fig_8B, Supp_Fig_8D, Supp_Fig_8F),
  nrow = 3,
  ncol = 2,
  print = FALSE,
  surv.plot.height = 0.78,
  risk.table.height = 0.22)

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_8.png"),
                plot = Supp_Figure_8, width = 12, height = 14, units = "in", bg = "white")

###### Supp. Figure 9 : Forest Plots for EFS, DDFS, OS by arm ######

# Combine adjusted continuous-score HRs from separate treatment-arm models for EFS, DDFS and OS, with pooled interaction p-values.

Supp_Fig_9A_p = add_left_space(Supp_Fig_9A, 1)
Supp_Fig_9B_p = add_left_space(Supp_Fig_9B, 1)
Supp_Fig_9C_p = add_left_space(Supp_Fig_9C, 1)

png(file.path(output_dir, "Supp_Figure_9.png"), width = 11, height = 9, units = "in", res = 300, bg = "white")
gridExtra::grid.arrange(Supp_Fig_9A_p, Supp_Fig_9B_p, Supp_Fig_9C_p, ncol = 1, heights = c(1, 1, 1))
grid::grid.text("A", x = unit(0.03, "npc"), y = unit(0.98, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
grid::grid.text("B", x = unit(0.03, "npc"), y = unit(0.65, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
grid::grid.text("C", x = unit(0.03, "npc"), y = unit(0.32, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))
dev.off()

###### Supp. Figure 10 : APM score by Intrinsic subtypes ######

# Compare APM-group proportions across intrinsic subtypes in PANTHER and SCAN-B.

Supp_Figure_10 = Supp_Figure_10A / Supp_Figure_10B +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(filename = file.path(output_dir, "Supp_Figure_10.png"), plot = Supp_Figure_10, width = 8, height = 10, units = "in", dpi = 300, bg = "white")

###### Supp. Figure 11 : Deconvolution boxplots (MCPcounter) ######

# Combine MCP-counter score comparisons between APM groups in PANTHER and SCAN-B with a shared legend.

Supp_Figure_11A = Supp_Figure_11A + ggplot2::labs(y = "MCP-counter score")
Supp_Figure_11B = Supp_Figure_11B + ggplot2::labs(y = "MCP-counter score")

Supp_Figure_11 = Supp_Figure_11A / Supp_Figure_11B +
  patchwork::plot_layout(
    guides = "collect",
    axes = "collect",
    axis_titles = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.tag = ggplot2::element_text(size = 14, face = "bold", margin = ggplot2::margin(r = 10)))

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_11.png"), plot = Supp_Figure_11, width = 9, height = 11, units = "in", dpi = 300, bg = "white")

###### Supp. Figure 12 : Deconvolution boxplots (Danaher) ######

# Combine Danaher marker-score comparisons between APM groups in PANTHER and SCAN-B with a shared legend.

Supp_Figure_12A = Supp_Figure_12A + ggplot2::labs(y = "Danaher immune cell score")
Supp_Figure_12B = Supp_Figure_12B + ggplot2::labs(y = "Danaher immune cell score")

Supp_Figure_12 = Supp_Figure_12A / Supp_Figure_12B +
  patchwork::plot_layout(
    guides = "collect",
    axes = "collect",
    axis_titles = "collect") +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.tag = ggplot2::element_text(size = 14, face = "bold", margin = ggplot2::margin(r = 10)))

ggplot2::ggsave(filename = file.path(output_dir, "Supp_Figure_12.png"), plot = Supp_Figure_12, width = 13, height = 13, units = "in", dpi = 300, bg = "white")

###### Supp. Figure 13 : QC Plots ######

# Placeholder for the RNA-seq QC summary; the figure is made by workflows/rnaseq_qc_and_subtyping/scripts/04_supplementary_figure_S13.R.

###### Supp. Figure 14 : Workflow for the APM score development ######

# Placeholder for the APM score development; The diagram was created in Power Point.

###### ----- Supplementary Tables -----
###### Supp. Table 5: Immune-cell scores and differential treatment benefit (BCRFS, PANTHER) ######

# Format and export immune-score interaction results and joint cytotoxic/humoral comparisons from the earlier exploratory analyses.

imm_order = c("APM score", paste("MCPcounter:", mcpcounter_order), paste("Danaher:", danaher_order), "14-gene IgG")
supp_imm = merge(q1, q2, by = "score", all.x = TRUE)
stopifnot(!anyNA(match(supp_imm$score, imm_order)))
supp_imm = supp_imm[order(match(supp_imm$score, imm_order)), ]

# Format p-values and display unavailable values as an en dash.
fmt_or_dash = function(p) ifelse(is.na(p), "\u2013", format_p_regular(p))

## Part A
# Summarize score correlations, within-arm HRs, individual interaction tests, BH-adjusted p-values and conditional APM interaction tests.
supp_imm_tab = data.frame(
  Method = ifelse(supp_imm$score == "APM score", "Reference",
                  ifelse(grepl("^MCPcounter", supp_imm$score), "MCP-counter",
                         ifelse(grepl("^Danaher", supp_imm$score), "Danaher", "Signature"))),
  Score = sub("^(MCPcounter|Danaher): ", "", supp_imm$score),
  rho = ifelse(supp_imm$score == "APM score", "\u2013", formatC(supp_imm$rho_with_APM, format = "f", digits = 2)),
  hr_fec = supp_imm$HR_per_SD_FEC_D,
  hr_tdd = supp_imm$HR_per_SD_tddEC_D,
  p_int = format_p_regular(supp_imm$p_interaction),
  p_bh = fmt_or_dash(supp_imm$p_interaction_BH),
  p_apm = fmt_or_dash(supp_imm$p_APM_interaction_given_immune))

colnames(supp_imm_tab) = c("Method", "Score", "Spearman \u03c1 with APM score", "HR per SD, FEC/D", "HR per SD, tDD EC/D",
                           "p interaction", "p interaction, BH", "APM p interaction, given score interaction")

## Part B
# Summarize conditional interaction tests from joint cytotoxic/CD8 and B-cell/IgG models.
supp_imm_pairs = data.frame(
  cyto = sub("^MCPcounter", "MCP-counter", q3$cytotoxic),
  humo = sub("^MCPcounter", "MCP-counter", q3$humoral),
  rho = formatC(q3$rho, format = "f", digits = 2),
  p_cyto = format_p_regular(q3$p_cytotoxic_interaction_given_humoral),
  p_humo = format_p_regular(q3$p_humoral_interaction_given_cytotoxic))

colnames(supp_imm_pairs) = c("Cytotoxic / CD8\u207a score", "B-cell / IgG score", "Spearman \u03c1",
                             "p cytotoxic interaction, given B-cell/IgG", "p B-cell/IgG interaction, given cytotoxic")

write_table(supp_imm_tab, file.path(output_dir, "Supp_Table_5_immune_interactions_A.csv"))
write_table(supp_imm_pairs, file.path(output_dir, "Supp_Table_5_immune_interactions_B.csv"))

gt::gt(supp_imm_tab, groupname_col = "Method") |>
  gt::tab_header(title = "A. Immune-cell scores and differential benefit from tDD EC/D vs FEC/D (BCRFS)") |>
  gt::cols_align(align = "center", columns = -Score) |>
  gt::cols_width(Score ~ gt::px(190), c(`HR per SD, FEC/D`, `HR per SD, tDD EC/D`) ~ gt::px(165), everything() ~ gt::px(120)) |>
  gt::tab_style(gt::cell_text(whitespace = "nowrap"), gt::cells_body(columns = c(`HR per SD, FEC/D`, `HR per SD, tDD EC/D`))) |>
  gt::tab_source_note(paste0("HRs per SD from separate adjusted Cox models in each arm; p interaction from adjusted likelihood-ratio tests (n = ", nrow(d_imm), "). BH: Benjamini\u2013Hochberg across the ", length(immune_vars), " immune scores. Last column: APM score interaction in a joint model that also includes the immune score's own treatment interaction.")) |>
  gt::gtsave(file.path(output_dir, "Supp_Table_5_immune_interactions_A.png"), vwidth = 1200, expand = 20)

gt::gt(supp_imm_pairs) |>
  gt::tab_header(title = "B. Cytotoxic vs humoral immune scores in joint interaction models (BCRFS)") |>
  gt::cols_align(align = "center", columns = -c(1, 2)) |>
  gt::cols_width(1:2 ~ gt::px(270), 3 ~ gt::px(100), everything() ~ gt::px(210)) |>
  gt::tab_source_note(paste0("Each row is one adjusted Cox model with both treatment interactions; p-values are nominal likelihood-ratio tests for each interaction given the other (n = ", nrow(d_imm), ").")) |>
  gt::gtsave(file.path(output_dir, "Supp_Table_5_immune_interactions_B.png"), vwidth = 1300, expand = 20)

###### Session Info ######

# Save the R version, platform and loaded-package versions for reproducibility.

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
