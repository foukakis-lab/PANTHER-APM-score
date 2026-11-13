## =============================================================================
## PANTHER: APM discovery pipeline with biologically restricted partners,
## and its whole-pipeline permutation test
##
## The pipeline (BCRFS; continuous TREAT x score interaction, 1-df LRT, adjusted
## model, the same complete cases throughout):
##
##   Step 1  Gate. Score the five multi-gene BC360 APM subgroups with GSVA and test each for a continuous TREAT x score interaction. 
##           Keep the subgroups with p < GATE_P (0.05). If none passes, the pipeline selects nothing.
##   Step 2  Class. Each kept subgroup is MHC-I or MHC-II (its BC360 main group).
##           The APM80 main groups that act directly on that class (see Methods):
##             MHC-I  <- Regulators, Repressors, Post-transcriptional MHC-I regulation
##             MHC-II <- Regulators, Repressors
##           Partner candidates = the subgroups of those main groups with >= 2 genes.
##           (All three post-transcriptional subgroups are single genes, so in practice MHC-I and MHC-II anchors get the same 7 candidates:
##           RFX, IFN5, IFNGR/JAK2, IRF1/IRF2, PRC2, DNMT, HDAC.)
##   Step 3  Combination. For every kept subgroup (anchor) and every candidate, in both directions:  
##           APM score = GSVA(anchor) + s * GSVA(partner),  s = +1 / -1
##           (the same signed GSVA as the main script: an unscaled difference or sum, then standardised; HR per SD). 
##           The combination with the smallest interaction p (largest chi2) over all anchors and candidates is selected.
##           Partner gate (PARTNER_GATE = TRUE): the anchors alone also compete, 
##           so a partner is kept only if its combination has a smaller interaction p than every passing anchor alone (ties: no partner). 
##           Otherwise the signature is the best anchor alone.
##   Statistic  T = chi2 of the selected signature; T = 0 when no anchor passes.
##
## Permutation: treatment labels are shuffled (optionally within strata) and steps
## 1-3 are re-run on each permuted data set. GSVA scores do not depend on the labels, so they are computed once.
##   whole-pipeline p = (1 + #{T_null >= T_obs}) / (1 + N_PERM)
## Null runs in which no anchor passes the gate count as T = 0 (no discovery).
##
## The p-value is conditional on choices made before this pipeline was written down, some after seeing outcomes.
## The permutation covers only the steps above. The partner gate was added after the first runs; 
## it gives the null runs one more option per anchor, so it can only raise the whole-pipeline p.
##
## Whatever the pipeline selects on the observed data is the result.
##
## Requires (in the session, from the main script): tmm_normalized_tpm_symbols (gene symbol x sample,
## log2(TPM+1)); all_clinical (sample_id, TREAT, timeBCRFS, BCRFS.status and the adjustment variables);
## output_dir. Optional cross-checks: FIXED_APM80_SUBGROUPS, BC360_APM_SUBGROUPS.
## Packages: survival, GSVA, ggplot2, future, future.apply.

## =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(future)
  library(future.apply)
})

## -----------------------------------------------------------------------------
## 0. Settings (fix before looking at the output)
## -----------------------------------------------------------------------------

if (!exists('OUT_DIR_OVERRIDE'))
  OUT_DIR = file.path(output_dir, 'APM_discovery_permutation') else
    OUT_DIR = OUT_DIR_OVERRIDE

ENDPOINT     = 'BCRFS'
TREAT_LEVELS = c('FEC-D', 'tddEC-D')      # first = reference
GATE_P       = 0.05                       # step 1: an anchor needs interaction p < GATE_P
ANCHOR_RULE  = 'all_passing'              # 'all_passing' : every BC360 subgroup with p < GATE_P is an anchor (pre-specified)
                                          # 'best_passing': only the passing subgroup with the lowest p
PARTNER_RULE = 'direct_actors'            # 'direct_actors': step 2 as above (pre-specified)
                                          # 'all_multi'    : sensitivity; every multi-gene APM80 subgroup sharing no gene with the anchor
PARTNER_GATE = TRUE                       # keep a partner only if it lowers the interaction p of the anchors alone (FALSE: always one partner)
MIN_GENES    = 2                          # candidates need >= MIN_GENES genes
SIGNS        = c(1, -1)                   # both biological directions
GSVA_BACKGROUND = 'cohort'                # all non-constant genes of the PANTHER matrix, as in the main script
GSVA_KCDF    = 'Gaussian'
N_PERM       = 5000                       # Monte Carlo SE at p ~ 0.1: about 0.004 (1000 permutations: about 0.009)
SEED         = 2027
N_WORKERS    = max(1, parallel::detectCores() - 1)
CHUNK        = 50                         # permutations per saved chunk (resume-safe)
PERM_STRATA  = 'ERpos'                    # shuffle arms within ER strata, matching the ER-stratified Cox model (strata(ERpos)); NULL = over all patients
OTHER_ENDPOINTS = c('EFS', 'DDFS', 'OS')  # descriptive only, for the selected combination

DIRECT_ACTORS = list(MHC_I  = c('Regulators', 'Repressors', 'Post-transcriptional MHC-I regulation'),
                     MHC_II = c('Regulators', 'Repressors'))

ANCHOR_RULE  = match.arg(ANCHOR_RULE, c('all_passing', 'best_passing'))
PARTNER_RULE = match.arg(PARTNER_RULE, c('direct_actors', 'all_multi'))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists('adjust_terms')) {
  adjust_terms = 'strata(ERpos) + newHER2 + scale(tumdiam) + scale(posnodes) + histgrade'
  message('adjust_terms not in the session; using: ', adjust_terms)
}

## -----------------------------------------------------------------------------
## 1. Gene-set definitions (frozen membership, as in the two panel figures)
## -----------------------------------------------------------------------------
BC360_ANCHORS = list(
  BC_MHCI_surface = list(label = 'Surface MHC-I chains (BC360)', class = 'MHC_I',
                         genes = c('HLA-A', 'HLA-B', 'HLA-C')),
  
  BC_MHCI_TAP     = list(label = 'Peptide transport and loading (BC360)', class = 'MHC_I',
                         genes = c('TAP1', 'TAP2', 'TAPBP')),
  
  BC_MHCI_iprot   = list(label = 'Immunoproteasome (BC360)', class = 'MHC_I',
                         genes = c('PSMB9', 'PSMB10')),
  
  BC_MHCII_chains = list(label = 'HLA6: HLA-DR/-DQ/-DP chains (BC360)', class = 'MHC_II',
                         genes = c('HLA-DRA', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DPA1', 'HLA-DPB1')),
  
  BC_MHCII_DMDO   = list(label = 'HLA-DM/HLA-DO (BC360)', class = 'MHC_II',
                         genes = c('HLA-DMA', 'HLA-DMB', 'HLA-DOB')))

sg = function(main, label, genes) list(main = main, label = label, genes = genes)
APM80 = list(
  MHCI_1 = sg('MHC-I', 'Surface MHC-I chains and B2M', c('HLA-A', 'HLA-B', 'HLA-C', 'B2M')),
  MHCI_2 = sg('MHC-I', 'Peptide transport and loading', c('TAP1', 'TAP2', 'TAPBP', 'TAPBPL')),
  MHCI_3 = sg('MHC-I', 'Immunoproteasome', c('PSMB8', 'PSMB9', 'PSMB10')),
  MHCI_4 = sg('MHC-I', 'Constitutive proteasome', c('PSMB5', 'PSMB6', 'PSMB7')),
  MHCI_5 = sg('MHC-I', 'Proteasome activators', c('PSME1', 'PSME2')),
  MHCI_6 = sg('MHC-I', 'ER peptide trimming', c('ERAP1', 'ERAP2')),
  MHCI_7 = sg('MHC-I', 'ER folding and loading', c('CALR', 'CANX', 'PDIA3')),
  MHCII_1 = sg('MHC-II', 'HLA9: HLA-DR/-DQ/-DP chains', c('HLA-DRA', 'HLA-DRB1', 'HLA-DRB5', 'HLA-DQA1', 'HLA-DQA2',
                                                          'HLA-DQB1', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPB1')),
  MHCII_2 = sg('MHC-II', 'HLA-DM/HLA-DO', c('HLA-DMA', 'HLA-DMB', 'HLA-DOA', 'HLA-DOB')),
  MHCII_3 = sg('MHC-II', 'Invariant chain', 'CD74'),
  MHCII_4 = sg('MHC-II', 'Endolysosomal proteases', c('CTSS', 'CTSL', 'CTSB', 'CTSD', 'LGMN')),
  MHCII_5 = sg('MHC-II', 'Lysosomal thiol reductase', 'IFI30'),
  POST_1 = sg('Post-transcriptional MHC-I regulation', 'NBR1', 'NBR1'),
  POST_2 = sg('Post-transcriptional MHC-I regulation', 'MAL2', 'MAL2'),
  POST_3 = sg('Post-transcriptional MHC-I regulation', 'SPPL3', 'SPPL3'),
  XP_1 = sg('Cross presentation', 'Trafficking and redox', c('LNPEP', 'RAB7A', 'SEC22B', 'CYBB')),
  XP_2 = sg('Cross presentation', 'Antigen export / dislocation', c('DERL1', 'DERL2', 'VCP', 'SEC61A1')),
  REG_1 = sg('Regulators', 'NLRC5', 'NLRC5'),
  REG_2 = sg('Regulators', 'CIITA', 'CIITA'),
  REG_3 = sg('Regulators', 'RFX', c('RFX5', 'RFXANK', 'RFXAP')),
  REG_4 = sg('Regulators', 'IFNG', 'IFNG'),
  REG_5 = sg('Regulators', 'IFN5', c('JAK1', 'TYK2', 'STAT1', 'STAT2', 'IRF9')),
  REG_6 = sg('Regulators', 'IFNGR/JAK2', c('IFNGR1', 'IFNGR2', 'JAK2')),
  REG_7 = sg('Regulators', 'IRF1/IRF2', c('IRF1', 'IRF2')),
  NC_1 = sg('Nonclassical evasion', 'Nonclassical HLA-I', c('HLA-E', 'HLA-G', 'HLA-F')),
  REP_1 = sg('Repressors', 'PRC2', c('EZH2', 'SUZ12', 'EED')),
  REP_2 = sg('Repressors', 'DNMT', c('DNMT1', 'DNMT3A', 'DNMT3B')),
  REP_3 = sg('Repressors', 'HDAC', c('HDAC1', 'HDAC2', 'HDAC3')))
stopifnot(length(APM80) == 28L, length(unique(unlist(lapply(APM80, `[[`, 'genes')))) == 80L)

## cross-checks against the definitions in the session, if present
if (exists('FIXED_APM80_SUBGROUPS')) {
  ses = FIXED_APM80_SUBGROUPS
  hit = vapply(APM80, function(s) any(vapply(ses, setequal, logical(1), s$genes)), logical(1))
  if (!all(hit) || length(ses) != 28L)
    stop('APM80 subgroups differ from FIXED_APM80_SUBGROUPS: ', paste(names(APM80)[!hit], collapse = ', '))
  message('APM80 definitions match FIXED_APM80_SUBGROUPS (28 subgroups).')
}

if (exists('BC360_APM_SUBGROUPS')) {
  hit = vapply(BC360_ANCHORS, function(a) any(vapply(BC360_APM_SUBGROUPS, setequal, logical(1), a$genes)), logical(1))
  if (!all(hit)) stop('BC360 anchors not found in BC360_APM_SUBGROUPS: ', paste(names(BC360_ANCHORS)[!hit], collapse = ', '))
  multi = Filter(function(g) length(g) >= 2, BC360_APM_SUBGROUPS)
  if (length(multi) != length(BC360_ANCHORS)) stop('BC360_APM_SUBGROUPS has ', length(multi), ' multi-gene subgroups, not 5.')
  message('BC360 anchors match BC360_APM_SUBGROUPS (5 multi-gene subgroups).')
}

## step 2: partner candidates per anchor
multi_ids = names(APM80)[vapply(APM80, function(s) length(s$genes) >= MIN_GENES, logical(1))]
candidates_for = function(a) {
  A = BC360_ANCHORS[[a]]
  if (PARTNER_RULE == 'direct_actors') {
    groups = DIRECT_ACTORS[[A$class]]
    ids = names(APM80)[vapply(APM80, function(s) s$main %in% groups, logical(1))]
  } else ids = names(APM80)
  ids = intersect(ids, multi_ids)
  ids[vapply(ids, function(k) !length(intersect(APM80[[k]]$genes, A$genes)), logical(1))]
}
CANDS = lapply(setNames(names(BC360_ANCHORS), names(BC360_ANCHORS)), candidates_for)
cand_tab = do.call(rbind, lapply(names(BC360_ANCHORS), function(a) {
  A = BC360_ANCHORS[[a]]
  pool = if (PARTNER_RULE == 'direct_actors') names(APM80)[vapply(APM80, function(s) s$main %in% DIRECT_ACTORS[[A$class]], logical(1))] else names(APM80)
  data.frame(anchor = a, anchor_class = A$class, apm80_subgroup = pool, main_group = vapply(APM80[pool], `[[`, '', 'main'),
             label = vapply(APM80[pool], `[[`, '', 'label'), n_genes = lengths(lapply(APM80[pool], `[[`, 'genes')),
             candidate = pool %in% CANDS[[a]],
             reason = ifelse(pool %in% CANDS[[a]], 'candidate',
                             ifelse(lengths(lapply(APM80[pool], `[[`, 'genes')) < MIN_GENES, 'single gene', 'shares genes with anchor')),
             stringsAsFactors = FALSE, row.names = NULL)
}))
lab80 = vapply(APM80, `[[`, '', 'label')
labBC = vapply(BC360_ANCHORS, `[[`, '', 'label')
shortBC = c(BC_MHCI_surface = 'MHC-I chains', BC_MHCI_TAP = 'TAP', BC_MHCI_iprot = 'PSMB9/10',
            BC_MHCII_chains = 'HLA6', BC_MHCII_DMDO = 'DM/DO')
combo_label = function(a, p, s) paste(shortBC[a], ifelse(s > 0, '+', '-'), lab80[p])
sig_label = function(a, p, s) if (is.na(p) || !nzchar(p)) unname(shortBC[a]) else unname(combo_label(a, p, s))   # anchor alone or combination
message('Partner rule: ', PARTNER_RULE, '; candidates per anchor class:')
for (cl in names(DIRECT_ACTORS)) {
  a1 = names(BC360_ANCHORS)[vapply(BC360_ANCHORS, `[[`, '', 'class') == cl][1]
  message('  ', cl, ': ', paste(lab80[CANDS[[a1]]], collapse = ', '))
}

## -----------------------------------------------------------------------------
## 2. GSVA scores (expression only; computed once)
## -----------------------------------------------------------------------------
for (obj in c('tmm_normalized_tpm_symbols', 'all_clinical'))
  if (!exists(obj)) stop(obj, ' not found in the session.')
expr = as.matrix(tmm_normalized_tpm_symbols)
ids = intersect(as.character(all_clinical$sample_id), colnames(expr))
message('PANTHER: ', length(ids), ' RNA-matched patients')
X = expr[, ids, drop = FALSE]
if (max(X, na.rm = TRUE) > 50) warning('Expression max > 50: the matrix may not be on the log2 scale.')
X = X[rowSums(!is.finite(X)) == 0, , drop = FALSE]
X = X[apply(X, 1, stats::sd) > 1e-8, , drop = FALSE]

SETS = c(lapply(BC360_ANCHORS, `[[`, 'genes'), lapply(APM80[multi_ids], `[[`, 'genes'))
present = lapply(SETS, intersect, rownames(X))
cov_tab = data.frame(set = names(SETS), n_defined = lengths(SETS), n_present = lengths(present),
                     missing = vapply(names(SETS), function(k) paste(setdiff(SETS[[k]], present[[k]]), collapse = ' '), ''),
                     row.names = NULL)
if (any(cov_tab$n_present < cov_tab$n_defined))
  warning('Genes missing or constant in the GSVA background:\n', paste(capture.output(print(cov_tab[cov_tab$n_present < cov_tab$n_defined, ])), collapse = '\n'))
if (any(cov_tab$n_present < 2)) stop('A set has fewer than 2 scorable genes; see the coverage table.')
message(sprintf('GSVA (%s background, %d genes x %d patients, kcdf %s) ...', GSVA_BACKGROUND, nrow(X), ncol(X), GSVA_KCDF))
t0 = Sys.time()
es = if (exists('gsvaParam', envir = asNamespace('GSVA'), inherits = FALSE)) {
  GSVA::gsva(GSVA::gsvaParam(X, present, kcdf = GSVA_KCDF, minSize = 1, maxSize = Inf, maxDiff = TRUE), verbose = FALSE)
} else GSVA::gsva(X, present, method = 'gsva', kcdf = GSVA_KCDF, min.sz = 1, max.sz = Inf, mx.diff = TRUE, verbose = FALSE)
S = t(as.matrix(es))[ids, names(SETS), drop = FALSE]
message(sprintf('  done in %.0f s', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

## -----------------------------------------------------------------------------
## 3. Analysis set (complete cases for the adjusted model)
## -----------------------------------------------------------------------------
cl = as.data.frame(all_clinical)
cl$sample_id = as.character(cl$sample_id)
cl = cl[match(ids, cl$sample_id), ]
time_col = paste0('time', ENDPOINT); status_col = paste0(ENDPOINT, '.status')
adj_vars = all.vars(as.formula(paste('~', adjust_terms)))
need = unique(c('sample_id', time_col, status_col, 'TREAT', adj_vars))
miss = setdiff(need, names(cl)); if (length(miss)) stop('Not in all_clinical: ', paste(miss, collapse = ', '))
keep = complete.cases(cl[, need])
D = cl[keep, need]
D$TREAT = factor(as.character(D$TREAT), levels = TREAT_LEVELS)
if (anyNA(D$TREAT)) stop('TREAT has levels other than ', paste(TREAT_LEVELS, collapse = ', '))
D$.time = D[[time_col]]; D$.status = D[[status_col]]
if (!all(D$.status %in% c(0, 1))) stop(status_col, ' must be 0/1')
rows = match(D$sample_id, ids)     # rows of S for the analysis patients
message(sprintf('%s: %d patients, %d events (%s %d / %s %d)', ENDPOINT, nrow(D), sum(D$.status),
                TREAT_LEVELS[1], sum(D$.status[D$TREAT == TREAT_LEVELS[1]]),
                TREAT_LEVELS[2], sum(D$.status[D$TREAT == TREAT_LEVELS[2]])))

F0 = as.formula(paste('Surv(.time, .status) ~ TREAT + SCORE +', adjust_terms))
F1 = as.formula(paste('Surv(.time, .status) ~ TREAT * SCORE +', adjust_terms))
INT = paste0('TREAT', TREAT_LEVELS[2], ':SCORE')

## a score over all RNA-matched patients, standardised over them (as in the main script)
zs = function(x) (x - mean(x)) / stats::sd(x)
anchor_z = sapply(names(BC360_ANCHORS), function(a) zs(S[, a])[rows])
combo_grid = do.call(rbind, lapply(names(BC360_ANCHORS), function(a)
  expand.grid(anchor = a, partner = CANDS[[a]], sign = SIGNS, stringsAsFactors = FALSE)))
combo_z = sapply(seq_len(nrow(combo_grid)), function(i)
  zs(S[, combo_grid$anchor[i]] + combo_grid$sign[i] * S[, combo_grid$partner[i]])[rows])
message(nrow(combo_grid), ' possible combinations (anchors x candidates x 2 directions).')

## -----------------------------------------------------------------------------
## 4. The pipeline
## -----------------------------------------------------------------------------
int_test = function(d, z) {
  d$SCORE = z
  tryCatch({
    f0 = suppressWarnings(coxph(F0, data = d))
    f1 = suppressWarnings(coxph(F1, data = d))
    c(chi2 = max(0, 2 * (f1$loglik[2] - f0$loglik[2])), beta = unname(coef(f1)[INT]), se = sqrt(vcov(f1)[INT, INT]))
  }, error = function(e) c(chi2 = NA_real_, beta = NA_real_, se = NA_real_))
}

run_pipeline = function(d, trace = FALSE) {
  ## step 1: gate
  st1 = t(vapply(colnames(anchor_z), function(a) int_test(d, anchor_z[, a]), numeric(3)))
  p1 = pchisq(st1[, 'chi2'], 1, lower.tail = FALSE)
  pass = rownames(st1)[is.finite(p1) & p1 < GATE_P]
  if (ANCHOR_RULE == 'best_passing' && length(pass) > 1) pass = pass[which.max(st1[pass, 'chi2'])]
  out = list(step1 = cbind(st1, p = p1), passed = pass, T = 0, sel = NA_integer_, sel_anchor = NA_character_,
             T_combo = 0, T_anchor = 0, grid = NULL)
  if (!length(pass)) return(out)
  ## steps 2-3: candidates of the passing anchors, both directions
  idx = which(combo_grid$anchor %in% pass)
  res = t(vapply(idx, function(i) int_test(d, combo_z[, i]), numeric(3)))
  chi = res[, 'chi2']; chi[!is.finite(chi)] = -Inf
  j = which.max(chi)
  if (is.finite(chi[j])) { out$T = chi[j]; out$sel = idx[j]; out$sel_anchor = combo_grid$anchor[idx[j]]; out$T_combo = chi[j] }
  ## partner gate: the best anchor alone wins unless a combination has a larger chi2 (ties: no partner)
  a_best = pass[which.max(st1[pass, 'chi2'])]
  out$T_anchor = st1[a_best, 'chi2']
  if (PARTNER_GATE && out$T_anchor >= out$T) { out$T = out$T_anchor; out$sel = NA_integer_; out$sel_anchor = a_best }
  if (trace) out$grid = cbind(combo_grid[idx, ], res, p = pchisq(res[, 'chi2'], 1, lower.tail = FALSE))
  out
}

## -----------------------------------------------------------------------------
## 5. Observed run
## -----------------------------------------------------------------------------
t0 = Sys.time()
obs = run_pipeline(D, trace = TRUE)
sec_obs = as.numeric(difftime(Sys.time(), t0, units = 'secs'))
T_obs = obs$T

s1 = data.frame(anchor = rownames(obs$step1), label = labBC[rownames(obs$step1)],
                class = vapply(BC360_ANCHORS[rownames(obs$step1)], `[[`, '', 'class'),
                genes = vapply(BC360_ANCHORS[rownames(obs$step1)], function(a) paste(a$genes, collapse = ';'), ''),
                chi2 = obs$step1[, 'chi2'], p = obs$step1[, 'p'], int_HR_per_SD = exp(obs$step1[, 'beta']),
                passes_gate = rownames(obs$step1) %in% obs$passed, row.names = NULL)
write.csv(s1, file.path(OUT_DIR, '01_step1_BC360_gate.csv'), row.names = FALSE)
write.csv(cand_tab, file.path(OUT_DIR, '02_step2_partner_candidates.csv'), row.names = FALSE)
write.csv(cov_tab, file.path(OUT_DIR, '00b_gene_coverage.csv'), row.names = FALSE)

if (!length(obs$passed)) {
  message('\nObserved: no BC360 subgroup passes p < ', GATE_P, '. The pipeline selects nothing (T = 0).')
  sel_label = 'none'
} else {
  g = obs$grid
  g$combination = combo_label(g$anchor, g$partner, g$sign)
  g$direction = ifelse(g$sign > 0, 'plus', 'minus')
  g$partner_label = lab80[g$partner]; g$partner_main = vapply(APM80[g$partner], `[[`, '', 'main')
  g$int_HR_per_SD = exp(g$beta)
  i = obs$sel                           # NA when the partner gate keeps the anchor alone
  a_sel = obs$sel_anchor
  partner_kept = !is.na(i)
  best_combo = g$combination[which.max(ifelse(is.finite(g$chi2), g$chi2, -Inf))]
  g$selected = partner_kept & seq_len(nrow(g)) == which.max(ifelse(is.finite(g$chi2), g$chi2, -Inf))
  g = g[order(-g$chi2), c('combination', 'anchor', 'partner', 'partner_label', 'partner_main', 'direction',
                          'chi2', 'p', 'int_HR_per_SD', 'selected')]
  write.csv(g, file.path(OUT_DIR, '03_step3_combinations.csv'), row.names = FALSE)
  sel_label = if (partner_kept) sig_label(a_sel, combo_grid$partner[i], combo_grid$sign[i]) else sig_label(a_sel, NA, NA)

  ## selected signature: effect sizes
  D$SIG = if (partner_kept) combo_z[, i] else anchor_z[, a_sel]
  f1 = coxph(as.formula(paste('Surv(.time, .status) ~ TREAT * SIG +', adjust_terms)), data = D)
  iname = paste0('TREAT', TREAT_LEVELS[2], ':SIG')
  b = coef(f1); V = vcov(f1)
  hr = function(w) { est = sum(w * b[names(w)]); se = sqrt(drop(t(w) %*% V[names(w), names(w)] %*% w)); exp(est + c(0, -1.96, 1.96) * se) }
  eff = rbind(interaction = hr(setNames(1, iname)), slope_ref = hr(setNames(1, 'SIG')), slope_other = hr(setNames(c(1, 1), c('SIG', iname))))
  a_row = s1[s1$anchor == a_sel, ]
  sig_tab = data.frame(selected = sel_label, anchor = a_sel, partner = if (partner_kept) combo_grid$partner[i] else '',
                       direction = if (partner_kept) ifelse(combo_grid$sign[i] > 0, 'plus', 'minus') else 'none',
                       partner_genes = if (partner_kept) paste(APM80[[combo_grid$partner[i]]]$genes, collapse = ';') else '',
                       partner_gate = PARTNER_GATE, partner_kept = partner_kept,
                       best_combination = best_combo, best_combination_chi2 = obs$T_combo, best_anchor_alone_chi2 = obs$T_anchor,
                       n = nrow(D), events = sum(D$.status), chi2 = T_obs, p_nominal = pchisq(T_obs, 1, lower.tail = FALSE),
                       int_HR_per_SD = eff[1, 1], int_low95 = eff[1, 2], int_high95 = eff[1, 3],
                       score_HR_per_SD_ref_arm = eff[2, 1], score_HR_per_SD_other_arm = eff[3, 1],
                       anchor_alone_chi2 = a_row$chi2, anchor_alone_p = a_row$p, stringsAsFactors = FALSE)
  names(sig_tab) = sub('ref_arm', TREAT_LEVELS[1], sub('other_arm', TREAT_LEVELS[2], names(sig_tab)))
  write.csv(sig_tab, file.path(OUT_DIR, '04_selected_combination.csv'), row.names = FALSE)

  ## other endpoints, descriptive (nominal; they share events with BCRFS)
  other = do.call(rbind, lapply(OTHER_ENDPOINTS, function(ep) {
    tc = paste0('time', ep); sc = paste0(ep, '.status')
    if (!all(c(tc, sc) %in% names(cl))) return(NULL)
    e = cbind(D[, c('sample_id', 'TREAT', 'SIG', adj_vars)], cl[match(D$sample_id, cl$sample_id), c(tc, sc)])
    e = e[complete.cases(e), ]; e$.time = e[[tc]]; e$.status = e[[sc]]
    r = int_test(e, e$SIG)
    data.frame(endpoint = ep, n = nrow(e), events = sum(e$.status), chi2 = r[['chi2']],
               p_nominal = pchisq(r[['chi2']], 1, lower.tail = FALSE), int_HR_per_SD = exp(r[['beta']]))
  }))
  if (!is.null(other)) write.csv(other, file.path(OUT_DIR, '04b_selected_other_endpoints_descriptive.csv'), row.names = FALSE)

  message(sprintf('\nObserved run (%.1f s): passing anchor(s): %s\n  best combination %s (chi2 %.2f) vs best anchor alone (chi2 %.2f): %s\n  selected: %s\n  chi2 %.2f (nominal p %.2g), interaction HR per SD %.2f (%.2f-%.2f)',
                  sec_obs, paste(shortBC[obs$passed], collapse = ', '), best_combo, obs$T_combo, obs$T_anchor,
                  if (!PARTNER_GATE) 'partner gate off' else if (partner_kept) 'partner kept' else 'no partner',
                  sel_label, T_obs, pchisq(T_obs, 1, lower.tail = FALSE), eff[1, 1], eff[1, 2], eff[1, 3]))
}

## -----------------------------------------------------------------------------
## 6. Whole-pipeline permutation
## -----------------------------------------------------------------------------
perm_file = file.path(OUT_DIR, '05_permutation_replicates.csv')
summary_file = file.path(OUT_DIR, '06_permutation_summary.csv')
SETTINGS_KEY = paste(SEED, if (is.null(PERM_STRATA)) 'none' else PERM_STRATA, GATE_P, ANCHOR_RULE, PARTNER_RULE,
                     GSVA_BACKGROUND, MIN_GENES, paste(SIGNS, collapse = '/'), nrow(D),
                     if (PARTNER_GATE) 'partnerGate' else 'onePartner', sep = '|')

if (N_PERM > 0) {
  set.seed(SEED)    # all permutations drawn up front, so a resumed run uses the same ones
  strata = if (is.null(PERM_STRATA)) rep(1, nrow(D)) else D[[PERM_STRATA]]
  if (anyNA(strata)) stop('PERM_STRATA has missing values')
  PERMS = lapply(seq_len(N_PERM), function(b) {
    idx = seq_len(nrow(D))
    for (s in unique(strata)) { w = which(strata == s); idx[w] = w[sample.int(length(w))] }
    idx
  })
  done = if (file.exists(perm_file)) read.csv(perm_file, stringsAsFactors = FALSE) else NULL
  if (!is.null(done)) {
    if (!all(done$settings == SETTINGS_KEY)) stop(perm_file, ' was written with other settings; move it away to start again.')
    message(nrow(done), ' permutations already done; resuming.')
  }
  todo = setdiff(seq_len(N_PERM), done$perm)

  one_perm = function(b) {
    dp = D; dp$TREAT = D$TREAT[PERMS[[b]]]
    r = run_pipeline(dp)
    i = r$sel; a = r$sel_anchor
    data.frame(perm = b, settings = SETTINGS_KEY, n_passing = length(r$passed),
               passing = if (length(r$passed)) paste(shortBC[r$passed], collapse = ';') else '',
               min_anchor_p = min(r$step1[, 'p'], na.rm = TRUE),
               selected = if (is.na(a)) 'none' else if (is.na(i)) sig_label(a, NA, NA)
                          else sig_label(a, combo_grid$partner[i], combo_grid$sign[i]),
               anchor = if (is.na(a)) '' else a, partner = if (is.na(i)) '' else combo_grid$partner[i],
               sign = if (is.na(i)) NA_real_ else combo_grid$sign[i],
               chi2_best_combination = r$T_combo, chi2_best_anchor = r$T_anchor, chi2 = r$T, stringsAsFactors = FALSE)
  }

  if (length(todo)) {
    message(sprintf('Permutation: %d to run on %d workers; observed run took %.1f s, so roughly %.0f min.',
                    length(todo), N_WORKERS, sec_obs, length(todo) * sec_obs / N_WORKERS / 60))
    if (N_WORKERS > 1) plan(multisession, workers = N_WORKERS) else plan(sequential)
    for (chunk in split(todo, ceiling(seq_along(todo) / CHUNK))) {
      res = do.call(rbind, future_lapply(chunk, one_perm, future.seed = TRUE, future.packages = 'survival',
                                         future.globals = list(D = D, PERMS = PERMS, SETTINGS_KEY = SETTINGS_KEY,
                                                               run_pipeline = run_pipeline, int_test = int_test,
                                                               anchor_z = anchor_z, combo_z = combo_z, combo_grid = combo_grid,
                                                               combo_label = combo_label, sig_label = sig_label, shortBC = shortBC, lab80 = lab80,
                                                               F0 = F0, F1 = F1, INT = INT, GATE_P = GATE_P, ANCHOR_RULE = ANCHOR_RULE,
                                                               PARTNER_GATE = PARTNER_GATE)))
      done = rbind(done, res)
      write.csv(done, perm_file, row.names = FALSE)
      message(sprintf('  %d / %d done', nrow(done), N_PERM))
    }
    plan(sequential)
  }

  done = done[order(done$perm), ]
  null = done$chi2
  B = length(null)
  n_ge = if (T_obs > 0) sum(null >= T_obs) else B
  gate = done$n_passing > 0
  sel_freq = sort(table(done$selected[gate]), decreasing = TRUE)
  anchor_freq = table(factor(unlist(strsplit(done$passing[gate], ';')), levels = shortBC))
  perm_summary = data.frame(
    endpoint = ENDPOINT, n = nrow(D), events = sum(D$.status), gate_p = GATE_P, anchor_rule = ANCHOR_RULE,
    partner_rule = PARTNER_RULE, gsva_background = GSVA_BACKGROUND, observed_selected = sel_label,
    observed_chi2 = T_obs, observed_p_nominal = if (T_obs > 0) pchisq(T_obs, 1, lower.tail = FALSE) else NA,
    n_permutations = B, n_null_ge_observed = n_ge, whole_pipeline_p = (1 + n_ge) / (1 + B),
    mc_se = sqrt(((1 + n_ge) / (1 + B)) * (1 - (1 + n_ge) / (1 + B)) / B),
    mc_95ci_low = stats::binom.test(n_ge, B)$conf.int[1],     # exact (Clopper-Pearson) CI for the exceedance rate
    mc_95ci_high = stats::binom.test(n_ge, B)$conf.int[2],
    null_fraction_gate_passed = mean(gate),
    partner_gate = PARTNER_GATE,
    null_fraction_no_partner_if_gate_passed = if (any(gate)) mean(is.na(done$partner[gate]) | done$partner[gate] == '') else NA,
    null_fraction_final_nominal_p05 = mean(null >= qchisq(0.95, 1)),
    null_median_chi2_if_gate_passed = if (any(gate)) median(null[gate]) else NA,
    null_q95_chi2_all = unname(quantile(null, 0.95)), null_q99_chi2_all = unname(quantile(null, 0.99)),
    null_anchor_pass_counts = paste(names(anchor_freq), as.integer(anchor_freq), sep = ': ', collapse = '; '),
    null_top_selections = paste(head(names(sel_freq), 8), head(as.integer(sel_freq), 8), sep = ': ', collapse = '; '),
    perm_strata = if (is.null(PERM_STRATA)) 'none' else PERM_STRATA, seed = SEED, stringsAsFactors = FALSE)
  write.csv(perm_summary, summary_file, row.names = FALSE)

  ## figure: null distribution (runs that pass the gate) and the share that selects nothing
  pv = perm_summary$whole_pipeline_p
  pd = data.frame(chi2 = null[gate])
  pfig = ggplot(pd, aes(chi2)) +
    geom_histogram(bins = 40, fill = '#c3c2b7', colour = 'white', linewidth = 0.2) +
    { if (T_obs > 0) geom_vline(xintercept = T_obs, colour = '#e34948', linewidth = 0.9) } +
    annotate('text', x = if (T_obs > 0) T_obs else 0, y = Inf, vjust = 1.5, hjust = if (T_obs > 0) 1.08 else -0.1, size = 3.4,
             label = if (T_obs > 0) sprintf('observed %.1f\nwhole-pipeline p %s', T_obs,
                                            if (n_ge == 0) sprintf('<= %.2g', 1 / (1 + B)) else paste('=', sprintf('%.3f', pv)))
                     else 'observed: no anchor passed the gate') +
    labs(x = sprintf('Selected-combination TREAT x score chi2 (%s)', ENDPOINT), y = 'Permutations',
         title = paste0('Whole-pipeline permutation: BC360 anchor gate + direct-actor APM80 partner',
                        if (PARTNER_GATE) ' (partner gate)' else ''),
         subtitle = sprintf('%d permutations; %.0f%% passed the gate (shown); %.0f%% selected nothing (chi2 = 0).\nObserved: %s',
                            B, 100 * mean(gate), 100 * mean(!gate), sel_label)) +
    theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = 'bold'))
  ggsave(file.path(OUT_DIR, 'Fig_permutation_null.png'), pfig, width = 8, height = 4.8, dpi = 300, bg = 'white')

  cat('\n== Whole-pipeline permutation ==\n')
  print(t(perm_summary), quote = FALSE)
}

## manifest
write.csv(data.frame(run_timestamp = format(Sys.time(), tz = 'UTC', usetz = TRUE), endpoint = ENDPOINT,
                     n = nrow(D), events = sum(D$.status), adjust_terms = adjust_terms, gate_p = GATE_P,
                     anchor_rule = ANCHOR_RULE, partner_rule = PARTNER_RULE, partner_gate = PARTNER_GATE, min_genes = MIN_GENES,
                     gsva_background = GSVA_BACKGROUND, gsva_genes = nrow(X), gsva_kcdf = GSVA_KCDF,
                     gsva_version = as.character(utils::packageVersion('GSVA')), n_perm = N_PERM, seed = SEED,
                     perm_strata = if (is.null(PERM_STRATA)) 'none' else PERM_STRATA,
                     observed_selected = sel_label, observed_chi2 = T_obs),
          file.path(OUT_DIR, '00_run_manifest.csv'), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(OUT_DIR, 'sessionInfo.txt'))

cat('\n== Step 1: BC360 gate (GSVA) ==\n')
print(transform(s1[, c('label', 'class', 'chi2', 'p', 'int_HR_per_SD', 'passes_gate')],
                chi2 = round(chi2, 2), p = signif(p, 2), int_HR_per_SD = round(int_HR_per_SD, 2)), row.names = FALSE)
if (length(obs$passed)) {
  cat('\n== Step 3: combinations, best first ==\n')
  print(transform(head(g[, c('combination', 'direction', 'chi2', 'p', 'int_HR_per_SD', 'selected')], 14),
                  chi2 = round(chi2, 2), p = signif(p, 2), int_HR_per_SD = round(int_HR_per_SD, 2)), row.names = FALSE)
}
cat('\nOutputs in ', OUT_DIR, '\n', sep = '')
