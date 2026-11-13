## =============================================================================
## PANTHER APM score: synthetic input data for test runs of the analysis script
##
## Writes simulated versions of the non-public inputs of the main script into data/, where the main script reads them (see data/README.md). 
## The four files in the repository were made with this script (seed 2027); it is included to show how. Running it again overwrites them.
## Files written (under data/):
##   PANTHER/PANTHER_clinical.rds                  one row per randomized patient; sample_id set for the 444 tumours with RNA-seq data
##   PANTHER/PANTHER_expression_log2_TMM_TPM.rds   log2(TPM + 1), gene symbols x the 444 samples
##   TEX/pd.unique.mets.RData                      TEX relapse-biopsy clinical data
##   TEX/TEX_clinical incl microarray SEB.xlsx     SampleIDRelapse -> CEL barcode of the GSE56493 arrays
##
## The public inputs (BC360 gene list, Danaher annotation, SCAN-B, MATADOR, TEX arrays) are not generated: download them into data/ first. 
## Three are used here when available:
##   - the BC360 and Danaher gene lists, and the MCP-counter gene list (read from GitHub), are added to
##     the synthetic transcriptome, so every signature in the main script finds its genes;
##   - GSE56493 gives the real CEL barcodes, so the synthetic TEX clinical rows link to the public TEX
##     arrays exactly as the real ones do. Without it the two TEX files are not written.
##
## PANTHER expression values and all clinical characteristics and outcomes are simulated.
## Synthetic TEX records are randomly linked to publicly available GEO array barcodes.
## No private participant records or original clinical-to-array links are read or copied.
## A treatment x (HLA6 - IFN5) interaction is planted in BCRFS so that the signature definition (Part 1 of the main script) passes its anchor gate and selects a signature. 
## Results from these files have no biological meaning.
##
## Packages: here, openxlsx (TEX workbook); optional: readxl, curl, GEOquery, Biobase.
## =============================================================================

## -----------------------------------------------------------------------------
## 0. Settings
## -----------------------------------------------------------------------------
here::i_am("workflows/apm_analysis/scripts/PANTHER_make_synthetic_data.R")   # finds the repository root
data_dir  = here::here("data")   # public inputs are read from here; synthetic files are written here
SEED      = 2027
N_TRIAL   = 2003               # randomized patients (all PANTHER)
N_SWEDISH = 780                # Swedish patients: patno 10001.. (Table 1 uses patno < 12000); others 20001..
N_RNA     = 444                # Swedish patients with RNA-seq data that passed QC
N_GENES   = 5000               # synthetic transcriptome size (signature + annotation genes + background)
N_TEX     = 140                # TEX clinical rows; up to the number of GSE56493 arrays are linked
PLANT_LOGHR = c(`FEC-D` = 0.5, `tddEC-D` = -0.7)   # BCRFS log-HR per SD of the latent HLA6 - IFN5 contrast

set.seed(SEED)

paths = list(
  clinical     = "PANTHER/PANTHER_clinical.rds",
  expression   = "PANTHER/PANTHER_expression_log2_TMM_TPM.rds",
  tex_clinical = "TEX/pd.unique.mets.RData",
  tex_barcodes = "TEX/TEX_clinical incl microarray SEB.xlsx",
  # public inputs, read only if present
  bc360        = "annotation/LBL-10541-03-Breast-Cancer-360-Gene-List.xlsx",
  danaher      = "annotation/gannot-cip-junecelltype.csv",
  tex_geo_dir  = "TEX/GEO")

out = function(key) { f = file.path(data_dir, paths[[key]]); dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); f }
inp = function(key) file.path(data_dir, paths[[key]])

## -----------------------------------------------------------------------------
## 1. Gene universe
## -----------------------------------------------------------------------------
genes = list(
  hla6        = c("HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1"),
  mhc2_other  = c("HLA-DRB5", "HLA-DQA2", "HLA-DQB2", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CD74", "CIITA", "CTSS", "CTSL", 
                  "CTSB", "CTSD", "LGMN", "IFI30", "RFX5", "RFXANK", "RFXAP"),
  ifn5        = c("JAK1", "TYK2", "STAT1", "STAT2", "IRF9"),
  ifn_other   = c("IFNG", "IFNGR1", "IFNGR2", "JAK2", "IRF1", "IRF2", "NLRC5"),
  mhc1        = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "TAPBPL", "PSMB8", "PSMB9", "PSMB10",
                  "PSMB5", "PSMB6", "PSMB7", "PSME1", "PSME2", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "HLA-E", "HLA-F", "HLA-G"),
  apm_other   = c("NBR1", "MAL2", "SPPL3", "LNPEP", "RAB7A", "SEC22B", "CYBB", "DERL1", "DERL2", "VCP", "SEC61A1",
                  "EZH2", "SUZ12", "EED", "DNMT1", "DNMT3A", "DNMT3B", "HDAC1", "HDAC2", "HDAC3"),
  bc360_other = c("CD1E", "CD8A", "THBS1"),
  igg         = c("CD27", "CD79A", "CXCL8", "IGKC", "IGLC6", "IGLV3-25", "IL2RG", "JCHAIN", "LAX1", "NTN3", "PIM2", "POU2AF1", "TNFRSF17"),
  set_er      = c("ABAT", "ADCY1", "AZGP1", "CA12", "CD2", "CD3D", "DNAJC12", "KCNE4", "MAPT", "MRPS30", "NAT1", "NPY1R",
                  "PDZK1", "QDPR", "SCUBE2", "SLC39A6", "STC2", "BCL2"),
  reference   = c("AK2", "APPBP2", "ATP5J2", "DARS", "LDHA", "TRIM2", "UBE2Z", "UGP2", "VDAC2", "WIPF2", "ACTB", "GAPDH", "GUSB", "RPLP0", "TFRC"),
  esr1 = "ESR1", pgr = "PGR",
  her2        = c("ERBB2", "GRB7"),
  prolif      = c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS", "UBE2C", "ANLN",
                  "CCNE1", "CDC6", "CENPF", "EXO1", "KIF2C", "MELK", "ORC6", "UBE2T", "AURKA", "MYBL2"),
  other_rs    = c("MMP11", "CTSV", "CD68", "GSTM1", "BAG1"),
  immune      = c("CD3E", "CD3G", "CD247", "CD8B", "GZMA", "GZMB", "GZMK", "PRF1", "NKG7", "KLRK1", "CD19", "MS4A1", "CD22",
                  "CD14", "CD163", "CSF1R", "FCGR3A", "CD68", "PTPRC", "LCK", "ZAP70", "FOXP3", "CTLA4", "PDCD1", "CD274",
                  "LAG3", "TIGIT", "XCL1", "XCL2", "CCL5", "CXCL9", "CXCL10", "CXCL13"),
  stroma      = c("COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRB", "PECAM1", "VWF", "CDH5", "THBS1"))

# Gene lists from the public annotation files, if present in data/
public_genes = character(0)
if (file.exists(inp("bc360")) && requireNamespace("readxl", quietly = TRUE)) {
  bc = readxl::read_excel(inp("bc360"), sheet = 3, skip = 1)
  public_genes = c(public_genes, bc$Gene[!is.na(bc$Gene)])
}
danaher_genes = character(0)
if (file.exists(inp("danaher"))) {
  ga = utils::read.csv(inp("danaher"), row.names = 1)
  danaher_genes = rownames(ga)[!is.na(ga$Cell.Type.TCGA2)]
  public_genes = c(public_genes, rownames(ga))
}
mcp_genes = tryCatch({
  m = utils::read.table(curl::curl("https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),
                        sep = "\t", header = TRUE, colClasses = "character", check.names = FALSE)
  m$`HUGO symbols`[m$`HUGO symbols` == "C"] = "SDPR"
  m$`HUGO symbols`[m$`HUGO symbols` == "DIPK2B"] = "CXorf36"
  m$`HUGO symbols`
}, error = function(e) character(0))
public_genes = unique(c(public_genes, mcp_genes))

signature_genes = unique(unlist(genes, use.names = FALSE))
symbols = unique(c(signature_genes, public_genes))
symbols = symbols[!is.na(symbols) & nzchar(symbols)]
symbols = c(symbols, sprintf("SYNGENE%05d", seq_len(max(0, N_GENES - length(symbols)))))
G = length(symbols)
message(sprintf("Synthetic transcriptome: %d genes (%d signature genes, %d more from public gene lists)",
                G, length(signature_genes), length(setdiff(public_genes, signature_genes))))

## -----------------------------------------------------------------------------
## 2. All randomized patients: baseline characteristics
## -----------------------------------------------------------------------------
patno = c(10000 + seq_len(N_SWEDISH), 20000 + seq_len(N_TRIAL - N_SWEDISH))
n = length(patno)

age = round(pmin(65, pmax(18, rnorm(n, 50, 8))))
ERpos = ifelse(runif(n) < 0.75, "Yes", "No")
newHER2 = ifelse(runif(n) < 0.20, "Yes", "No")
histgrade = sample(c("I", "II", "III"), n, TRUE, prob = c(0.15, 0.45, 0.40))
histgrade[ERpos == "No"] = sample(c("II", "III"), sum(ERpos == "No"), TRUE, prob = c(0.2, 0.8))
posnodes = ifelse(runif(n) < 0.82, pmin(40, 1 + rgeom(n, 0.25)), 0)
tumdiam = round(pmin(120, pmax(3, exp(rnorm(n, log(25), 0.45)))))
# node-negative patients follow the eligibility rule (ER-negative, >= 20 mm, grade 3) or are younger than 35
nn = posnodes == 0
young = nn & runif(n) < 0.15
age[young] = sample(25:34, sum(young), TRUE)
ERpos[nn & !young] = "No"; tumdiam[nn & !young] = pmax(tumdiam[nn & !young], 20); histgrade[nn & !young] = "III"

mensstat = ifelse(age < 48, "Premenopausal", ifelse(age < 55, "Postmenopausal (0,5-5 years)", "Postmenopausal (> 5 years)"))
mensstat[runif(n) < 0.03] = "Uncertain"
optype = sample(c("Mastectomy", "Breast-conserving surgery"), n, TRUE, prob = c(0.6, 0.4))
TREAT = sample(rep(c("FEC-D", "tddEC-D"), length.out = n))
prolifvalue = round(pmin(95, pmax(1, rnorm(n, 25, 15))))   # Ki-67 (%)

# A few missing values, as in a trial database
histgrade[sample(n, 12)] = NA
newHER2[sample(n, 4)] = NA

# RNA-seq cohort: N_RNA Swedish patients; one of them without tumour grade (as in the real data)
rna_idx = sort(sample(which(patno < 12000), N_RNA))
histgrade[rna_idx][is.na(histgrade[rna_idx])] = "III"
newHER2[rna_idx][is.na(newHER2[rna_idx])] = "No"
histgrade[rna_idx[sample(N_RNA, 1)]] = NA
sample_id = rep(NA_character_, n)
sample_id[rna_idx] = sprintf("SYN.%04d.%d", seq_len(N_RNA), patno[rna_idx])

## -----------------------------------------------------------------------------
## 3. Expression of the RNA-seq cohort: log2(TPM + 1), gene symbols x samples
## -----------------------------------------------------------------------------
er_i   = as.numeric(ERpos[rna_idx] == "Yes")
her2_i = as.numeric(newHER2[rna_idx] == "Yes")
g3_i   = as.numeric(histgrade[rna_idx] %in% "III")

F = cbind(                                   # latent factors per tumour
  immune = rnorm(N_RNA),
  stroma = rnorm(N_RNA),
  er     = 2.2 * er_i + rnorm(N_RNA, 0, 0.5),
  her2   = 3.0 * her2_i + rnorm(N_RNA, 0, 0.4),
  prolif = 0.6 * g3_i - 0.6 * er_i + rnorm(N_RNA, 0, 0.8))
F = cbind(F,
  mhc2  = 0.6 * F[, "immune"] + 0.8 * rnorm(N_RNA),
  ifn   = 0.3 * F[, "immune"] + 0.95 * rnorm(N_RNA),
  mhc1  = 0.5 * F[, "immune"] + 0.85 * rnorm(N_RNA),
  bcell = 0.6 * F[, "immune"] + 0.8 * rnorm(N_RNA))

L = matrix(0, G, ncol(F), dimnames = list(symbols, colnames(F)))   # gene loadings on the factors
set_load = function(g, f, v) { g = intersect(g, symbols); L[g, f] <<- v }
set_load(genes$hla6, "mhc2", 0.9)
set_load(genes$mhc2_other, "mhc2", 0.15); set_load(genes$mhc2_other, "immune", 0.3)
set_load(genes$ifn5, "ifn", 0.8)
set_load(genes$ifn_other, "ifn", 0.3); set_load(genes$ifn_other, "immune", 0.5)
set_load(genes$mhc1, "mhc1", 0.7); set_load(genes$mhc1, "ifn", 0.1)
set_load(c(genes$immune, genes$bc360_other, mcp_genes, danaher_genes, "CD2", "CD3D"), "immune", 0.7)
set_load(genes$igg, "bcell", 1.0)
set_load(genes$set_er, "er", 0.6)
set_load(genes$esr1, "er", 1.2)
set_load(genes$her2, "her2", 1.0)
set_load(genes$prolif, "prolif", 0.8)
set_load(genes$stroma, "stroma", 0.8)

mu = pmin(13, pmax(0.5, rnorm(G, 5, 1.8))); names(mu) = symbols
high = intersect(c(genes$hla6, genes$mhc1, "CD74", "B2M"), symbols)
mu[high] = runif(length(high), 8, 11)
ref = intersect(genes$reference, symbols)
mu[ref] = runif(length(ref), 7, 10)
sigma = runif(G, 0.2, 0.6); names(sigma) = symbols
sigma[ref] = 0.12

X = matrix(mu, G, N_RNA) + L %*% t(F) + matrix(rnorm(G * N_RNA), G, N_RNA) * sigma
X["PGR", ] = ifelse(er_i == 1 & runif(N_RNA) < 0.75, rnorm(N_RNA, 6.5, 0.8), rnorm(N_RNA, 2, 0.8))   # bimodal PGR
X = pmax(X, 0)
tpm = 2^X - 1
X = log2(sweep(tpm, 2, colSums(tpm), "/") * 1e6 + 1)          # each sample sums to 1e6 TPM
dimnames(X) = list(symbols, sample_id[rna_idx])
X = X[order(rowSums(X), decreasing = TRUE), , drop = FALSE]   # row order as in the real matrix
saveRDS(as.data.frame(X), out("expression"))

## -----------------------------------------------------------------------------
## 4. Outcomes (all patients; planted interaction in the RNA-seq cohort)
## -----------------------------------------------------------------------------
# The planted contrast is computed from the simulated HLA6 and IFN5 expression (mean gene z-score of each set)
set_z = function(g) as.numeric(scale(colMeans(t(scale(t(X[g, , drop = FALSE]))))))
latent = rep(0, n)
latent[rna_idx] = as.numeric(scale(set_z(genes$hla6) - set_z(genes$ifn5)))

lp = 0.35 * log1p(posnodes) + 0.4 * (log(tumdiam) - log(25)) + 0.35 * (histgrade %in% "III") +
     0.35 * (ERpos == "No") + 0.15 * (newHER2 %in% "Yes") - 0.10 * (TREAT == "tddEC-D") +
     ifelse(TREAT == "FEC-D", PLANT_LOGHR[["FEC-D"]], PLANT_LOGHR[["tddEC-D"]]) * latent
t_bc      = stats::rexp(n, 0.035 * exp(lp - mean(lp)))   # breast cancer recurrence (years)
t_other   = stats::rexp(n, 0.006)                       # death without a breast cancer event
t_cens    = runif(n, 5, 12)                             # end of follow-up
distant   = runif(n) < 0.75                             # recurrence includes distant disease
t_contra  = stats::rexp(n, 0.004)                       # contralateral / second primary (EFS only)
t_bcdeath = t_bc + runif(n, 0.5, 4)                     # death after distant recurrence

# Time in days and 0/1 status; competing deaths are censored for BCRFS and DDFS
endpoint = function(t_event, t_compete, t_end) {
  tt = pmin(t_event, t_compete, t_end)
  list(time = round(tt * 365.25), status = as.numeric(t_event <= pmin(t_compete, t_end)))
}
bcrfs = endpoint(t_bc, t_other, t_cens)
ddfs  = endpoint(ifelse(distant, t_bc, Inf), t_other, t_cens)
efs   = endpoint(pmin(t_bc, t_contra, t_other), Inf, t_cens)
os    = endpoint(pmin(ifelse(distant & runif(n) < 0.6, t_bcdeath, Inf), t_other), Inf, t_cens)

## -----------------------------------------------------------------------------
## 5. PANTHER clinical table: the columns and types of the real PANTHER_clinical.rds
## -----------------------------------------------------------------------------
sub_prob = function(er, her2) {
  if (her2 == 1) c(LumA = 0.10, LumB = 0.30, Her2 = 0.50, Basal = 0.10)
  else if (er == 1) c(LumA = 0.45, LumB = 0.45, Her2 = 0.05, Basal = 0.05)
  else c(LumA = 0.03, LumB = 0.07, Her2 = 0.15, Basal = 0.75)
}
subtype = rep(NA_character_, n)
subtype[rna_idx] = vapply(seq_len(N_RNA), function(i) sample(c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like"), 1,
                                                             prob = sub_prob(er_i[i], her2_i[i])), "")
rna = seq_len(n) %in% rna_idx
status_rna = function(s) ifelse(rna, s, NA_real_)   # 0/1 status is set for the RNA-seq cohort only, as in the real file

clinical = data.frame(
  patno = as.numeric(patno), sample_id = sample_id, TREAT = TREAT, age = age,
  mensstat = factor(mensstat, levels = c("Premenopausal", "Postmenopausal (0,5-5 years)", "Postmenopausal (> 5 years)", "Uncertain")),
  optype = optype,
  tumdiam = tumdiam, tumdiam_cat = factor(ifelse(tumdiam <= 20, "≤ 2 cm", "> 2 cm"), levels = c("≤ 2 cm", "> 2 cm")),
  posnodes = posnodes,
  posnodes_cat = factor(cut(posnodes, c(-Inf, 0, 3, 9, Inf), labels = c("N0", "N1-3", "N4-9", "N10+")), levels = c("N0", "N1-3", "N4-9", "N10+")),
  histgrade = factor(histgrade, levels = c("I", "II", "III")),
  ERpos = factor(ERpos, levels = c("No", "Yes")), newHER2 = factor(newHER2, levels = c("No", "Yes")),
  Clinical_Subtype = factor(ifelse(newHER2 %in% "Yes", "HER2+", ifelse(is.na(newHER2), NA, ifelse(ERpos == "Yes", "ER+/HER2-", "TNBC"))),
                            levels = c("ER+/HER2-", "HER2+", "TNBC")),
  subtypecd_ssp.subtype = factor(subtype, levels = c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like")),
  prolifvalue = prolifvalue,
  timeBCRFS = bcrfs$time, BCRFS.status = status_rna(bcrfs$status),
  timeEFS = efs$time, EFS.status = status_rna(efs$status),
  timeDDFS = ddfs$time, DDFS.status = status_rna(ddfs$status),
  timeOS = os$time, OS.status = status_rna(os$status),
  stringsAsFactors = FALSE)
clinical = clinical[order(!rna, clinical$patno), ]   # RNA-seq cohort first, as in the real file
rownames(clinical) = NULL
saveRDS(clinical, out("clinical"))

## -----------------------------------------------------------------------------
## 6. TEX clinical data linked to the public GSE56493 arrays
## -----------------------------------------------------------------------------
# CEL barcodes from the GEO titles ("<n>_<barcode>") of the GPL10379 arrays
tex_geo_barcodes = function(geo_dir) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) return(NULL)
  dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)
  local_file = list.files(geo_dir, pattern = "GSE56493.*GPL10379.*series_matrix.*", full.names = TRUE)
  geo = tryCatch(
    if (length(local_file)) GEOquery::getGEO(filename = local_file[1], getGPL = FALSE)
    else {
      gl = GEOquery::getGEO("GSE56493", GSEMatrix = TRUE, getGPL = FALSE, destdir = geo_dir)
      gl[[which(vapply(gl, Biobase::annotation, character(1)) == "GPL10379")]]
    }, error = function(e) NULL)
  if (is.null(geo)) return(NULL)
  sub("^[0-9]+_", "", Biobase::pData(geo)$title)
}

write_tex = function(barcodes) {
  stopifnot(requireNamespace("openxlsx", quietly = TRUE))
  ids = sprintf("SYN%04d", seq_len(N_TEX))                      # synthetic SampleIDRelapse
  site_labels = c("Lymph node axilla", "Axillary lymphnode", "Lymph node fossa supraclav", "Lymph node neck",
                  "Skin thoracic wall", "Skin, breast", "Scar, skin", "Breast", "Right breast",
                  "Liver", "Lung", "Bone spine", "Ascites")
  site_prob = c(0.16, 0.08, 0.08, 0.04, 0.12, 0.05, 0.03, 0.10, 0.06, 0.12, 0.08, 0.06, 0.02)
  treat = c(sample(c("ET", "TEX"), N_TEX - 12, TRUE), rep("Not randomised", 12))
  tex = data.frame(SampleIDRelapse = ids, ClinTrialTreat = sample(treat), ER_PT = ifelse(runif(N_TEX) < 0.7, 1, 0),
                   RelapseSite = sample(site_labels, N_TEX, TRUE, site_prob), stringsAsFactors = FALSE)
  t_ev = stats::rexp(N_TEX, 0.35 * exp(0.4 * (tex$ER_PT == 0)))
  t_c = runif(N_TEX, 1, 8)
  tex$RelToDeath_cens = round(pmin(t_ev, t_c), 3)
  tex$BCSS_cens = as.integer(t_ev <= t_c)
  miss = sample(N_TEX, 6); tex$ER_PT[miss[1:3]] = NA; tex$RelapseSite[miss[4:6]] = NA
  pd.unique.mets = tex
  save(pd.unique.mets, file = out("tex_clinical"))

  # Barcode workbook: the first dummy patients get the arrays, in random order
  k = min(length(barcodes), N_TEX)
  link = data.frame(SampleIDRelapse = ids[seq_len(k)], BARCODE = paste0("@", sample(barcodes)[seq_len(k)]))
  wb = openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "Sheet1"); openxlsx::writeData(wb, "Sheet1", link)
  openxlsx::saveWorkbook(wb, out("tex_barcodes"), overwrite = TRUE)
}

barcodes = tex_geo_barcodes(inp("tex_geo_dir"))
if (is.null(barcodes)) {
  warning("GSE56493 not available (GEOquery missing, or no internet and no saved series matrix): TEX files not written.")
} else write_tex(barcodes)

## -----------------------------------------------------------------------------
## 7. Note in the data folder
## -----------------------------------------------------------------------------
writeLines(c(
  "SYNTHETIC DATA - simulated, not patient data.",
  sprintf("Made by PANTHER_make_synthetic_data.R on %s (seed %d).", format(Sys.Date()), SEED),
  "Covers the non-public inputs: PANTHER clinical table and expression matrix; TEX clinical data and barcode map.",
  "A treatment x (HLA6 - IFN5) interaction is planted in BCRFS so that the signature definition runs end to end.",
  "Results from these files have no biological meaning."),
  file.path(data_dir, "README_SYNTHETIC_DATA.txt"))

message("Synthetic data written to ", normalizePath(data_dir))
