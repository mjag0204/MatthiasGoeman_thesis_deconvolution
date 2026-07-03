# =============================================================================
# Thesis Analysis Script
# Epigenetic Deconvolution in Alzheimer's Disease — Entorhinal Cortex
#
# Dataset: GSE105109 (Smith et al. 2019, Clinical Epigenetics)
# Author:  Matthias Goeman
# Last update: 30/06/2026
#
# PIPELINE OVERVIEW
# -----------------
# Section 0:  PACKAGES
# Section 1:  LOAD DATA
# Section 2:  QC — ORIGINAL PAPER PIPELINE
# Section 3:  CELL-TYPE DECONVOLUTION
# Section 4:  HELPER FUNCTIONS
# Section 5:  PAPER DESIGN REPLICATION + ANK1 VALIDATION (original pipeline)
# Section 6:  RELOAD DATA
# Section 7:  QC — EXTENDED PIPELINE
# Section 8:  PAPER RESULT VERIFICATION + ANK1 VALIDATION (extended pipeline)
# Section 9:  EXTENDED ANALYSIS
# Section 10: STATISTICAL OBJECTIVE: P-VALUE DISTRIBUTION
# Section 11: SENSITIVITY 5mC VS 5hmC
# Section 12: BIOLOGICAL OBJECTIVE: PROBE DISCOVERY
# Section 13: SAVE ALL RESULTS
# =============================================================================


# =============================================================================
# SECTION 0. PACKAGES
# =============================================================================

# Install (run once)
# install.packages("BiocManager")
# BiocManager::install(c("minfi", "wateRmelon", "EpiDISH",
#   "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#   "IlluminaHumanMethylation450kmanifest", "limma", "missMethyl", "HiBED"))
# install.packages(c("GEOquery","ggplot2","tidyverse","pheatmap",
#   "ggrepel","VennDiagram","gridExtra","quadprog","MASS","reshape2","R.utils"))
# BiocManager::install("methylumi")

library(minfi)
library(wateRmelon)
library(EpiDISH)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(limma)
library(missMethyl)
library(GEOquery)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(gridExtra)
library(HiBED)
library(MASS)
library(dplyr)
library(deconvR)
library(reshape2)
library(MLML2R)


# =============================================================================
# SECTION 1. LOAD DATA
# =============================================================================

base_dir <- "C:/Users/matth/OneDrive/Documents/UGhent/Masterthesis"

# --- (1) Unpack IDATs --------------------------------------------------------
tar_file <- file.path(base_dir, "GSE105109_RAW.tar")
tar_out  <- file.path(base_dir, "GSE105109_RAW")
dir.create(tar_out, showWarnings = FALSE)
untar(tar_file, exdir = tar_out)

gz_files <- list.files(tar_out, pattern = "\\.idat\\.gz$",
                       full.names = TRUE, recursive = TRUE)
out_dir  <- file.path(base_dir, "GSE105109_idat")
dir.create(out_dir, showWarnings = FALSE)

for (f in gz_files) {
  R.utils::gunzip(filename = f,
                  destname  = file.path(out_dir, basename(sub("\\.gz$", "", f))),
                  overwrite = TRUE, remove = FALSE)
}

# --- (2) Metadata ------------------------------------------------------------
gse   <- getGEO("GSE105109", getGPL = FALSE)[[1]]
pheno <- pData(gse)

roman_to_number <- c("I"=1, "II"=2, "III"=3, "IV"=4, "V"=5, "VI"=6)
pheno_braak_raw <- gsub("braak stage: ", "", pheno$characteristics_ch1.3,
                        ignore.case = TRUE)
pheno_braak_all <- ifelse(pheno_braak_raw == "0", 0,
                          roman_to_number[pheno_braak_raw])
names(pheno_braak_all) <- rownames(pheno)

na_rows    <- which(is.na(pheno_braak_all))
gsm_remove <- rownames(pheno)[na_rows]
pheno      <- pheno[!rownames(pheno) %in% gsm_remove, ]
pheno_braak_all <- pheno_braak_all[!is.na(pheno_braak_all)]

idat_files <- list.files(out_dir, pattern = "\\.idat$", full.names = TRUE)
file.remove(idat_files[grepl(paste(gsm_remove, collapse = "|"),
                             basename(idat_files))])

# --- (3) Keep only Entorhinal Cortex samples ---------------------------------
is_ec         <- grepl("entorhinal cortex", pheno$title, ignore.case = TRUE)
gsm_remove_cb <- rownames(pheno)[!is_ec]
pheno         <- pheno[is_ec, ]
pheno_braak_all <- pheno_braak_all[rownames(pheno)]

idat_files <- list.files(out_dir, pattern = "\\.idat$", full.names = TRUE)
file.remove(idat_files[grepl(paste(gsm_remove_cb, collapse = "|"),
                             basename(idat_files))])
cat("EC samples retained:", nrow(pheno), "\n")

# --- (4) Harmonised sample IDs -----------------------------------------------
num  <- as.numeric(sub(".*_(\\d+)$", "\\1", pheno$title))
num3 <- sprintf("%03d", num)
pheno$sample_ID <- paste0("sample", num3, "_EC")
pheno$assay <- ifelse(grepl("_oxbs_", pheno$title, ignore.case = TRUE), "OxBS",
                      ifelse(grepl("_bs_", pheno$title, ignore.case = TRUE), "BS", NA))

conversion_table <- pheno[, c("geo_accession", "sample_ID", "assay")]
gsm_to_sample    <- setNames(conversion_table$sample_ID,
                             conversion_table$geo_accession)
sample_to_gsm    <- setNames(conversion_table$geo_accession,
                             conversion_table$sample_ID)

pheno_BS   <- pheno[pheno$assay == "BS",   ]
pheno_OxBS <- pheno[pheno$assay == "OxBS", ]

# --- (5) Read IDATs ----------------------------------------------------------
RGset_all    <- read.metharray.exp(base = out_dir, extended = TRUE)
rg_gsm_clean <- sub("_(.*)", "", colnames(RGset_all))
colnames(RGset_all) <- rg_gsm_clean

RGset_BS   <- RGset_all[, colnames(RGset_all) %in% rownames(pheno_BS)]
RGset_OxBS <- RGset_all[, colnames(RGset_all) %in% rownames(pheno_OxBS)]
cat("RGset_BS samples:",   ncol(RGset_BS),   "| expected ~95\n")
cat("RGset_OxBS samples:", ncol(RGset_OxBS), "| expected ~95\n")

extract_covar <- function(pheno_df, char_col, pattern, type = "character") {
  x <- trimws(gsub(pattern, "", pheno_df[[char_col]], ignore.case = TRUE))
  if (type == "numeric") as.numeric(x) else x
}

sex_BS   <- factor(tolower(extract_covar(pheno_BS,   "characteristics_ch1.1", "gender: ")))
sex_OxBS <- factor(tolower(extract_covar(pheno_OxBS, "characteristics_ch1.1", "gender: ")))
names(sex_BS)   <- rownames(pheno_BS)
names(sex_OxBS) <- rownames(pheno_OxBS)

age_BS   <- extract_covar(pheno_BS,   "characteristics_ch1.2", "age at death: *", type = "numeric")
age_OxBS <- extract_covar(pheno_OxBS, "characteristics_ch1.2", "age at death: *", type = "numeric")
names(age_BS)   <- rownames(pheno_BS)
names(age_OxBS) <- rownames(pheno_OxBS)

pheno_braak_BS   <- pheno_braak_all[rownames(pheno_BS)]
pheno_braak_OxBS <- pheno_braak_all[rownames(pheno_OxBS)]

# --- (6) Blacklists ----------------------------------------------------------
cr_file    <- file.path(base_dir, "48639-non-specific-probes-Illumina450k.csv")
multi_file <- file.path(base_dir, "HumanMethylation450_15017482_v.1.1_hg19_bowtie_multimap.txt")

cr_probes    <- trimws(as.character(read.csv(cr_file,     header = TRUE,  stringsAsFactors = FALSE)[, 1]))
multi_probes <- trimws(as.character(read.table(multi_file, header = FALSE, stringsAsFactors = FALSE)[, 1]))

anno     <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_snp <- anno[, c("CpG_maf", "SBE_maf")]

cat(sprintf("[Blacklists]  Chen XR: %d | BOWTIE2: %d\n",
            length(cr_probes), length(multi_probes)))


# =============================================================================
# SECTION 2. QC — ORIGINAL PAPER PIPELINE
# =============================================================================
# Steps (matching Smith et al. 2019):
#   (1) Sex concordance check
#   (2) Blacklist filtering: Chen XR + BOWTIE2 + SNP probes (MAF > 5%)
#   (3) dasen normalisation
#   (4) pfilter
# Applied separately to BS and OxBS data, followed by paired matching.
# =============================================================================

run_QC <- function(RGset, sex_vec, age_vec, braak_vec, pheno_df, assay_label) {
  
  qc_log   <- list()
  log_step <- function(step, probes, samples)
    qc_log[[length(qc_log) + 1]] <<- data.frame(
      Step = step, Probes = probes, Samples = samples, stringsAsFactors = FALSE)
  
  log_step("Input", nrow(RGset), ncol(RGset))
  
  # ── (1) Sex concordance check ───────────────────────────────────────────────
  sex_predicted  <- getSex(mapToGenome(RGset))
  sex_pred_clean <- ifelse(sex_predicted$predictedSex == "M", "m", "f")
  sex_reported   <- as.character(sex_vec[colnames(RGset)])
  mismatch       <- sex_pred_clean != sex_reported
  cat(sprintf("[Sex check]  Mismatches: %d\n", sum(mismatch)))
  
  sex_df <- data.frame(xMed = sex_predicted$xMed, yMed = sex_predicted$yMed,
                       Reported = sex_reported, Mismatch = mismatch)
  
  p_sex <- ggplot(sex_df, aes(x = xMed, y = yMed, colour = Reported, shape = Mismatch)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_colour_manual(values = c("m" = "steelblue", "f" = "tomato"),
                        labels = c("m" = "Reported male", "f" = "Reported female")) +
    scale_shape_manual(values = c("FALSE" = 19, "TRUE" = 4),
                       labels = c("FALSE" = "Pass", "TRUE" = "Mismatch")) +
    theme_minimal(base_size = 11) +
    labs(title    = paste0("Sex concordance check — ", assay_label,
                           " (original pipeline)"),
         subtitle = sprintf("%d mismatch(es) removed", sum(mismatch)),
         x        = "X chromosome median log2 intensity",
         y        = "Y chromosome median log2 intensity",
         colour   = "Reported sex", shape = "QC status")
  ggsave(file.path(base_dir, paste0("QC_sex_check_", assay_label, ".jpeg")),
         p_sex, width = 7, height = 5, dpi = 300, device = "jpeg")
  
  RGset     <- RGset[,    !mismatch]
  sex_vec   <- sex_vec[   colnames(RGset)]
  age_vec   <- age_vec[   colnames(RGset)]
  braak_vec <- braak_vec[ colnames(RGset)]
  pheno_df  <- pheno_df[  colnames(RGset), ]
  cat(sprintf("[Sex check]  Samples remaining: %d\n", ncol(RGset)))
  log_step("After sex check", nrow(RGset), ncol(RGset))
  
  cat("[Pre-extract]  Computing detection p-values and beadcounts ...\n")
  det_p  <- detectionP(RGset)
  bc_mat <- beadcount(RGset)
  
  # ── (2) Blacklist filtering ─────────────────────────────────────────────────
  cat("[Blacklist]  Filtering probes ...\n")
  MSet_raw <- preprocessRaw(RGset)
  log_step("preprocessRaw", nrow(MSet_raw), ncol(MSet_raw))
  
  MSet_raw <- MSet_raw[!rownames(MSet_raw) %in% cr_probes, ]
  log_step("After Chen XR blacklist", nrow(MSet_raw), ncol(MSet_raw))
  
  MSet_raw <- MSet_raw[!rownames(MSet_raw) %in% multi_probes, ]
  log_step("After BOWTIE2 blacklist", nrow(MSet_raw), ncol(MSet_raw))
  
  anno_sub   <- anno_snp[rownames(MSet_raw), ]
  snp_probes <- rownames(anno_sub)[which(anno_sub$CpG_maf > 0.05 |
                                           anno_sub$SBE_maf > 0.05)]
  MSet_raw   <- MSet_raw[!rownames(MSet_raw) %in% snp_probes, ]
  log_step("After SNP filter (MAF > 5%)", nrow(MSet_raw), ncol(MSet_raw))
  
  # ── (3) dasen normalisation ─────────────────────────────────────────────────
  cat("[dasen]  Running normalisation ...\n")
  b_raw     <- getBeta(MSet_raw)
  MSet_norm <- dasen(MSet_raw)
  b_values  <- getBeta(MSet_norm)
  M_values  <- getM(MSet_norm)
  log_step("After dasen normalisation", nrow(b_values), ncol(b_values))
  
  # ── (4) pfilter ─────────────────────────────────────────────────────────────
  cat("[pfilter]  Aligning QC matrices ...\n")
  common_probes  <- Reduce(intersect, list(rownames(b_values), rownames(det_p), rownames(bc_mat)))
  common_samples <- Reduce(intersect, list(colnames(b_values), colnames(det_p), colnames(bc_mat)))
  
  pf_result <- pfilter(
    mn = b_values[common_probes, common_samples],
    bn = b_values[common_probes, common_samples],
    pn = det_p[   common_probes, common_samples],
    bc = bc_mat[  common_probes, common_samples],
    perc = 5, pthresh = 1, perCount = 5, pnthresh = 0.05, logical.return = FALSE)
  
  b_values <- as.matrix(pf_result$bn)
  log_step("After pfilter", nrow(b_values), ncol(b_values))
  
  b_clip   <- pmax(pmin(b_values, 1 - 1e-6), 1e-6)
  M_values <- log2(b_clip / (1 - b_clip))
  
  surviving <- colnames(b_values)
  sex_vec   <- sex_vec[   surviving]
  age_vec   <- age_vec[   surviving]
  braak_vec <- braak_vec[ surviving]
  pheno_df  <- pheno_df[  surviving, ]
  
  qc_table <- do.call(rbind, qc_log)
  cat(sprintf("\n[%s FINAL]  Samples: %d | Probes: %d\n",
              assay_label, ncol(b_values), nrow(b_values)))
  cat("           Paper target: ~91 BS | ~367,480 probes\n")
  cat("\n[QC tracking table —", assay_label, "(original pipeline)]\n")
  print(qc_table)
  
  list(b_raw = b_raw, b_values = b_values, M_values = M_values,
       pheno = pheno_df, sex = sex_vec, age = age_vec,
       braak = braak_vec, qc_table = qc_table)
}

# Run QC on BS data -----------------------------------------------------------
cat("Running QC on BS data (original pipeline)...\n")
qc_BS          <- run_QC(RGset_BS, sex_BS, age_BS, pheno_braak_BS, pheno_BS, "BS")
b_values_BS    <- qc_BS$b_values;  M_values_BS    <- qc_BS$M_values
pheno_BS       <- qc_BS$pheno;     sex_BS         <- qc_BS$sex
age_BS         <- qc_BS$age;       pheno_braak_BS <- qc_BS$braak
qc_table_BS    <- qc_BS$qc_table

# Run QC on OxBS data ---------------------------------------------------------
cat("\nRunning QC on OxBS data (original pipeline)...\n")
qc_OxBS          <- run_QC(RGset_OxBS, sex_OxBS, age_OxBS, pheno_braak_OxBS, pheno_OxBS, "OxBS")
b_values_OxBS    <- qc_OxBS$b_values;  M_values_OxBS    <- qc_OxBS$M_values
pheno_OxBS       <- qc_OxBS$pheno;     sex_OxBS         <- qc_OxBS$sex
age_OxBS         <- qc_OxBS$age;       pheno_braak_OxBS <- qc_OxBS$braak
qc_table_OxBS    <- qc_OxBS$qc_table

cat(sprintf("\n[QC Summary — original pipeline]\n"))
cat(sprintf("  BS:   %d samples | %d probes\n", ncol(b_values_BS),   nrow(b_values_BS)))
cat(sprintf("  OxBS: %d samples | %d probes\n", ncol(b_values_OxBS), nrow(b_values_OxBS)))

# Beta distribution plot: raw vs dasen-normalised -----------------------------
cat("Generating beta distribution comparison plot (original pipeline)...\n")
set.seed(123)
df_distribution <- rbind(
  data.frame(Beta  = as.vector(qc_BS$b_raw[  sample(nrow(qc_BS$b_raw),   5000), ]),
             Assay = "BS",   Stage = "Raw"),
  data.frame(Beta  = as.vector(qc_OxBS$b_raw[sample(nrow(qc_OxBS$b_raw), 5000), ]),
             Assay = "OxBS", Stage = "Raw"),
  data.frame(Beta  = as.vector(b_values_BS[  sample(nrow(b_values_BS),   5000), ]),
             Assay = "BS",   Stage = "dasen-normalised"),
  data.frame(Beta  = as.vector(b_values_OxBS[sample(nrow(b_values_OxBS), 5000), ]),
             Assay = "OxBS", Stage = "dasen-normalised"))
df_distribution <- df_distribution[!is.na(df_distribution$Beta), ]

p_global_dist <- ggplot(df_distribution, aes(x = Beta, fill = Stage, colour = Stage)) +
  geom_density(alpha = 0.25, linewidth = 0.7, adjust = 1.2) +
  facet_wrap(~ Assay, nrow = 2) +
  scale_colour_manual(values = c("Raw" = "#e41a1c", "dasen-normalised" = "#377eb8")) +
  scale_fill_manual(  values = c("Raw" = "#e41a1c", "dasen-normalised" = "#377eb8")) +
  theme_bw(base_size = 12) +
  theme(strip.background = element_rect(fill = "grey95"),
        strip.text       = element_text(face = "bold"),
        legend.position  = "bottom") +
  labs(title    = "Effect of dasen normalisation on beta-value distributions",
       subtitle = "Random sample of 5,000 probes per stage (original pipeline)",
       x        = "Beta value (0 = unmethylated, 1 = methylated)",
       y        = "Density",
       fill = "Processing stage", colour = "Processing stage")
ggsave(file.path(base_dir, "QC_global_beta_normalisation_comparison.jpeg"),
       plot = p_global_dist, width = 8, height = 6, dpi = 300, device = "jpeg")
cat("[Saved]  QC_global_beta_normalisation_comparison.jpeg\n")

# Paired matching and probe intersection --------------------------------------
colnames(b_values_BS)   <- gsm_to_sample[colnames(b_values_BS)]
colnames(M_values_BS)   <- gsm_to_sample[colnames(M_values_BS)]
rownames(pheno_BS)      <- pheno_BS$sample_ID
names(sex_BS)           <- gsm_to_sample[names(sex_BS)]
names(age_BS)           <- gsm_to_sample[names(age_BS)]
names(pheno_braak_BS)   <- gsm_to_sample[names(pheno_braak_BS)]

colnames(b_values_OxBS) <- gsm_to_sample[colnames(b_values_OxBS)]
colnames(M_values_OxBS) <- gsm_to_sample[colnames(M_values_OxBS)]
rownames(pheno_OxBS)    <- pheno_OxBS$sample_ID
names(sex_OxBS)         <- gsm_to_sample[names(sex_OxBS)]
names(age_OxBS)         <- gsm_to_sample[names(age_OxBS)]
names(pheno_braak_OxBS) <- gsm_to_sample[names(pheno_braak_OxBS)]

matched_samples <- intersect(colnames(b_values_BS), colnames(b_values_OxBS))
cat(sprintf("[Pairing]  Matched donors: %d\n", length(matched_samples)))

for (obj in c("b_values_BS", "M_values_BS", "b_values_OxBS", "M_values_OxBS")) {
  m <- get(obj); assign(obj, m[, sort(matched_samples)])
}
pheno_BS         <- pheno_BS[        sort(matched_samples), ]
sex_BS           <- sex_BS[          sort(matched_samples)]
age_BS           <- age_BS[          sort(matched_samples)]
pheno_braak_BS   <- pheno_braak_BS[  sort(matched_samples)]
pheno_OxBS       <- pheno_OxBS[      sort(matched_samples), ]
sex_OxBS         <- sex_OxBS[        sort(matched_samples)]
age_OxBS         <- age_OxBS[        sort(matched_samples)]
pheno_braak_OxBS <- pheno_braak_OxBS[sort(matched_samples)]

matched_probes <- intersect(rownames(b_values_BS), rownames(b_values_OxBS))
cat(sprintf("[Probes]   Shared probes: %d\n", length(matched_probes)))

b_values_BS   <- b_values_BS[  sort(matched_probes), ]
M_values_BS   <- M_values_BS[  sort(matched_probes), ]
b_values_OxBS <- b_values_OxBS[sort(matched_probes), ]
M_values_OxBS <- M_values_OxBS[sort(matched_probes), ]

stopifnot(identical(colnames(b_values_BS), colnames(b_values_OxBS)))
stopifnot(identical(rownames(b_values_BS), rownames(b_values_OxBS)))
cat("[Sanity check]  BS and OxBS matrices aligned. Ready for decomposition.\n")

# Subtraction decomposition ---------------------------------------------------
b_values_BS_sub   <- as.matrix(b_values_BS[,   sort(matched_samples)])
b_values_OxBS_sub <- as.matrix(b_values_OxBS[, sort(matched_samples)])

BETA_THRESH      <- 0.1
keep_subtraction <- rowMeans(b_values_BS_sub,   na.rm = TRUE) > BETA_THRESH &
  rowMeans(b_values_OxBS_sub, na.rm = TRUE) > BETA_THRESH
cat(sprintf("[Beta filter %.2f]  Probes retained: %d\n",
            BETA_THRESH, sum(keep_subtraction)))

b_final_BS   <- b_values_BS_sub[  keep_subtraction, ]
b_final_OxBS <- b_values_OxBS_sub[keep_subtraction, ]

mC_matrix  <- b_final_OxBS
hmC_matrix <- b_final_BS - b_final_OxBS
uC_matrix  <- 1 - b_final_BS

prop_to_M <- function(prop, offset = 0.0001) {
  prop <- pmax(pmin(as.matrix(prop), 1 - offset), offset)
  log2(prop / (1 - prop))
}

M_5mC   <- prop_to_M(mC_matrix)
M_5hmC  <- prop_to_M(hmC_matrix)
M_uC    <- prop_to_M(uC_matrix)
M_total <- prop_to_M(b_values_BS)

hmC_detected    <- rowMeans(hmC_matrix > 0, na.rm = TRUE) > 0.5
hmC_matrix_filt <- hmC_matrix[hmC_detected, ]
M_5hmC_filt     <- M_5hmC[    hmC_detected, ]
cat(sprintf("[5hmC filter]  Probes present in >50%% of samples: %d\n",
            sum(hmC_detected)))

# Covariate alignment ---------------------------------------------------------
braak_mlml <- pheno_braak_BS[matched_samples]
age_mlml   <- age_BS[        matched_samples]
sex_mlml   <- sex_BS[        matched_samples]

cat("braak_mlml — range:", paste(range(braak_mlml, na.rm = TRUE), collapse = " to "),
    "| NAs:", sum(is.na(braak_mlml)), "\n")
cat("age_mlml   — range:", paste(range(age_mlml,   na.rm = TRUE), collapse = " to "),
    "| NAs:", sum(is.na(age_mlml)), "\n")
cat("sex_mlml   — table:", paste(names(table(sex_mlml)), table(sex_mlml),
                                 sep = "=", collapse = ", "), "\n")


# =============================================================================
# SECTION 3. CELL-TYPE DECONVOLUTION (HiBED, 5 methods)
# =============================================================================
# Reference: HiBED Layer 2B (5 brain cell types).
# Paper replication: NNLS neuron proportion as single covariate.
# Extended analysis: all cell types minus Endothelial and Stromal (reference).
# =============================================================================

data("HiBED_Libraries")
ref_se  <- HiBED_Libraries[3]$Library_Layer2B
ref_mat <- assay(ref_se)
cat("HiBED cell types:", paste(colnames(ref_mat), collapse = ", "), "\n")

b_bulk      <- as.matrix(b_values_BS)
common_cpgs <- intersect(rownames(b_bulk), rownames(ref_mat))
cat("Shared probes (bulk vs HiBED):", length(common_cpgs), "\n")

bulk_sub <- b_bulk[  common_cpgs, ]
ref_sub  <- ref_mat[ common_cpgs, ]

bulk_dec <- bulk_sub %>% as.data.frame() %>%
  tibble::rownames_to_column("IDs") %>% dplyr::select(IDs, everything())
ref_dec  <- ref_sub  %>% as.data.frame() %>%
  tibble::rownames_to_column("IDs") %>% dplyr::select(IDs, everything())

extract_props <- function(x) as.data.frame(x$proportions)

cat("Running NNLS (deconvR)...\n"); props_nnls <- extract_props(deconvolute(reference = ref_dec, bulk = bulk_dec, model = "nnls"))
cat("Running SVR  (deconvR)...\n"); props_svr  <- extract_props(deconvolute(reference = ref_dec, bulk = bulk_dec, model = "svr"))
cat("Running QP   (deconvR)...\n"); props_qp   <- extract_props(deconvolute(reference = ref_dec, bulk = bulk_dec, model = "qp"))
cat("Running RPC  (EpiDISH)...\n"); props_rpc  <- as.data.frame(epidish(beta.m = bulk_sub, ref.m = ref_sub, method = "RPC")$estF)
cat("Running CP   (EpiDISH)...\n"); props_cp   <- as.data.frame(epidish(beta.m = bulk_sub, ref.m = ref_sub, method = "CP")$estF)

cat("\nMean row sums (should be ~1):\n")
for (nm in c("nnls", "svr", "qp", "rpc", "cp"))
  cat(sprintf("  %-4s: %.3f\n", toupper(nm), mean(rowSums(get(paste0("props_", nm))))))

# Stacked bar plots -----------------------------------------------------------
plot_stacked_bars <- function(props_df, method_name, sample_ids) {
  df <- props_df[sample_ids, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Sample") %>%
    tidyr::pivot_longer(-Sample, names_to = "CellType", values_to = "Proportion") %>%
    dplyr::mutate(Sample = factor(Sample, levels = sample_ids))
  ggplot(df, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_y_continuous(limits = c(0, 1.001), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          panel.grid      = element_blank(),
          legend.position = "right") +
    labs(title    = paste0("Cell-type proportions — ", method_name),
         subtitle = paste0(length(sample_ids), " samples (HiBED Layer 2B)"),
         x = "Sample", y = "Proportion", fill = "Cell type")
}

all_sids <- rownames(props_nnls)

jpeg(file.path(base_dir, "results_deconv_stacked_bars.jpeg"),
     width = 14, height = 12, units = "in", res = 300, quality = 95)
grid.arrange(
  plot_stacked_bars(props_nnls, "NNLS (deconvR)", all_sids),
  plot_stacked_bars(props_svr,  "SVR  (deconvR)", all_sids),
  plot_stacked_bars(props_qp,   "QP   (deconvR)", all_sids),
  plot_stacked_bars(props_rpc,  "RPC  (EpiDISH)", all_sids),
  plot_stacked_bars(props_cp,   "CP   (EpiDISH)", all_sids),
  ncol = 1)
dev.off()
cat("[Saved]  results_deconv_stacked_bars.jpeg\n")

# Neuron proportion -----------------------------------------------------------
find_neuron_col <- function(props_df)
  grep("NeuN_pos$|neuron|NeuN(?!.*neg)", colnames(props_df),
       ignore.case = TRUE, perl = TRUE, value = TRUE)[1]

neuron_col_nnls <- find_neuron_col(props_nnls)
neuron_col_rpc  <- find_neuron_col(props_rpc)
cat(sprintf("\nNeuron column (NNLS): %s\n", neuron_col_nnls))
cat(sprintf("Neuron column (RPC):  %s\n",  neuron_col_rpc))

neuron_prop <- props_nnls[matched_samples, neuron_col_nnls]
names(neuron_prop) <- matched_samples
cat(sprintf("Neuron proportion — mean: %.3f | SD: %.3f | range: %.3f-%.3f\n",
            mean(neuron_prop), sd(neuron_prop),
            min(neuron_prop),  max(neuron_prop)))

# Pairwise neuron proportion alignment ----------------------------------------
neuron_comparison <- data.frame(
  NNLS = props_nnls[matched_samples, find_neuron_col(props_nnls)],
  SVR  = props_svr[ matched_samples, find_neuron_col(props_svr)],
  QP   = props_qp[  matched_samples, find_neuron_col(props_qp)],
  RPC  = props_rpc[ matched_samples, find_neuron_col(props_rpc)],
  CP   = props_cp[  matched_samples, find_neuron_col(props_cp)]
)

cat("\nNeuron proportion correlations:\n")
print(round(cor(neuron_comparison), 3))

jpeg(file.path(base_dir, "results_neuron_proportion_correlation.jpeg"),
     width = 7, height = 7, units = "in", res = 300, quality = 95)
pairs(neuron_comparison,
      main = "Neuron proportion: pairwise method comparison",
      pch  = 16, col = rgb(0, 0, 0, 0.3))
dev.off()
cat("[Saved]  results_neuron_proportion_correlation.jpeg\n")

# Multi-cell-type proportions (reference category dropped) --------------------
ref_col <- "Endothelial and Stromal"

drop_reference <- function(props_df, sids, method_name) {
  m <- as.matrix(props_df[sids, , drop = FALSE])
  if (!ref_col %in% colnames(m)) {
    cat(sprintf("  WARNING: '%s' not found in %s — dropping last column: %s\n",
                ref_col, method_name, colnames(m)[ncol(m)]))
    return(m[, -ncol(m), drop = FALSE])
  }
  cat(sprintf("  %-6s dropping reference: '%s'\n", method_name, ref_col))
  m[, colnames(m) != ref_col, drop = FALSE]
}

cat("\nDropping reference cell type:", ref_col, "\n")
cells_nnls <- drop_reference(props_nnls, matched_samples, "NNLS")
cells_svr  <- drop_reference(props_svr,  matched_samples, "SVR")
cells_qp   <- drop_reference(props_qp,   matched_samples, "QP")
cells_rpc  <- drop_reference(props_rpc,  matched_samples, "RPC")
cells_cp   <- drop_reference(props_cp,   matched_samples, "CP")

cat("Cell types retained in models:\n")
cat("  deconvR:", paste(colnames(cells_nnls), collapse = ", "), "\n")
cat("  EpiDISH:", paste(colnames(cells_rpc),  collapse = ", "), "\n")


# =============================================================================
# SECTION 4. HELPER FUNCTIONS
# =============================================================================

method_labels <- c(
  no   = "No deconvolution",
  nnls = "NNLS (deconvR)",
  svr  = "SVR  (deconvR)",
  qp   = "QP   (deconvR)",
  rpc  = "RPC  (EpiDISH)",
  cp   = "CP   (EpiDISH)"
)
method_colours <- c(
  "No deconvolution" = "black",
  "NNLS (deconvR)"   = "#E41A1C",
  "SVR  (deconvR)"   = "#FF7F00",
  "QP   (deconvR)"   = "#4DAF4A",
  "RPC  (EpiDISH)"   = "#377EB8",
  "CP   (EpiDISH)"   = "#984EA3"
)

P_EW    <- 2.4e-7   # experiment-wide threshold (Saffari et al. 2018)
P_STR   <- 0.05     # BH-adjusted p-value threshold (extended analysis)
LFC_STR <- 0.20     # minimum |logFC| threshold (extended analysis)

# ANK1 probes reported in Smith et al. 2019
ANK1_PROBES <- c("cg05066959", "cg11823178")

run_ewas <- function(M_mat, design_mat, coef_name) {
  ok      <- apply(M_mat, 1, function(x) all(is.finite(x)))
  if (sum(!ok) > 0) cat("    Non-finite probes removed:", sum(!ok), "\n")
  M_clean <- M_mat[ok, , drop = FALSE]
  fit     <- eBayes(lmFit(M_clean, design_mat))
  tt      <- topTable(fit, coef = coef_name, number = Inf,
                      adjust.method = "BH", sort.by = "p")
  tt$sig_ew  <- tt$P.Value   < P_EW
  tt$sig_str <- tt$adj.P.Val < P_STR & abs(tt$logFC) > LFC_STR
  tt
}

run_paper_ewas <- function(M_mat, braak_v, age_v, sex_v, neuron_v, label = "") {
  sids     <- colnames(M_mat)
  braak_v  <- braak_v[ sids]; age_v    <- age_v[   sids]
  sex_v    <- sex_v[   sids]; neuron_v <- neuron_v[sids]
  cat(sprintf("  [%s]  %d samples | %d probes\n", label, ncol(M_mat), nrow(M_mat)))
  design     <- model.matrix(~ braak_v + age_v + sex_v + neuron_v)
  braak_coef <- grep("braak_v", colnames(design), value = TRUE)[1]
  cat("  Braak coefficient:", braak_coef, "\n")
  run_ewas(M_mat, design, braak_coef)
}

run_all_designs <- function(M_mat, braak_vec, age_vec, sex_vec,
                            cells_nnls_m, cells_svr_m, cells_qp_m,
                            cells_rpc_m,  cells_cp_m,
                            coef_pattern, label_prefix) {
  sids    <- colnames(M_mat)
  cn_nnls <- cells_nnls_m[sids, , drop = FALSE]
  cn_svr  <- cells_svr_m[ sids, , drop = FALSE]
  cn_qp   <- cells_qp_m[  sids, , drop = FALSE]
  cn_rpc  <- cells_rpc_m[ sids, , drop = FALSE]
  cn_cp   <- cells_cp_m[  sids, , drop = FALSE]
  
  d_no   <- model.matrix(~ braak_vec + age_vec + sex_vec)
  d_nnls <- model.matrix(~ braak_vec + age_vec + sex_vec + cn_nnls)
  d_svr  <- model.matrix(~ braak_vec + age_vec + sex_vec + cn_svr)
  d_qp   <- model.matrix(~ braak_vec + age_vec + sex_vec + cn_qp)
  d_rpc  <- model.matrix(~ braak_vec + age_vec + sex_vec + cn_rpc)
  d_cp   <- model.matrix(~ braak_vec + age_vec + sex_vec + cn_cp)
  
  cn <- grep(coef_pattern, colnames(d_no), value = TRUE)
  if (length(cn) == 0) stop(label_prefix, ": pattern '", coef_pattern,
                            "' not found in: ", paste(colnames(d_no), collapse = ", "))
  if (length(cn) > 1) cn <- cn[1]
  cat("  [", label_prefix, "] coefficient:", cn, "\n")
  
  cat("  No deconvolution...\n"); r_no   <- run_ewas(M_mat, d_no,   cn)
  cat("  NNLS...\n");             r_nnls <- run_ewas(M_mat, d_nnls, cn)
  cat("  SVR...\n");              r_svr  <- run_ewas(M_mat, d_svr,  cn)
  cat("  QP...\n");               r_qp   <- run_ewas(M_mat, d_qp,   cn)
  cat("  RPC...\n");              r_rpc  <- run_ewas(M_mat, d_rpc,  cn)
  cat("  CP...\n");               r_cp   <- run_ewas(M_mat, d_cp,   cn)
  
  list(no = r_no, nnls = r_nnls, svr = r_svr, qp = r_qp, rpc = r_rpc, cp = r_cp)
}

pearson_r_table <- function(res_list) {
  shared <- Reduce(intersect, lapply(res_list, rownames))
  base   <- res_list$no[shared, "logFC"]
  sapply(res_list[c("nnls","svr","qp","rpc","cp")], function(m)
    round(cor(base, m[shared, "logFC"], use = "complete.obs"), 4))
}

eff_scatter <- function(tt1, tt2, lab1, lab2, title_str) {
  shared <- intersect(rownames(tt1), rownames(tt2))
  df     <- data.frame(x = tt1[shared, "logFC"], y = tt2[shared, "logFC"])
  r      <- round(cor(df$x, df$y, use = "complete.obs"), 3)
  ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.15, size = 0.4, colour = "steelblue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
    theme_minimal(base_size = 11) +
    labs(title    = title_str,
         subtitle = paste0("Pearson r = ", r),
         x = paste("logFC —", lab1), y = paste("logFC —", lab2))
}

sig_summary_str <- function(tt, label) {
  data.frame(Model     = label,
             N_probes  = nrow(tt),
             Str_total = sum(tt$sig_str),
             Str_hyper = sum(tt$sig_str & tt$logFC > 0),
             Str_hypo  = sum(tt$sig_str & tt$logFC < 0),
             stringsAsFactors = FALSE)
}

summarise_str <- function(res_list, prefix) {
  bind_rows(lapply(names(method_labels), function(nm)
    sig_summary_str(res_list[[nm]], paste0(prefix, "_", toupper(nm)))))
}

build_pval_df <- function(res_list) {
  shared <- Reduce(intersect, lapply(res_list, rownames))
  bind_rows(lapply(names(method_labels), function(nm) {
    data.frame(Method = method_labels[nm],
               rawP   = res_list[[nm]][shared, "P.Value"],
               adjP   = res_list[[nm]][shared, "adj.P.Val"],
               logFC  = res_list[[nm]][shared, "logFC"],
               stringsAsFactors = FALSE)
  })) %>% dplyr::mutate(Method = factor(Method, levels = method_labels))
}

volcano_plot_str <- function(tt, title_str) {
  top20    <- rownames(head(tt[order(tt$P.Value), ], 20))
  tt$label <- ifelse(tt$sig_str & rownames(tt) %in% top20, rownames(tt), NA)
  ggplot(tt, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(colour = sig_str), alpha = 0.4, size = 0.6) +
    geom_hline(yintercept = -log10(P_STR),  linetype = "dashed", colour = "grey50") +
    geom_vline(xintercept = c(-LFC_STR, LFC_STR), linetype = "dashed", colour = "grey50") +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5,
                             max.overlaps = 30, na.rm = TRUE) +
    scale_colour_manual(values = c("FALSE" = "grey80", "TRUE" = "red"),
                        name   = "Significant") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank()) +
    labs(title    = title_str,
         subtitle = paste0(sum(tt$sig_str),
                           " CpGs significant (adj.P < 0.05, |logFC| > 0.20)"),
         x = "logFC (M-value)", y = "-log10(adj.P)")
}

# ANK1 validation helper — reused in Sections 5 and 8 ------------------------
# Accepts pre-computed paper EWAS tables (tt_5mC, tt_5hmC) and the proportion
# matrices from the same pipeline. bin_s and braak_bin_f define the Braak
# binary grouping for the boxplot.
run_ank1_validation <- function(tt_5mC, tt_5hmC, mC_mat, hmC_mat,
                                braak_bin_f, bin_s, pipeline_label) {
  
  cat(sprintf("\n===== ANK1 validation — %s =====\n", pipeline_label))
  
  # Console report: EWAS statistics for both probes
  for (mod_label in c("5mC", "5hmC")) {
    tt   <- if (mod_label == "5mC") tt_5mC else tt_5hmC
    hits <- ANK1_PROBES[ANK1_PROBES %in% rownames(tt)]
    cat(sprintf("\n%s — ANK1 probes:\n", mod_label))
    if (length(hits) == 0) { cat("  Neither probe found.\n"); next }
    out      <- tt[hits, c("logFC", "AveExpr", "P.Value", "adj.P.Val", "sig_ew")]
    out$rank <- match(hits, rownames(tt))
    print(out)
  }
  
  # Boxplot: 2 probes x 2 modifications (2x2 facet grid)
  ank1_long <- bind_rows(lapply(ANK1_PROBES, function(probe) {
    rows <- list()
    if (probe %in% rownames(mC_mat))
      rows[["5mC"]]  <- data.frame(Probe = probe, Modification = "5mC",
                                   Beta  = as.numeric(mC_mat[ probe, bin_s]),
                                   Braak_grp = braak_bin_f, stringsAsFactors = FALSE)
    if (probe %in% rownames(hmC_mat))
      rows[["5hmC"]] <- data.frame(Probe = probe, Modification = "5hmC",
                                   Beta  = as.numeric(hmC_mat[probe, bin_s]),
                                   Braak_grp = braak_bin_f, stringsAsFactors = FALSE)
    bind_rows(rows)
  }))
  
  if (nrow(ank1_long) == 0) {
    cat("[ANK1]  Neither probe found in proportion matrices — skipping plot.\n")
    return(invisible(NULL))
  }
  
  ank1_long$Probe        <- factor(ank1_long$Probe,        levels = ANK1_PROBES)
  ank1_long$Modification <- factor(ank1_long$Modification, levels = c("5mC", "5hmC"))
  
  p_ank1 <- ggplot(ank1_long, aes(x = Braak_grp, y = Beta, fill = Braak_grp)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.2) +
    facet_grid(Probe ~ Modification, scales = "free_y") +
    scale_fill_manual(values = c("low" = "steelblue", "high" = "tomato")) +
    theme_minimal(base_size = 12) +
    theme(legend.position  = "none",
          strip.text        = element_text(face = "bold"),
          strip.background  = element_rect(fill = "grey95")) +
    labs(title    = paste0("ANK1 validation — cg05066959 and cg11823178",
                           " (", pipeline_label, ")"),
         subtitle = "Expected: 5mC increases, 5hmC decreases in high Braak stage",
         x        = "Braak group (low = 0-II, high = V-VI)",
         y        = "Estimated proportion")
  
  fname <- file.path(base_dir,
                     paste0("results_ANK1_validation_",
                            gsub(" ", "_", pipeline_label), ".jpeg"))
  jpeg(fname, width = 8, height = 7, units = "in", res = 300, quality = 95)
  print(p_ank1)
  dev.off()
  cat(sprintf("[Saved]  %s\n", basename(fname)))
}


# =============================================================================
# SECTION 5. PAPER DESIGN REPLICATION + ANK1 VALIDATION (original pipeline)
# =============================================================================
# Design: M ~ Braak (continuous 0-VI) + age + sex + neuron proportion
# Significance threshold: experiment-wide P < 2.4e-7 (Saffari et al. 2018)
# =============================================================================

cat("\n===== SECTION 5: Paper design replication (original pipeline) =====\n")
cat("Design: M ~ Braak + age + sex + neuron proportion\n")
cat("Significance threshold: experiment-wide P < 2.4e-7\n\n")

cat("-- 5mC (OxBS) --\n")
tt_5mC_paper  <- run_paper_ewas(M_5mC,       braak_mlml, age_mlml, sex_mlml,
                                neuron_prop, "5mC (OxBS)")
cat("\n-- 5hmC (filtered) --\n")
tt_5hmC_paper <- run_paper_ewas(M_5hmC_filt, braak_mlml, age_mlml, sex_mlml,
                                neuron_prop, "5hmC (filtered)")
cat("\n-- uC --\n")
tt_uC_paper   <- run_paper_ewas(M_uC,        braak_mlml, age_mlml, sex_mlml,
                                neuron_prop, "uC")

# Volcano plots — experiment-wide threshold -----------------------------------
volcano_ew <- function(tt, title_str) {
  top20    <- rownames(head(tt[order(tt$P.Value), ], 20))
  tt$label <- ifelse(tt$sig_ew & rownames(tt) %in% top20, rownames(tt), NA)
  ggplot(tt, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(aes(colour = sig_ew), alpha = 0.4, size = 0.6) +
    geom_hline(yintercept = -log10(P_EW), linetype = "dashed", colour = "grey50") +
    ggrepel::geom_text_repel(aes(label = label), size = 2.5,
                             max.overlaps = 30, na.rm = TRUE) +
    scale_colour_manual(values = c("FALSE" = "grey80", "TRUE" = "red"),
                        name   = "Significant") +
    theme_minimal(base_size = 11) +
    labs(title    = title_str,
         subtitle = paste0(sum(tt$sig_ew), " CpGs significant (P < 2.4e-7)"),
         x = "logFC (M-value)", y = "-log10(P)")
}

jpeg(file.path(base_dir, "results_fig1_volcano_paper_replication.jpeg"),
     width = 14, height = 5, units = "in", res = 300, quality = 95)
grid.arrange(
  volcano_ew(tt_5mC_paper,  "5mC (OxBS) — original pipeline"),
  volcano_ew(tt_5hmC_paper, "5hmC (filtered) — original pipeline"),
  volcano_ew(tt_uC_paper,   "uC — original pipeline"),
  ncol = 3)
dev.off()
cat("[Saved]  results_fig1_volcano_paper_replication.jpeg\n")

# Replication report: cg10696062 (primary DMP reported in Smith et al. 2019) --
probe_of_interest <- "cg10696062"

compute_delta <- function(prop_mat, braak_v, probe_id) {
  if (!probe_id %in% rownames(prop_mat)) return(NA_real_)
  vals <- prop_mat[probe_id, ]
  sids <- names(vals)
  b0   <- sids[braak_v[sids] == 0]
  b6   <- sids[braak_v[sids] == 6]
  if (length(b0) == 0 || length(b6) == 0) return(NA_real_)
  (mean(vals[b6], na.rm = TRUE) - mean(vals[b0], na.rm = TRUE)) * 100
}

cat(sprintf("\n===== Replication report: %s (original pipeline) =====\n",
            probe_of_interest))
for (lst in list(list(tt = tt_5mC_paper,  mat = mC_matrix,       label = "5mC (OxBS)"),
                 list(tt = tt_5hmC_paper, mat = hmC_matrix_filt, label = "5hmC (filtered)"),
                 list(tt = tt_uC_paper,   mat = uC_matrix,       label = "uC"))) {
  cat(sprintf("\n--- %s ---\n", lst$label))
  if (!probe_of_interest %in% rownames(lst$tt)) { cat("  Not found\n"); next }
  row   <- lst$tt[probe_of_interest, ]
  delta <- compute_delta(lst$mat, braak_mlml, probe_of_interest)
  cat(sprintf("  P      = %.3e\n  adj.P  = %.3e\n  logFC  = %+.4f\n  sig_ew = %s\n  Delta  = %+.2f%%\n",
              row$P.Value, row$adj.P.Val, row$logFC,
              ifelse(row$sig_ew, "YES ***", "no"), delta))
}

cat("\n--- Experiment-wide significant hits (original pipeline) ---\n")
cat(sprintf("  5mC: %d  |  5hmC: %d  |  uC: %d\n",
            sum(tt_5mC_paper$sig_ew), sum(tt_5hmC_paper$sig_ew),
            sum(tt_uC_paper$sig_ew)))

write.csv(tt_5mC_paper[  tt_5mC_paper$sig_ew,  ], file.path(base_dir, "partA_sig_EW_5mC.csv"))
write.csv(tt_5hmC_paper[ tt_5hmC_paper$sig_ew, ], file.path(base_dir, "partA_sig_EW_5hmC.csv"))
write.csv(tt_uC_paper[   tt_uC_paper$sig_ew,   ], file.path(base_dir, "partA_sig_EW_uC.csv"))
cat("[Saved]  partA_sig_EW_*.csv\n")

# ANK1 validation — original pipeline -----------------------------------------
# Binary Braak split on braak_mlml for the boxplot (matches later Section 9)
braak_bin_tmp     <- ifelse(braak_mlml <= 2, "low",
                            ifelse(braak_mlml >= 5, "high", NA))
braak_bin_fac_tmp <- factor(braak_bin_tmp[!is.na(braak_bin_tmp)], levels = c("low","high"))
bin_sids_tmp      <- names(braak_bin_fac_tmp)

run_ank1_validation(
  tt_5mC       = tt_5mC_paper,
  tt_5hmC      = tt_5hmC_paper,
  mC_mat       = mC_matrix[,  bin_sids_tmp, drop = FALSE],
  hmC_mat      = hmC_matrix_filt[, bin_sids_tmp, drop = FALSE],
  braak_bin_f  = braak_bin_fac_tmp,
  bin_s        = bin_sids_tmp,
  pipeline_label = "original pipeline"
)


# =============================================================================
# SECTION 6. RELOAD DATA
# =============================================================================
# Reloads raw data from disk for the extended QC pipeline (NOOB + dasen).
# Metadata and covariate extraction steps are identical to Section 1.
# =============================================================================

gse   <- getGEO("GSE105109", getGPL = FALSE)[[1]]
pheno <- pData(gse)

roman_to_number <- c("I"=1, "II"=2, "III"=3, "IV"=4, "V"=5, "VI"=6)
pheno_braak_raw <- gsub("braak stage: ", "", pheno$characteristics_ch1.3,
                        ignore.case = TRUE)
pheno_braak_all <- ifelse(pheno_braak_raw == "0", 0,
                          roman_to_number[pheno_braak_raw])
names(pheno_braak_all) <- rownames(pheno)

na_rows    <- which(is.na(pheno_braak_all))
gsm_remove <- rownames(pheno)[na_rows]
pheno      <- pheno[!rownames(pheno) %in% gsm_remove, ]
pheno_braak_all <- pheno_braak_all[!is.na(pheno_braak_all)]

idat_files <- list.files(out_dir, pattern = "\\.idat$", full.names = TRUE)
file.remove(idat_files[grepl(paste(gsm_remove, collapse = "|"), basename(idat_files))])

is_ec         <- grepl("entorhinal cortex", pheno$title, ignore.case = TRUE)
gsm_remove_cb <- rownames(pheno)[!is_ec]
pheno         <- pheno[is_ec, ]
pheno_braak_all <- pheno_braak_all[rownames(pheno)]

idat_files <- list.files(out_dir, pattern = "\\.idat$", full.names = TRUE)
file.remove(idat_files[grepl(paste(gsm_remove_cb, collapse = "|"), basename(idat_files))])
cat("EC samples retained:", nrow(pheno), "\n")

num  <- as.numeric(sub(".*_(\\d+)$", "\\1", pheno$title))
num3 <- sprintf("%03d", num)
pheno$sample_ID <- paste0("sample", num3, "_EC")
pheno$assay <- ifelse(grepl("_oxbs_", pheno$title, ignore.case = TRUE), "OxBS",
                      ifelse(grepl("_bs_", pheno$title, ignore.case = TRUE), "BS", NA))

conversion_table <- pheno[, c("geo_accession", "sample_ID", "assay")]
gsm_to_sample    <- setNames(conversion_table$sample_ID,  conversion_table$geo_accession)
sample_to_gsm    <- setNames(conversion_table$geo_accession, conversion_table$sample_ID)

pheno_BS   <- pheno[pheno$assay == "BS",   ]
pheno_OxBS <- pheno[pheno$assay == "OxBS", ]

RGset_all    <- read.metharray.exp(base = out_dir, extended = TRUE)
rg_gsm_clean <- sub("_(.*)", "", colnames(RGset_all))
colnames(RGset_all) <- rg_gsm_clean

RGset_BS   <- RGset_all[, colnames(RGset_all) %in% rownames(pheno_BS)]
RGset_OxBS <- RGset_all[, colnames(RGset_all) %in% rownames(pheno_OxBS)]
cat("RGset_BS samples:",   ncol(RGset_BS),   "| expected ~95\n")
cat("RGset_OxBS samples:", ncol(RGset_OxBS), "| expected ~95\n")

sex_BS   <- factor(tolower(extract_covar(pheno_BS,   "characteristics_ch1.1", "gender: ")))
sex_OxBS <- factor(tolower(extract_covar(pheno_OxBS, "characteristics_ch1.1", "gender: ")))
names(sex_BS)   <- rownames(pheno_BS)
names(sex_OxBS) <- rownames(pheno_OxBS)

age_BS   <- extract_covar(pheno_BS,   "characteristics_ch1.2", "age at death: *", type = "numeric")
age_OxBS <- extract_covar(pheno_OxBS, "characteristics_ch1.2", "age at death: *", type = "numeric")
names(age_BS)   <- rownames(pheno_BS)
names(age_OxBS) <- rownames(pheno_OxBS)

pheno_braak_BS   <- pheno_braak_all[rownames(pheno_BS)]
pheno_braak_OxBS <- pheno_braak_all[rownames(pheno_OxBS)]

cr_probes    <- trimws(as.character(read.csv(cr_file,     header = TRUE,  stringsAsFactors = FALSE)[, 1]))
multi_probes <- trimws(as.character(read.table(multi_file, header = FALSE, stringsAsFactors = FALSE)[, 1]))
anno     <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_snp <- anno[, c("CpG_maf", "SBE_maf")]
cat(sprintf("[Blacklists]  Chen XR: %d | BOWTIE2: %d\n",
            length(cr_probes), length(multi_probes)))


# =============================================================================
# SECTION 7. QC — EXTENDED PIPELINE (NOOB + dasen)
# =============================================================================
# Extends the original QC by adding preprocessNoob background correction
# before dasen normalisation. All other steps are identical to Section 2.
# =============================================================================

run_QC_ext <- function(RGset, sex_vec, age_vec, braak_vec, pheno_df, assay_label) {
  
  qc_log   <- list()
  log_step <- function(step, probes, samples)
    qc_log[[length(qc_log) + 1]] <<- data.frame(
      Step = step, Probes = probes, Samples = samples, stringsAsFactors = FALSE)
  log_step("Input", nrow(RGset), ncol(RGset))
  
  # ── (1) Sex concordance check ───────────────────────────────────────────────
  sex_predicted  <- getSex(mapToGenome(RGset))
  sex_pred_clean <- ifelse(sex_predicted$predictedSex == "M", "m", "f")
  sex_reported   <- as.character(sex_vec[colnames(RGset)])
  mismatch       <- sex_pred_clean != sex_reported
  cat(sprintf("[Sex check]  Mismatches: %d\n", sum(mismatch)))
  
  sex_df <- data.frame(xMed = sex_predicted$xMed, yMed = sex_predicted$yMed,
                       Reported = sex_reported, Mismatch = mismatch)
  p_sex <- ggplot(sex_df, aes(x = xMed, y = yMed, colour = Reported, shape = Mismatch)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_colour_manual(values = c("m" = "steelblue", "f" = "tomato"),
                        labels = c("m" = "Reported male", "f" = "Reported female")) +
    scale_shape_manual(values = c("FALSE" = 19, "TRUE" = 4),
                       labels = c("FALSE" = "Pass", "TRUE" = "Mismatch")) +
    theme_minimal(base_size = 11) +
    labs(title    = paste0("Sex concordance check — ", assay_label,
                           " (extended pipeline)"),
         subtitle = sprintf("%d mismatch(es) removed", sum(mismatch)),
         x        = "X chromosome median log2 intensity",
         y        = "Y chromosome median log2 intensity",
         colour   = "Reported sex", shape = "QC status")
  ggsave(file.path(base_dir, paste0("QC_ext_sex_check_", assay_label, ".jpeg")),
         p_sex, width = 7, height = 5, dpi = 300, device = "jpeg")
  
  RGset     <- RGset[,    !mismatch]
  sex_vec   <- sex_vec[   colnames(RGset)]
  age_vec   <- age_vec[   colnames(RGset)]
  braak_vec <- braak_vec[ colnames(RGset)]
  pheno_df  <- pheno_df[  colnames(RGset), ]
  log_step("After sex check", nrow(RGset), ncol(RGset))
  
  cat("[Pre-extract]  Computing detection p-values and beadcounts ...\n")
  det_p  <- detectionP(RGset)
  bc_mat <- beadcount(RGset)
  
  # ── (2) Blacklist filtering ─────────────────────────────────────────────────
  MSet_raw <- preprocessRaw(RGset)
  log_step("preprocessRaw", nrow(MSet_raw), ncol(MSet_raw))
  
  MSet_raw <- MSet_raw[!rownames(MSet_raw) %in% cr_probes, ]
  log_step("After Chen XR blacklist", nrow(MSet_raw), ncol(MSet_raw))
  MSet_raw <- MSet_raw[!rownames(MSet_raw) %in% multi_probes, ]
  log_step("After BOWTIE2 blacklist", nrow(MSet_raw), ncol(MSet_raw))
  
  anno_sub   <- anno_snp[rownames(MSet_raw), ]
  snp_probes <- rownames(anno_sub)[which(anno_sub$CpG_maf > 0.05 |
                                           anno_sub$SBE_maf > 0.05)]
  MSet_raw   <- MSet_raw[!rownames(MSet_raw) %in% snp_probes, ]
  cat(sprintf("[Pre-NOOB]  Probes: %d | Samples: %d\n",
              nrow(MSet_raw), ncol(MSet_raw)))
  log_step("After SNP filter (MAF > 5%)", nrow(MSet_raw), ncol(MSet_raw))
  
  # ── (3) NOOB background correction ─────────────────────────────────────────
  cat("[NOOB]  Running preprocessNoob ...\n")
  MSet_noob <- preprocessNoob(RGset)
  MSet_noob <- MSet_noob[rownames(MSet_raw), ]
  log_step("After NOOB", nrow(MSet_noob), ncol(MSet_noob))
  
  # ── (4) dasen normalisation ─────────────────────────────────────────────────
  cat("[dasen]  Running normalisation ...\n")
  b_noob  <- getBeta(MSet_noob)
  df_noob <- data.frame(
    beta  = as.vector(b_noob[sample(nrow(b_noob), min(5000, nrow(b_noob))), ]),
    stage = "NOOB")
  
  MSet_norm <- dasen(MSet_noob)
  b_values  <- getBeta(MSet_norm)
  M_values  <- getM(MSet_norm)
  
  df_dasen <- data.frame(
    beta  = as.vector(b_values[sample(nrow(b_values), min(5000, nrow(b_values))), ]),
    stage = "NOOB + dasen")
  
  p_dist <- ggplot(rbind(df_noob, df_dasen), aes(x = beta, colour = stage)) +
    geom_density(linewidth = 0.8, adjust = 1.2) +
    scale_colour_manual(values = c("NOOB" = "darkorange", "NOOB + dasen" = "steelblue")) +
    theme_minimal(base_size = 11) +
    labs(title    = paste0("Beta-value distribution — ", assay_label,
                           " (extended pipeline)"),
         subtitle = "Random sample of 5,000 probes per normalisation stage",
         x = "Beta value", y = "Density", colour = "Normalisation stage")
  ggsave(file.path(base_dir, paste0("QC_ext_beta_dist_", assay_label, ".jpeg")),
         p_dist, width = 7, height = 4, dpi = 300, device = "jpeg")
  log_step("After dasen normalisation", nrow(b_values), ncol(b_values))
  
  # ── (5) pfilter ─────────────────────────────────────────────────────────────
  cat("[pfilter]  Aligning QC matrices ...\n")
  common_probes  <- Reduce(intersect, list(rownames(b_values), rownames(det_p), rownames(bc_mat)))
  common_samples <- Reduce(intersect, list(colnames(b_values), colnames(det_p), colnames(bc_mat)))
  
  pf_result <- pfilter(
    mn = b_values[common_probes, common_samples],
    bn = b_values[common_probes, common_samples],
    pn = det_p[   common_probes, common_samples],
    bc = bc_mat[  common_probes, common_samples],
    perc = 5, pthresh = 1, perCount = 5, pnthresh = 0.05, logical.return = FALSE)
  
  b_values <- as.matrix(pf_result$bn)
  cat(sprintf("[pfilter]  Kept %d probes | %d samples\n",
              nrow(b_values), ncol(b_values)))
  log_step("After pfilter", nrow(b_values), ncol(b_values))
  
  b_clip   <- pmax(pmin(b_values, 1 - 1e-6), 1e-6)
  M_values <- log2(b_clip / (1 - b_clip))
  
  surviving <- colnames(b_values)
  sex_vec   <- sex_vec[   surviving]
  age_vec   <- age_vec[   surviving]
  braak_vec <- braak_vec[ surviving]
  pheno_df  <- pheno_df[  surviving, ]
  
  cat(sprintf("\n[%s FINAL]  Samples: %d | Probes: %d\n",
              assay_label, ncol(b_values), nrow(b_values)))
  cat("           Paper target: ~91 BS | ~367,480 probes\n")
  
  qc_table <- do.call(rbind, qc_log)
  cat("\n[QC tracking table —", assay_label, "(extended pipeline)]\n")
  print(qc_table)
  
  list(b_values = b_values, M_values = M_values, pheno = pheno_df,
       sex = sex_vec, age = age_vec, braak = braak_vec, qc_table = qc_table)
}

# Run extended QC on BS data --------------------------------------------------
cat("Running QC on BS data (extended pipeline)...\n")
qc_BS           <- run_QC_ext(RGset_BS, sex_BS, age_BS, pheno_braak_BS, pheno_BS, "BS")
b_values_BS     <- qc_BS$b_values;  M_values_BS    <- qc_BS$M_values
pheno_BS        <- qc_BS$pheno;     sex_BS         <- qc_BS$sex
age_BS          <- qc_BS$age;       pheno_braak_BS <- qc_BS$braak
qc_table_BS_ext <- qc_BS$qc_table

# Run extended QC on OxBS data ------------------------------------------------
cat("\nRunning QC on OxBS data (extended pipeline)...\n")
qc_OxBS           <- run_QC_ext(RGset_OxBS, sex_OxBS, age_OxBS, pheno_braak_OxBS, pheno_OxBS, "OxBS")
b_values_OxBS     <- qc_OxBS$b_values;  M_values_OxBS    <- qc_OxBS$M_values
pheno_OxBS        <- qc_OxBS$pheno;     sex_OxBS         <- qc_OxBS$sex
age_OxBS          <- qc_OxBS$age;       pheno_braak_OxBS <- qc_OxBS$braak
qc_table_OxBS_ext <- qc_OxBS$qc_table

cat(sprintf("\n[QC Summary — extended pipeline]\n"))
cat(sprintf("  BS:   %d samples | %d probes\n", ncol(b_values_BS),   nrow(b_values_BS)))
cat(sprintf("  OxBS: %d samples | %d probes\n", ncol(b_values_OxBS), nrow(b_values_OxBS)))

# Paired matching and probe intersection --------------------------------------
colnames(b_values_BS)   <- gsm_to_sample[colnames(b_values_BS)]
colnames(M_values_BS)   <- gsm_to_sample[colnames(M_values_BS)]
rownames(pheno_BS)      <- pheno_BS$sample_ID
names(sex_BS)           <- gsm_to_sample[names(sex_BS)]
names(age_BS)           <- gsm_to_sample[names(age_BS)]
names(pheno_braak_BS)   <- gsm_to_sample[names(pheno_braak_BS)]

colnames(b_values_OxBS) <- gsm_to_sample[colnames(b_values_OxBS)]
colnames(M_values_OxBS) <- gsm_to_sample[colnames(M_values_OxBS)]
rownames(pheno_OxBS)    <- pheno_OxBS$sample_ID
names(sex_OxBS)         <- gsm_to_sample[names(sex_OxBS)]
names(age_OxBS)         <- gsm_to_sample[names(age_OxBS)]
names(pheno_braak_OxBS) <- gsm_to_sample[names(pheno_braak_OxBS)]

matched_samples <- intersect(colnames(b_values_BS), colnames(b_values_OxBS))
cat(sprintf("[Pairing]  Matched donors: %d\n", length(matched_samples)))

for (obj in c("b_values_BS", "M_values_BS", "b_values_OxBS", "M_values_OxBS")) {
  m <- get(obj); assign(obj, m[, sort(matched_samples)])
}
pheno_BS         <- pheno_BS[        sort(matched_samples), ]
sex_BS           <- sex_BS[          sort(matched_samples)]
age_BS           <- age_BS[          sort(matched_samples)]
pheno_braak_BS   <- pheno_braak_BS[  sort(matched_samples)]
pheno_OxBS       <- pheno_OxBS[      sort(matched_samples), ]
sex_OxBS         <- sex_OxBS[        sort(matched_samples)]
age_OxBS         <- age_OxBS[        sort(matched_samples)]
pheno_braak_OxBS <- pheno_braak_OxBS[sort(matched_samples)]

matched_probes <- intersect(rownames(b_values_BS), rownames(b_values_OxBS))
cat(sprintf("[Probes]   Shared probes: %d\n", length(matched_probes)))

b_values_BS   <- b_values_BS[  sort(matched_probes), ]
M_values_BS   <- M_values_BS[  sort(matched_probes), ]
b_values_OxBS <- b_values_OxBS[sort(matched_probes), ]
M_values_OxBS <- M_values_OxBS[sort(matched_probes), ]

stopifnot(identical(colnames(b_values_BS), colnames(b_values_OxBS)))
stopifnot(identical(rownames(b_values_BS), rownames(b_values_OxBS)))
cat("[Sanity check]  BS and OxBS matrices aligned.\n")

# Subtraction decomposition ---------------------------------------------------
b_values_BS_sub   <- as.matrix(b_values_BS[,   sort(matched_samples)])
b_values_OxBS_sub <- as.matrix(b_values_OxBS[, sort(matched_samples)])

BETA_THRESH      <- 0.02
keep_subtraction <- rowMeans(b_values_BS_sub,   na.rm = TRUE) > BETA_THRESH &
  rowMeans(b_values_OxBS_sub, na.rm = TRUE) > BETA_THRESH
cat(sprintf("[Beta filter %.2f]  Probes retained: %d\n",
            BETA_THRESH, sum(keep_subtraction)))

b_final_BS   <- b_values_BS_sub[  keep_subtraction, ]
b_final_OxBS <- b_values_OxBS_sub[keep_subtraction, ]

mC_matrix  <- b_final_OxBS
hmC_matrix <- b_final_BS - b_final_OxBS
uC_matrix  <- 1 - b_final_BS

M_5mC   <- prop_to_M(mC_matrix)
M_5hmC  <- prop_to_M(hmC_matrix)
M_uC    <- prop_to_M(uC_matrix)
M_total <- prop_to_M(b_values_BS)

hmC_detected    <- rowMeans(hmC_matrix > 0, na.rm = TRUE) > 0.5
hmC_matrix_filt <- hmC_matrix[hmC_detected, ]
M_5hmC_filt     <- M_5hmC[    hmC_detected, ]
cat(sprintf("[5hmC filter]  Probes present in >50%% of samples: %d\n",
            sum(hmC_detected)))

# Covariate alignment ---------------------------------------------------------
braak_mlml <- pheno_braak_BS[matched_samples]
age_mlml   <- age_BS[        matched_samples]
sex_mlml   <- sex_BS[        matched_samples]

cat("braak_mlml — range:", paste(range(braak_mlml, na.rm = TRUE), collapse = " to "),
    "| NAs:", sum(is.na(braak_mlml)), "\n")
cat("age_mlml   — range:", paste(range(age_mlml,   na.rm = TRUE), collapse = " to "),
    "| NAs:", sum(is.na(age_mlml)), "\n")
cat("sex_mlml   — table:", paste(names(table(sex_mlml)), table(sex_mlml),
                                 sep = "=", collapse = ", "), "\n")


# =============================================================================
# SECTION 8. PAPER RESULT VERIFICATION + ANK1 VALIDATION (extended pipeline)
# =============================================================================

cat("\n===== SECTION 8: Paper result verification (extended pipeline) =====\n")
cat("Design: M ~ Braak + age + sex + neuron proportion\n")
cat("Significance threshold: experiment-wide P < 2.4e-7\n\n")

cat("-- 5mC (OxBS) --\n")
tt_5mC_paper  <- run_paper_ewas(M_5mC,       braak_mlml, age_mlml, sex_mlml,
                                neuron_prop, "5mC (OxBS)")
cat("\n-- 5hmC (filtered) --\n")
tt_5hmC_paper <- run_paper_ewas(M_5hmC_filt, braak_mlml, age_mlml, sex_mlml,
                                neuron_prop, "5hmC (filtered)")
cat("\n-- uC --\n")
tt_uC_paper   <- run_paper_ewas(M_uC,        braak_mlml, age_mlml, sex_mlml,
                                neuron_prop, "uC")

cat(sprintf("\n===== Replication report: %s (extended pipeline) =====\n",
            probe_of_interest))
for (lst in list(list(tt = tt_5mC_paper,  mat = mC_matrix,       label = "5mC (OxBS)"),
                 list(tt = tt_5hmC_paper, mat = hmC_matrix_filt, label = "5hmC (filtered)"),
                 list(tt = tt_uC_paper,   mat = uC_matrix,       label = "uC"))) {
  cat(sprintf("\n--- %s ---\n", lst$label))
  if (!probe_of_interest %in% rownames(lst$tt)) { cat("  Not found\n"); next }
  row   <- lst$tt[probe_of_interest, ]
  delta <- compute_delta(lst$mat, braak_mlml, probe_of_interest)
  cat(sprintf("  P      = %.3e\n  adj.P  = %.3e\n  logFC  = %+.4f\n  sig_ew = %s\n  Delta  = %+.2f%%\n",
              row$P.Value, row$adj.P.Val, row$logFC,
              ifelse(row$sig_ew, "YES ***", "no"), delta))
}

cat("\n--- Experiment-wide significant hits (extended pipeline) ---\n")
cat(sprintf("  5mC: %d  |  5hmC: %d  |  uC: %d\n",
            sum(tt_5mC_paper$sig_ew), sum(tt_5hmC_paper$sig_ew),
            sum(tt_uC_paper$sig_ew)))

# ANK1 validation — extended pipeline -----------------------------------------
# Binary Braak split is defined later in Section 9 but needed here for the
# boxplot; pre-computed on braak_mlml from the extended pipeline.
braak_bin_vec_tmp <- ifelse(braak_mlml <= 2, "low",
                            ifelse(braak_mlml >= 5, "high", NA))
braak_bin_fac_s8  <- factor(braak_bin_vec_tmp[!is.na(braak_bin_vec_tmp)],
                            levels = c("low", "high"))
bin_sids_s8       <- names(braak_bin_fac_s8)

run_ank1_validation(
  tt_5mC       = tt_5mC_paper,
  tt_5hmC      = tt_5hmC_paper,
  mC_mat       = mC_matrix[,  bin_sids_s8, drop = FALSE],
  hmC_mat      = hmC_matrix_filt[, bin_sids_s8, drop = FALSE],
  braak_bin_f  = braak_bin_fac_s8,
  bin_s        = bin_sids_s8,
  pipeline_label = "extended pipeline"
)


# =============================================================================
# SECTION 9. EXTENDED ANALYSIS
# =============================================================================
# Design: M ~ Braak_binary + age + sex + cell-type proportions (4 types)
# Binary Braak: low (0-II) vs high (V-VI); Braak III-IV excluded
# Five deconvolution methods compared; no-deconvolution as baseline
# Significance: BH adj.P < 0.05 AND |logFC| > 0.20
# =============================================================================

cat("\n===== SECTION 9: Extended analysis =====\n")

braak_bin_vec <- ifelse(braak_mlml <= 2, "low",
                        ifelse(braak_mlml >= 5, "high", NA))
keep_bin      <- !is.na(braak_bin_vec)
braak_bin_fac <- factor(braak_bin_vec[keep_bin], levels = c("low", "high"))
bin_sids      <- names(braak_bin_fac)
age_bin       <- age_mlml[bin_sids]
sex_bin       <- sex_mlml[bin_sids]

cat(sprintf("Binary Braak — low (0-II): %d | high (V-VI): %d\n",
            sum(braak_bin_fac == "low"), sum(braak_bin_fac == "high")))

M_5mC_bin  <- M_5mC[,        bin_sids]
M_5hmC_bin <- M_5hmC_filt[,  bin_sids]
M_uC_bin   <- M_uC[,         bin_sids]

cat("\n-- 5mC (binary Braak) --\n")
res_5mC  <- run_all_designs(M_5mC_bin,  braak_bin_fac, age_bin, sex_bin,
                            cells_nnls, cells_svr, cells_qp, cells_rpc, cells_cp,
                            coef_pattern = "high", label_prefix = "5mC")
cat("\n-- 5hmC (binary Braak) --\n")
res_5hmC <- run_all_designs(M_5hmC_bin, braak_bin_fac, age_bin, sex_bin,
                            cells_nnls, cells_svr, cells_qp, cells_rpc, cells_cp,
                            coef_pattern = "high", label_prefix = "5hmC")


# =============================================================================
# SECTION 10. STATISTICAL OBJECTIVE: P-VALUE DISTRIBUTION
# =============================================================================

cat("\n===== SECTION 10: P-value distribution comparison =====\n")

pval_df_5mC  <- build_pval_df(res_5mC)
pval_df_5hmC <- build_pval_df(res_5hmC)

# Raw p-value histograms ------------------------------------------------------
plot_pval_histogram <- function(pval_df, title_str) {
  n_per_method <- nrow(pval_df) / length(unique(pval_df$Method))
  ggplot(pval_df, aes(x = rawP, fill = Method, colour = Method)) +
    geom_histogram(bins = 50, alpha = 0.3, position = "identity", linewidth = 0.2) +
    geom_hline(yintercept = n_per_method / 50, linetype = "dashed",
               colour = "grey40", linewidth = 0.6) +
    scale_fill_manual(  values = method_colours) +
    scale_colour_manual(values = method_colours) +
    scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1)) +
    facet_wrap(~ Method, ncol = 2, scales = "free_y") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank()) +
    labs(title    = title_str,
         subtitle = "Dashed line = expected count under uniform null",
         x = "Raw p-value", y = "Count")
}

jpeg(file.path(base_dir, "results_fig3_pvalue_histograms_5mC.jpeg"),
     width = 10, height = 8, units = "in", res = 300, quality = 95)
print(plot_pval_histogram(pval_df_5mC, "Raw p-value histograms by design — 5mC"))
dev.off()

jpeg(file.path(base_dir, "results_fig4_pvalue_histograms_5hmC.jpeg"),
     width = 10, height = 8, units = "in", res = 300, quality = 95)
print(plot_pval_histogram(pval_df_5hmC, "Raw p-value histograms by design — 5hmC"))
dev.off()
cat("[Saved]  results_fig3/4_pvalue_histograms JPEGs\n")

# Adjusted p-value ECDF -------------------------------------------------------
plot_pval_ecdf <- function(pval_df, title_str) {
  ggplot(pval_df, aes(x = adjP, colour = Method, linetype = Method)) +
    stat_ecdf(linewidth = 0.8) +
    geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey50",
               linewidth = 0.5) +
    scale_colour_manual(values = method_colours) +
    scale_linetype_manual(values = c("solid","dashed","dotdash",
                                     "longdash","dotted","twodash")) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank()) +
    labs(title    = title_str,
         subtitle = "Vertical dashed line at adj.P = 0.05",
         x = "BH-adjusted p-value", y = "Cumulative proportion",
         colour = "Design", linetype = "Design")
}

jpeg(file.path(base_dir, "results_fig5_ecdf_5mC_5hmC.jpeg"),
     width = 12, height = 5, units = "in", res = 300, quality = 95)
grid.arrange(
  plot_pval_ecdf(pval_df_5mC,  "ECDF of adjusted p-values — 5mC"),
  plot_pval_ecdf(pval_df_5hmC, "ECDF of adjusted p-values — 5hmC"),
  ncol = 2)
dev.off()
cat("[Saved]  results_fig5_ecdf_5mC_5hmC.jpeg\n")

# KS tests vs no-deconvolution baseline ---------------------------------------
run_ks_tests <- function(res_list, mod_label) {
  shared <- Reduce(intersect, lapply(res_list, rownames))
  base_p <- res_list$no[shared, "adj.P.Val"]
  cat("\n", mod_label, "— KS test vs no-deconvolution baseline:\n")
  bind_rows(lapply(c("nnls","svr","qp","rpc","cp"), function(nm) {
    ks <- ks.test(base_p, res_list[[nm]][shared, "adj.P.Val"])
    cat(sprintf("  %-4s  D = %.4f  p = %.2e  significant = %s\n",
                toupper(nm), ks$statistic, ks$p.value,
                ifelse(ks$p.value < 0.05, "YES", "no")))
    data.frame(Method   = toupper(nm),
               D_stat   = round(ks$statistic, 4),
               p_value  = signif(ks$p.value, 3),
               Sig_diff = ks$p.value < 0.05)
  }))
}

ks_5mC  <- run_ks_tests(res_5mC,  "5mC")
ks_5hmC <- run_ks_tests(res_5hmC, "5hmC")

# Proportion of probes with adj.P < 0.05 --------------------------------------
prop_below_05 <- function(res_list, mod_label) {
  cat("\n", mod_label, "— proportion of probes with adj.P < 0.05:\n")
  bind_rows(lapply(names(method_labels), function(nm) {
    tt   <- res_list[[nm]]
    n    <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
    prop <- mean(tt$adj.P.Val < 0.05, na.rm = TRUE)
    cat(sprintf("  %-22s  n = %6d  (%.2f%%)\n", method_labels[nm], n, prop * 100))
    data.frame(Method = method_labels[nm], N = n, Pct = round(prop * 100, 2))
  }))
}

prop_5mC  <- prop_below_05(res_5mC,  "5mC")
prop_5hmC <- prop_below_05(res_5hmC, "5hmC")

cat("\n5mC  Pearson r of logFC vs no-deconvolution:\n")
print(pearson_r_table(res_5mC))
cat("\n5hmC Pearson r of logFC vs no-deconvolution:\n")
print(pearson_r_table(res_5hmC))

# Effect-size scatter plots ---------------------------------------------------
for (mod in list(list(res = res_5mC,  name = "5mC",  fig = "fig6"),
                 list(res = res_5hmC, name = "5hmC", fig = "fig7"))) {
  jpeg(file.path(base_dir,
                 paste0("results_", mod$fig, "_effectsize_", mod$name, ".jpeg")),
       width = 12, height = 8, units = "in", res = 300, quality = 95)
  grid.arrange(
    eff_scatter(mod$res$no, mod$res$nnls, "No deconvolution", "NNLS", "NNLS vs no deconvolution"),
    eff_scatter(mod$res$no, mod$res$svr,  "No deconvolution", "SVR",  "SVR vs no deconvolution"),
    eff_scatter(mod$res$no, mod$res$qp,   "No deconvolution", "QP",   "QP vs no deconvolution"),
    eff_scatter(mod$res$no, mod$res$rpc,  "No deconvolution", "RPC",  "RPC vs no deconvolution"),
    eff_scatter(mod$res$no, mod$res$cp,   "No deconvolution", "CP",   "CP vs no deconvolution"),
    ncol = 3,
    top  = grid::textGrob(
      paste0("Effect-size shift vs no deconvolution — ", mod$name),
      gp = grid::gpar(fontsize = 12, fontface = "bold")))
  dev.off()
  cat(sprintf("[Saved]  results_%s_effectsize_%s.jpeg\n", mod$fig, mod$name))
}


# =============================================================================
# SECTION 11. SENSITIVITY — 5mC vs 5hmC
# =============================================================================

cat("\n===== SECTION 11: Sensitivity analysis — 5mC vs 5hmC =====\n")

sens_df <- data.frame(
  Method    = c("NNLS", "SVR", "QP", "RPC", "CP"),
  KS_D_5mC  = ks_5mC$D_stat,
  KS_D_5hmC = ks_5hmC$D_stat,
  delta_KS  = ks_5mC$D_stat - ks_5hmC$D_stat)
cat("\nSensitivity table:\n"); print(sens_df)

sens_long <- tidyr::pivot_longer(
  dplyr::select(sens_df, Method, KS_D_5mC, KS_D_5hmC),
  cols = c(KS_D_5mC, KS_D_5hmC),
  names_to = "Modification", values_to = "KS_D") %>%
  dplyr::mutate(Modification = dplyr::recode(Modification,
                                             KS_D_5mC = "5mC", KS_D_5hmC = "5hmC"))

jpeg(file.path(base_dir, "results_fig8_KS_sensitivity.jpeg"),
     width = 7, height = 5, units = "in", res = 300, quality = 95)
print(
  ggplot(sens_long, aes(x = Method, y = KS_D, fill = Modification)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = round(KS_D, 3)),
              position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
    scale_fill_manual(values = c("5mC" = "steelblue", "5hmC" = "gold3")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank()) +
    labs(title    = "KS D-statistic vs no-deconvolution baseline — 5mC vs 5hmC",
         subtitle = "Higher D-statistic indicates greater distributional shift after deconvolution",
         x = "Deconvolution method", y = "KS D-statistic", fill = "Modification")
)
dev.off()
cat("[Saved]  results_fig8_KS_sensitivity.jpeg\n")


# =============================================================================
# SECTION 12. BIOLOGICAL OBJECTIVE: PROBE DISCOVERY
# =============================================================================

cat("\n===== SECTION 12: Biological probe discovery =====\n")

anno_450k <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
anno_cols <- c("chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group",
               "Relation_to_Island")

# 5hmC probe classification ---------------------------------------------------
summ_5hmC <- summarise_str(res_5hmC, "5hmC")
cat("\n5hmC stringent hits (adj.P <", P_STR, ", |logFC| >", LFC_STR, "):\n")
print(summ_5hmC)

hits_by_method_hmC     <- lapply(names(res_5hmC), function(nm)
  rownames(res_5hmC[[nm]])[res_5hmC[[nm]]$sig_str])
names(hits_by_method_hmC) <- names(res_5hmC)

hits_no_hmC         <- hits_by_method_hmC$no
all_deconv_hits_hmC <- unique(unlist(hits_by_method_hmC[c("nnls","svr","qp","rpc","cp")]))
artifacts_hmC <- hits_no_hmC[!hits_no_hmC %in% all_deconv_hits_hmC]
robust_hmC    <- hits_no_hmC[ hits_no_hmC %in% all_deconv_hits_hmC]
unmasked_hmC  <- all_deconv_hits_hmC[!all_deconv_hits_hmC %in% hits_no_hmC]

cat("\n--- 5hmC probe classification ---\n")
cat("Artifacts (disappear after deconvolution correction):  ", length(artifacts_hmC), "\n")
if (length(artifacts_hmC) > 0) print(artifacts_hmC)
cat("Robust (consistent across all designs):                ", length(robust_hmC),    "\n")
if (length(robust_hmC)    > 0) print(robust_hmC)
cat("Unmasked (appear only after deconvolution correction): ", length(unmasked_hmC),  "\n")
if (length(unmasked_hmC)  > 0) print(unmasked_hmC)

# 5mC probe classification ----------------------------------------------------
summ_5mC <- summarise_str(res_5mC, "5mC")
cat("\n5mC stringent hits (adj.P <", P_STR, ", |logFC| >", LFC_STR, "):\n")
print(summ_5mC)

hits_by_method_mC     <- lapply(names(res_5mC), function(nm)
  rownames(res_5mC[[nm]])[res_5mC[[nm]]$sig_str])
names(hits_by_method_mC) <- names(res_5mC)

hits_no_mC         <- hits_by_method_mC$no
all_deconv_hits_mC <- unique(unlist(hits_by_method_mC[c("nnls","svr","qp","rpc","cp")]))
artifacts_mC <- hits_no_mC[!hits_no_mC %in% all_deconv_hits_mC]
robust_mC    <- hits_no_mC[ hits_no_mC %in% all_deconv_hits_mC]
unmasked_mC  <- all_deconv_hits_mC[!all_deconv_hits_mC %in% hits_no_mC]

cat("\n--- 5mC probe classification ---\n")
cat("Artifacts (disappear after deconvolution correction):  ", length(artifacts_mC), "\n")
if (length(artifacts_mC) > 0) print(artifacts_mC)
cat("Robust (consistent across all designs):                ", length(robust_mC),    "\n")
if (length(robust_mC)    > 0) print(robust_mC)
cat("Unmasked (appear only after deconvolution correction): ", length(unmasked_mC),  "\n")
if (length(unmasked_mC)  > 0) print(unmasked_mC)

# Probe annotation report (5hmC primary discovery framework) ------------------
all_poi   <- unique(c(hits_no_hmC, all_deconv_hits_hmC))
artifacts <- artifacts_hmC
unmasked  <- unmasked_hmC
robust    <- robust_hmC

cat("\n=== Probe annotations (5hmC discovery framework) ===\n")
for (probe in all_poi) {
  cat(sprintf("\n--- %s ---\n", probe))
  ann <- anno_450k[probe, intersect(anno_cols, colnames(anno_450k)), drop = FALSE]
  cat(sprintf("Gene: %s | Region: %s | CpG island: %s | Chr: %s | Pos: %s\n",
              ann$UCSC_RefGene_Name, ann$UCSC_RefGene_Group,
              ann$Relation_to_Island, ann$chr, ann$pos))
  cat("Category:", ifelse(probe %in% artifacts, "ARTIFACT",
                          ifelse(probe %in% unmasked,  "UNMASKED", "ROBUST")), "\n")
  for (nm in names(res_5hmC)) {
    tt <- res_5hmC[[nm]]
    if (!probe %in% rownames(tt)) next
    cat(sprintf("  %-22s  logFC = %+.3f  adj.P = %.2e  significant = %s\n",
                method_labels[nm], tt[probe, "logFC"],
                tt[probe, "adj.P.Val"],
                ifelse(tt[probe, "sig_str"], "YES", "no")))
  }
}

# Volcano plots — 5hmC designs ------------------------------------------------
for (nm in names(res_5hmC)) {
  jpeg(file.path(base_dir, paste0("results_fig9_volcano_5hmC_", toupper(nm), ".jpeg")),
       width = 6, height = 5, units = "in", res = 300, quality = 95)
  print(volcano_plot_str(res_5hmC[[nm]],
                         paste0("Volcano plot — 5hmC, ", method_labels[nm])))
  dev.off()
}
cat("[Saved]  results_fig9_volcano_5hmC_*.jpeg\n")

# Volcano plots — 5mC designs -------------------------------------------------
for (nm in names(res_5mC)) {
  jpeg(file.path(base_dir, paste0("results_fig10_volcano_5mC_", toupper(nm), ".jpeg")),
       width = 6, height = 5, units = "in", res = 300, quality = 95)
  print(volcano_plot_str(res_5mC[[nm]],
                         paste0("Volcano plot — 5mC, ", method_labels[nm])))
  dev.off()
}
cat("[Saved]  results_fig10_volcano_5mC_*.jpeg\n")

# Boxplots of 5hmC proportions for discovered probes --------------------------
if (length(all_poi) > 0) {
  jpeg(file.path(base_dir, "results_fig11_probe_hmC_boxplots.jpeg"),
       width = 7, height = 4 * length(all_poi), units = "in", res = 300, quality = 95)
  for (probe in all_poi) {
    if (!probe %in% rownames(hmC_matrix)) next
    df_p <- data.frame(hmC       = as.numeric(hmC_matrix[probe, bin_sids]),
                       Braak_grp = braak_bin_fac)
    gene      <- anno_450k[probe, "UCSC_RefGene_Name"]
    cat_label <- ifelse(probe %in% artifacts,
                        "Artifact — disappears after deconvolution correction",
                        ifelse(probe %in% unmasked,
                               "Unmasked — appears after deconvolution correction",
                               "Robust — consistent across all designs"))
    print(ggplot(df_p, aes(x = Braak_grp, y = hmC, fill = Braak_grp)) +
            geom_boxplot(alpha = 0.6, outlier.shape = NA) +
            geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
            scale_fill_manual(values = c("low" = "steelblue", "high" = "tomato")) +
            theme_minimal(base_size = 12) +
            theme(legend.position = "none") +
            labs(title    = paste0(probe, ifelse(is.na(gene) | gene == "", "",
                                                 paste0(" (", gene, ")"))),
                 subtitle = cat_label,
                 x = "Braak group (low = 0-II, high = V-VI)",
                 y = "Estimated 5hmC proportion"))
  }
  dev.off()
  cat("[Saved]  results_fig11_probe_hmC_boxplots.jpeg\n")
}


# =============================================================================
# SECTION 13. SAVE ALL RESULTS
# =============================================================================

cat("\n===== SECTION 13: Saving all results =====\n")

# QC tracking tables
write.csv(qc_table_BS,        file.path(base_dir, "QC_tracking_original_BS.csv"),   row.names = FALSE)
write.csv(qc_table_OxBS,      file.path(base_dir, "QC_tracking_original_OxBS.csv"), row.names = FALSE)
write.csv(qc_table_BS_ext,    file.path(base_dir, "QC_tracking_extended_BS.csv"),   row.names = FALSE)
write.csv(qc_table_OxBS_ext,  file.path(base_dir, "QC_tracking_extended_OxBS.csv"), row.names = FALSE)

# Paper replication: experiment-wide significant probes
write.csv(tt_5mC_paper[  tt_5mC_paper$sig_ew,  ], file.path(base_dir, "results_sig_EW_5mC.csv"))
write.csv(tt_5hmC_paper[ tt_5hmC_paper$sig_ew, ], file.path(base_dir, "results_sig_EW_5hmC.csv"))
write.csv(tt_uC_paper[   tt_uC_paper$sig_ew,   ], file.path(base_dir, "results_sig_EW_uC.csv"))

# Extended analysis: stringent hits per design
for (nm in names(res_5hmC)) {
  write.csv(res_5hmC[[nm]][res_5hmC[[nm]]$sig_str, ],
            file.path(base_dir, paste0("results_str_5hmC_", toupper(nm), ".csv")))
  write.csv(res_5mC[[nm]][ res_5mC[[nm]]$sig_str,  ],
            file.path(base_dir, paste0("results_str_5mC_",  toupper(nm), ".csv")))
}

# Deconvolution proportions
write.csv(props_nnls, file.path(base_dir, "deconv_proportions_NNLS.csv"))
write.csv(props_svr,  file.path(base_dir, "deconv_proportions_SVR.csv"))
write.csv(props_qp,   file.path(base_dir, "deconv_proportions_QP.csv"))
write.csv(props_rpc,  file.path(base_dir, "deconv_proportions_RPC.csv"))
write.csv(props_cp,   file.path(base_dir, "deconv_proportions_CP.csv"))

# KS test and sensitivity tables
write.csv(ks_5mC,  file.path(base_dir, "results_KS_5mC.csv"),  row.names = FALSE)
write.csv(ks_5hmC, file.path(base_dir, "results_KS_5hmC.csv"), row.names = FALSE)
write.csv(sens_df, file.path(base_dir, "results_sensitivity_5mC_vs_5hmC.csv"), row.names = FALSE)

# Probe classification
probe_universe <- unique(c(all_poi, artifacts_mC, robust_mC, unmasked_mC))
probe_classification <- data.frame(
  probe        = probe_universe,
  category_hmC = ifelse(probe_universe %in% artifacts_hmC, "Artifact",
                        ifelse(probe_universe %in% unmasked_hmC, "Unmasked", "Robust")),
  category_mC  = ifelse(probe_universe %in% artifacts_mC,  "Artifact",
                        ifelse(probe_universe %in% unmasked_mC,  "Unmasked", "Robust"))
)
write.csv(probe_classification,
          file.path(base_dir, "results_probe_classification.csv"), row.names = FALSE)

cat("\n[DONE]  All results saved", "\n")