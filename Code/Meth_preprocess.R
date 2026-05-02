# -----------------------------------------------------------------------------
# NAc BMIQ + ComBat-Met pipeline
# Input:  swappedmislable_beta_final_filtered.Rdata  (raw, filtered betas, all regions)
# Output: BMIQ-normalized, batch-corrected M-value matrix (NAc only)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# DEFINE HELPER FUNCTIONS
# -----------------------------------------------------------------------------
bmiq_one <- function(beta_col, type_num) {
  keep <- !is.na(beta_col) & !is.na(type_num)
  tab  <- table(type_num[keep])
  if (length(tab) < 2L || any(tab == 0L)) return(beta_col)  # skip if one type missing
  nfit <- min(2000L, as.integer(tab[["1"]]), as.integer(tab[["2"]]))
  out  <- BMIQ(beta.v = beta_col[keep], design.v = type_num[keep],
               nL = 3, nfit = nfit, plots = FALSE)
  nb   <- if (!is.null(out$nbeta)) out$nbeta else out$beta
  res  <- beta_col; res[keep] <- nb; res
}

beta_to_m <- function(b, eps = 1e-6) {
  b <- pmin(pmax(b, eps), 1 - eps)
  log2(b / (1 - b))
}

# -----------------------------------------------------------------------------
# LIBRARIES
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(wateRmelon)
library(ComBatMet)
library(RPMM)

# -----------------------------------------------------------------------------
# 1. LOAD BETA MATRIX AND SUBSET TO NAc
# -----------------------------------------------------------------------------
setwd("~/Documents/Lab/Collaboration/DNMT/data")

beta_final <- readRDS("swappedmislable_beta_final_filtered.Rdata")
beta_final <- as.data.frame(beta_final)

cat("Full matrix:", nrow(beta_final), "probes x", ncol(beta_final), "samples\n")

# Keep only NAc samples (column names ending in "N")
beta_N <- beta_final[, grepl("N$", colnames(beta_final)), drop = FALSE]
rm(beta_final); gc()

cat("NAc subset :", nrow(beta_N), "probes x", ncol(beta_N), "samples\n")

# -----------------------------------------------------------------------------
# 2. LOAD AND ALIGN PROBE-TYPE ANNOTATION
# -----------------------------------------------------------------------------
anno_df <- read.csv("EPIC-8v2-0_A1.csv", skip = 7, stringsAsFactors = FALSE)

anno_df1 <- anno_df %>%
  transmute(
    IlmnID   = Name,
    CHR      = CHR,
    MAPINFO  = MAPINFO,
    Gene     = UCSC_RefGene_Name,
    Group    = UCSC_RefGene_Group,
    Island   = Relation_to_UCSC_CpG_Island,
    Reg      = Regulatory_Feature_Group,
    Type     = Infinium_Design_Type,
    Enhancer = X450k_Enhancer
  )
anno_df1 <- anno_df1[anno_df1$Type %in% c("I", "II"), ]
rm(anno_df); gc()

anno_aln <- anno_df1[match(rownames(beta_N), anno_df1$IlmnID), ]
design.v <- as.numeric(ifelse(anno_aln$Type == "I",  1,
                              ifelse(anno_aln$Type == "II", 2, NA_real_)))

cat("Probes matched:", sum(!is.na(design.v)), "/", nrow(beta_N), "\n")
cat("  Type I :", sum(design.v == 1, na.rm = TRUE), "\n")
cat("  Type II:", sum(design.v == 2, na.rm = TRUE), "\n")

# -----------------------------------------------------------------------------
# 3. BMIQ NORMALIZATION (within sample, probe-type bias)
# -----------------------------------------------------------------------------
cat("\nRunning BMIQ on", ncol(beta_N), "samples...\n")
betas_bmiq <- apply(beta_N, 2, bmiq_one, type_num = design.v)
rownames(betas_bmiq) <- rownames(beta_N)
rm(beta_N); gc()

# -----------------------------------------------------------------------------
# 4. LOAD SENTRIX METADATA AND BUILD BATCH VECTOR
# -----------------------------------------------------------------------------
clin_Kyle    <- readRDS("clin_Kyle.Rdata")
sentrix_meta <- read.csv("sentrix_meta.csv")

x <- as.character(sentrix_meta$sesame_id)
sentrix_meta$Sentrix_ID       <- sub("_.*$", "", x)
sentrix_meta$Sentrix_Position <- sub("^[^_]*_", "", x)

# Restrict metadata to NAc
sentrix_meta_N <- sentrix_meta[sentrix_meta$region %in% c("N"), ]

sm_aln <- sentrix_meta_N[match(colnames(betas_bmiq), sentrix_meta_N$Sample.ID), ]
batch  <- sm_aln$Sentrix_ID

cat("\nSamples mapped to Sentrix metadata:",
    sum(!is.na(batch)), "/", ncol(betas_bmiq), "\n")
cat("Unique Sentrix chips (batches):", length(unique(na.omit(batch))), "\n")

stopifnot(all(!is.na(batch)))   # fail loudly if any sample lacks a batch

# -----------------------------------------------------------------------------
# 5. ComBat-Met BATCH CORRECTION
# -----------------------------------------------------------------------------
mat_in <- as.matrix(betas_bmiq)
rm(betas_bmiq); gc()

beta_cbmet <- ComBat_met(
  vmat   = mat_in,
  dtype  = "b-value",
  batch  = batch,
  ncores = 2                  
)
rm(mat_in); gc()

saveRDS(beta_cbmet, "swappedmislable_bmiq_combatmet_beta_NAc.Rdata")

# -----------------------------------------------------------------------------
# 6. CONVERT TO M-VALUES FOR DOWNSTREAM ANALYSIS
# -----------------------------------------------------------------------------
m_cbmet <- beta_to_m(beta_cbmet)
rm(beta_cbmet); gc()

saveRDS(m_cbmet, "swappedmislable_bmiq_combatmet_m_NAc.Rdata")

cat("\nDone. Final NAc M-value matrix:",
    nrow(m_cbmet), "probes x", ncol(m_cbmet), "samples\n")