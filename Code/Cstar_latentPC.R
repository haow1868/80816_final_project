library(DESeq2)
##filter
setwd("/Users/haowang/Documents/Lab/Collaboration/DNMT/data/DE")
clin_Kyle <- readRDS(file="clin_Kyle.Rdata")
raw <- readRDS(file="swapped_KyleNAc_raw.Rdata") #60676
CPM.Cutoff = 2 #log2CPM>1 so CPM>2
lib.size <- colSums(raw)
MedianLibSize <- median(lib.size)
min.count = ceiling(CPM.Cutoff/1e6*MedianLibSize)
filter = edgeR::filterByExpr(raw, design = model.matrix(~clin_Kyle$group), 
                             min.count = min.count)
sum(filter)
Filtered_raw = raw[filter,]

rownames(clin_Kyle) <- paste0(clin_Kyle$Case,"N")
Filtered_raw <- Filtered_raw[, rownames(clin_Kyle)]
all(colnames(Filtered_raw) == rownames(clin_Kyle)) 
clin_Kyle$BMI <- as.numeric(clin_Kyle$BMI)
clin_Kyle$BMI[is.na(clin_Kyle$BMI)] <- mean(clin_Kyle$BMI, na.rm = TRUE)
names(clin_Kyle)[names(clin_Kyle) == "Age (yr)"] <- "Age"
names(clin_Kyle)[names(clin_Kyle) == "Subject RIN"] <- "RIN"
names(clin_Kyle)[names(clin_Kyle) == "Age (yr)"] <- "Age"
names(clin_Kyle)[names(clin_Kyle) == "Subject RIN"] <- "RIN"

clin_Kyle$group = factor(clin_Kyle$group)
clin_Kyle$Sex = factor(clin_Kyle$Sex)
clin_Kyle$Race = factor(clin_Kyle$Race)
rownames(clin_Kyle) <- paste0("S", clin_Kyle$Case,"N")

Filtered_raw_log = edgeR::cpm(Filtered_raw, log = TRUE, prior.count = 1)



##
## CStaR for MDD bulk RNA-seq (NAc) - latent PC variant
## ======================================================
##
## Uses latent Gaussian copula correlation estimation (Fan et al. 2017,
## adapted from xidongdxi/latentPC) to handle binary MDD response cleanly:
##   - Continuous gene-gene: closed-form Kendall's tau -> Pearson via sin(pi/2 * tau)
##   - Binary-continuous (MDD-gene): numerical inversion of BC bridge function
##   - Latent correlation matrix passed to gaussCItest (and idaFast)
##
## Assumes preprocessing has already been run:
##   - Filtered_raw    : raw counts, genes x samples, post filterByExpr
##   - Filtered_raw_log: log2-CPM, genes x samples
##   - clin_Kyle       : clinical data, aligned to samples
###############################################################################


## Source latentPC functions ==================================================
latentPC_path <- "/Users/haowang/Documents/Lab/Collaboration/DNMT/code/final_project/code"
source(file.path(latentPC_path, "utility.R"))    # tau_fun, BB_fun, BC_fun, CC_fun, etc.
source(file.path(latentPC_path, "latent_pc.R"))  # label_fun, latent_pc
## ============================================================================


## Dependencies ===============================================================
library(pcalg)
library(limma)
library(mvtnorm) 
library(pcaPP)      
library(Matrix)      
## ============================================================================


## Input parameters ===========================================================
alpha       <- 0.1
noStabRuns  <- 100
corCutoff   <- 0.3                
qSet        <- seq(500, 50, -50)
seedVal     <- 42
## ============================================================================


## Verify preprocessing state =================================================
stopifnot(exists("Filtered_raw_log"))
stopifnot(exists("clin_Kyle"))
stopifnot(all(colnames(Filtered_raw_log) == rownames(clin_Kyle)))
cat("Starting CStaR: n =", ncol(Filtered_raw_log),
    ", p =", nrow(Filtered_raw_log), "\n")
## ============================================================================


## Step 1: Define response and covariates =====================================
levels(clin_Kyle$group)
y <- as.numeric(clin_Kyle$group == "MDD")
stopifnot(!anyNA(y))
cat("Response: MDD =", sum(y), ", Control =", sum(1 - y), "\n")

covars <- data.frame(
  Age = as.numeric(clin_Kyle$Age),
  Sex = factor(clin_Kyle$Sex),
  BMI = as.numeric(clin_Kyle$BMI),
  PMI = as.numeric(clin_Kyle[["PMI (hr)"]]),
  pH  = as.numeric(clin_Kyle$pHa),
  RIN = as.numeric(clin_Kyle$RIN)
)

na_counts <- sapply(covars, function(x) sum(is.na(x)))
cat("NAs per covariate:\n"); print(na_counts)
for (cn in names(covars)) {
  x <- covars[[cn]]
  if (is.numeric(x)) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
  } else {
    mode_val <- names(sort(table(x), decreasing = TRUE))[1]
    x[is.na(x)] <- mode_val
  }
  covars[[cn]] <- x
}
stopifnot(!any(is.na(covars)))
## ============================================================================


## Step 2: Residualize expression against covariates ==========================
design_keep <- model.matrix(~ y)
design_cov  <- model.matrix(~ ., data = covars)[, -1]

expr_resid  <- removeBatchEffect(
  x          = as.matrix(Filtered_raw_log),
  covariates = design_cov,
  design     = design_keep
)
cat("Residualization complete.\n")
## ============================================================================


expr_final = expr_resid


## Step 3: Assemble CStaR input matrix ========================================
## CRITICAL for latent PC: do NOT scale the matrix.
##   - MDD must remain integer 0/1 so label_fun() classifies it as "binary".
##   - Genes can stay on residualized log-CPM scale; latent_pc uses Kendall's
##     tau which is invariant to monotone transformations.
gene_mat <- t(expr_final)              # samples x genes
dat <- cbind(MDD = y, gene_mat)        # MDD as 0/1, genes continuous
yInd <- 1
n <- nrow(dat)
p <- ncol(dat) - 1
cat("CStaR input: n =", n, ", p+1 =", ncol(dat), "\n")
## ============================================================================


## Step 4: Stability selection + latent PC + IDA loop =========================
set.seed(seedVal)
subSize <- floor(n / 2)
effMat  <- matrix(0, nrow = noStabRuns, ncol = p,
                  dimnames = list(NULL, colnames(dat)[-yInd]))

t.start <- proc.time()
for (k in 1:noStabRuns) {
  run_start <- Sys.time()
  cat(sprintf("[%s] Stability run %d / %d\n",
              format(run_start, "%H:%M:%S"), k, noStabRuns))
  
  ## Subsample (half-sample, per Meinshausen-Bühlmann)
  subIndices <- sample(n, subSize, replace = FALSE)
  subDat     <- dat[subIndices, ]
  
  ## Marginal correlation screen (point-biserial; applied per subsample)
  marCor  <- apply(subDat, 2, function(x) abs(cor(subDat[, yInd], x)))
  subDat  <- subDat[, marCor >= corCutoff]
  subP    <- ncol(subDat)
  
  if (!"MDD" %in% colnames(subDat)) {
    warning("Run ", k, ": response dropped by correlation screen. Skipping.")
    next
  }
  yInd_sub <- which(colnames(subDat) == "MDD")
  
  ## --- Auto-label variables (MDD as "binary", genes as "continuous") -----
  labels <- tryCatch(
    label_fun(subDat),
    error = function(e) {
      warning("label_fun failed run ", k, ": ", e$message); NULL
    }
  )
  if (is.null(labels)) next
  
  ## Sanity check on first run
  if (k == 1) {
    cat("  Variable types in run 1:\n"); print(table(labels))
  }
  
  ## --- Latent correlation matrix via Kendall's tau bridge functions ------
  sig <- tryCatch(
    latent_pc(subDat, labels),
    error = function(e) {
      warning("latent_pc failed run ", k, ": ", e$message); NULL
    }
  )
  if (is.null(sig)) next
  
  ## --- Project to nearest positive-definite correlation matrix ----------
  ## Required because bridge-inverted matrix is not guaranteed PD.
  sig <- tryCatch(
    as.matrix(Matrix::nearPD(sig, corr = TRUE, maxit = 10000)$mat),
    error = function(e) {
      warning("nearPD failed run ", k, ": ", e$message); NULL
    }
  )
  if (is.null(sig)) next
  
  colnames(sig) <- rownames(sig) <- colnames(subDat)
  
  ## --- PC algorithm with gaussCItest on latent correlation matrix --------
  suffStat <- list(C = sig, n = subSize)
  pcFit <- tryCatch(
    pc(suffStat,
       indepTest   = gaussCItest,
       p           = subP,
       alpha       = alpha,
       skel.method = "stable.fast",   # matches latentPC simulation
       verbose     = FALSE),
    error = function(e) { warning("PC failed run ", k, ": ", e$message); NULL }
  )
  if (is.null(pcFit)) next
  
  ## --- IDA via idaFast on latent correlation matrix ----------------------
  ## NOTE: idaFast accepts a covariance matrix; the latent correlation
  ## matrix serves this role since variables are effectively standardized
  ## in the latent Gaussian copula framework.
  effects_matrix <- suppressWarnings(tryCatch(
    idaFast(yInd_sub,                 # response (target)
            seq_len(subP)[-yInd_sub], # predictors (all genes)
            sig,                      # latent correlation matrix
            pcFit@graph),
    error = function(e) NA
  ))
  
  if (!all(is.na(effects_matrix))) {
    ## For each gene, take min(|effect|) across its multi-set
    for (j in seq_len(subP)) {
      if (j == yInd_sub) next
      gene_name <- colnames(subDat)[j]
      
      ## Map subDat column j to corresponding row in effects_matrix
      row_idx <- if (j < yInd_sub) j else j - 1
      multSet <- effects_matrix[row_idx, ]
      multSet <- multSet[!is.na(multSet)]
      
      if (length(multSet) > 0) {
        effMat[k, gene_name] <- min(abs(multSet))
      }
    }
  }
  
  run_elapsed <- as.numeric(difftime(Sys.time(), run_start, units = "mins"))
  cat(sprintf("  Run %d done in %.1f min. Nonzero effects: %d\n",
              k, run_elapsed, sum(effMat[k, ] > 0)))
}
cat("Stability loop complete. Elapsed:\n"); print(proc.time() - t.start)
## ============================================================================


## Step 5: Aggregate across q grid ============================================
PCER <- function(p, q, freq) {
  out <- q^2 / ((2 * freq - 1) * p^2)
  out[freq <= 0.5] <- 1
  out[out < 0]     <- 1
  out
}

P <- length(qSet)
stabMat  <- matrix(0, nrow = P, ncol = p,
                   dimnames = list(as.character(qSet), colnames(effMat)))
rankMat  <- stabMat
errorMat <- stabMat

for (t.q in seq_along(qSet)) {
  q <- qSet[t.q]
  topMat <- matrix(FALSE, nrow = noStabRuns, ncol = p)
  colnames(topMat) <- colnames(effMat)
  
  for (t.run in 1:noStabRuns) {
    effs <- effMat[t.run, ]
    if (all(effs == 0)) next
    lowerQbound <- sort(effs, decreasing = TRUE)[q]
    if (lowerQbound == 0) {
      lowerQbound <- min(effs[effs > 0])
    }
    topMat[t.run, ] <- effs >= lowerQbound
  }
  
  ## Divide by noStabRuns (NOT n -- corrects original Stekhoven bug)
  stabMat[t.q, ]  <- colSums(topMat) / noStabRuns
  rankMat[t.q, ]  <- rank(-stabMat[t.q, ], ties.method = "min")
  errorMat[t.q, ] <- PCER(p = p, q = q, freq = stabMat[t.q, ])
}

medianRank   <- apply(rankMat,  2, median)
medianEffect <- apply(effMat,   2, median)
medianPCER   <- apply(errorMat, 2, median)

stabRank <- data.frame(
  gene             = names(medianRank),
  medianRank       = medianRank,
  medianEffect     = medianEffect,
  medianPCER       = medianPCER,
  row.names        = NULL,
  stringsAsFactors = FALSE
)
stabRank <- stabRank[order(stabRank$medianRank, -stabRank$medianEffect), ]
## ============================================================================


## Step 6: Save output ========================================================
t.date <- format(Sys.time(), "%Y_%m_%d_%H%M")
out_file <- paste0("CStaR_MDD_ranking_latentPC_", t.date, ".txt")
write.table(stabRank, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Results written to:", out_file, "\n")

saveRDS(list(effMat   = effMat,
             stabMat  = stabMat,
             rankMat  = rankMat,
             errorMat = errorMat,
             stabRank = stabRank,
             params   = list(alpha = alpha, noStabRuns = noStabRuns,
                             corCutoff = corCutoff, qSet = qSet,
                             seed = seedVal,
                             ci_method = "latent_pc")),
        paste0("CStaR_MDD_intermediates_latentPC_", t.date, ".rds"))
## ============================================================================


## Diagnostic previews ========================================================
cat("\n--- Top 20 genes by CStaR median rank ---\n")
print(head(stabRank, 20))

cat("\n--- Genes with median selection frequency > 0.5 at smallest q ---\n")
small_q_row <- stabMat[nrow(stabMat), ]
high_stab   <- names(small_q_row)[small_q_row > 0.5]
cat("Count:", length(high_stab), "\n")
## ============================================================================