## ===========================================================================
## STEP 0: Load CStaR session objects
## ===========================================================================
## These objects were created in the CStaR step (Step 1 of the pipeline).
## We rebuild them here from saved files + the preprocessing already done.

library(limma)
library(edgeR)

setwd("/Users/haowang/Documents/Lab/Collaboration/DNMT/data/DE")

## --- 0.1 Load clinical data and recreate processing -------------------------
clin_Kyle <- readRDS(file = "clin_Kyle.Rdata")

## Re-do the renaming/cleaning from CStaR preprocessing (idempotent if needed)
rownames(clin_Kyle) <- paste0(clin_Kyle$Case, "N")
clin_Kyle$BMI <- as.numeric(clin_Kyle$BMI)
clin_Kyle$BMI[is.na(clin_Kyle$BMI)] <- mean(clin_Kyle$BMI, na.rm = TRUE)
names(clin_Kyle)[names(clin_Kyle) == "Age (yr)"]    <- "Age"
names(clin_Kyle)[names(clin_Kyle) == "Subject RIN"] <- "RIN"
clin_Kyle$group <- factor(clin_Kyle$group)
clin_Kyle$Sex   <- factor(clin_Kyle$Sex)
clin_Kyle$Race  <- factor(clin_Kyle$Race)
rownames(clin_Kyle) <- paste0("S", clin_Kyle$Case,"N")


## --- 0.2 Define binary MDD response -----------------------------------------
levels(clin_Kyle$group)   # confirm coding
y <- as.numeric(clin_Kyle$group == "MDD")  # adjust string if your level differs
stopifnot(!anyNA(y))


## --- 0.3 Build covariate frame (same six as CStaR) --------------------------
covars <- data.frame(
  Age = as.numeric(clin_Kyle$Age),
  Sex = factor(clin_Kyle$Sex),
  BMI = as.numeric(clin_Kyle$BMI),
  PMI = as.numeric(clin_Kyle[["PMI (hr)"]]),
  pH  = as.numeric(clin_Kyle$pHa),
  RIN = as.numeric(clin_Kyle$RIN)
)

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


## --- 0.4 Rebuild residualized expression matrix -----------------------------
## Same preprocessing as CStaR: filterByExpr -> log2-CPM -> covariate residualization
raw <- readRDS(file = "swapped_KyleNAc_raw.Rdata")
CPM.Cutoff    <- 2
lib.size      <- colSums(raw)
MedianLibSize <- median(lib.size)
min.count     <- ceiling(CPM.Cutoff / 1e6 * MedianLibSize)

filter <- edgeR::filterByExpr(raw,
                              design    = model.matrix(~clin_Kyle$group),
                              min.count = min.count)
Filtered_raw <- raw[filter, ]
Filtered_raw <- Filtered_raw[, rownames(clin_Kyle)]
stopifnot(all(colnames(Filtered_raw) == rownames(clin_Kyle)))

Filtered_raw_log <- edgeR::cpm(Filtered_raw, log = TRUE, prior.count = 1)

## Residualize against the six covariates while preserving MDD signal
design_keep <- model.matrix(~ y)
design_cov  <- model.matrix(~ ., data = covars)[, -1]

expr_resid <- removeBatchEffect(
  x          = as.matrix(Filtered_raw_log),
  covariates = design_cov,
  design     = design_keep
)
colnames(expr_resid) = paste0("S", colnames(expr_resid))
cat("Expression residualization complete. Dimensions:", dim(expr_resid), "\n")

rm(raw, Filtered_raw, Filtered_raw_log, filter); gc()
## ============================================================================

## ===========================================================================
## STEP 0.5: Build CpG annotation table (anno_long)
## ===========================================================================
## Reproduces the parsing from your earlier code: parses EPIC v2 manifest,
## filters to Type I/II probes, expands multi-gene CpGs into long format.

library(dplyr)
library(tidyr)
library(stringr)

## Path to the EPIC v2 manifest CSV
epic_manifest_path <- "~/Documents/Lab/Collaboration/DNMT/data/EPIC-8v2-0_A1.csv"

anno_df <- read.csv(epic_manifest_path, skip = 7, stringsAsFactors = FALSE)

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

## Keep only Type I and Type II probes
anno_df1 <- anno_df1[anno_df1$Type %in% c("I", "II"), ]

## Expand multi-gene CpG annotations into long format
anno_long <- anno_df1 %>%
  filter(!is.na(Gene), Gene != "") %>%
  separate_rows(Gene, Group, sep = ";") %>%
  mutate(Group = str_trim(Group)) %>%
  filter(Group != "") %>%
  distinct()

cat("anno_long built:", nrow(anno_long), "rows,",
    length(unique(anno_long$IlmnID)), "unique CpGs,",
    length(unique(anno_long$Gene)), "unique genes\n")

rm(anno_df, anno_df1); gc()
## ============================================================================


##
## Step 2 of CStaR pipeline: Local CpG-Gene-MDD causal networks
## ============================================================
##
## For each of 8 CStaR-prioritized genes, build a local multi-omic network
## (CpGs + gene expression + MDD) and run three causal discovery algorithms:
##   - Latent PC  (constraint-based, via Kendall's tau bridge functions)
##   - Latent FCI (handles unmeasured confounders, outputs PAG)
##   - GES        (score-based comparison)
##
## Inputs assumed to be available:
##   - swappedmislable_bmiq_combatmet_m_NAc.Rdata: M-value matrix (probes x 86)
##   - expr_resid: residualized log2-CPM matrix (genes x 86) from CStaR step
##   - clin_Kyle: clinical data with binary MDD coding
##   - covars: same 6-covariate frame used in CStaR (Age, Sex, BMI, PMI, pH, RIN)
##   - y: binary MDD vector (0/1) from CStaR step
##   - anno_long: long-format CpG-gene-region annotation from EPIC v2 manifest
##
## latentPC source files at:
##   /Users/haowang/Documents/Lab/Collaboration/DNMT/code/final_project/code
###############################################################################


## Source latentPC functions ==================================================
latentPC_path <- "/Users/haowang/Documents/Lab/Collaboration/DNMT/code/final_project/code"
source(file.path(latentPC_path, "utility.R"))
source(file.path(latentPC_path, "latent_pc.R"))
## ============================================================================


## Dependencies ===============================================================
library(pcalg)
library(limma)
library(mvtnorm)
library(pcaPP)
library(Matrix)
library(org.Hs.eg.db)   # for ENSG -> SYMBOL mapping
## ============================================================================


## Output directory for results ===============================================
results_dir <- "/Users/haowang/Documents/Lab/Collaboration/DNMT/code/final_project/local_networks_results"
dir.create(results_dir, showWarnings = FALSE)
## ============================================================================


## ===========================================================================
## STEP 1: Load and align methylation data, residualize against covariates
## ===========================================================================

## --- 1.1 Load M-value matrix ------------------------------------------------
m_mat <- readRDS("~/Documents/Lab/Collaboration/DNMT/data/swappedmislable_bmiq_combatmet_m_NAc.Rdata")
cat("M-value matrix:", nrow(m_mat), "probes x", ncol(m_mat), "samples\n")

## --- 1.2 Verify sample alignment with expression and clinical ---------------
## Both expr_resid and m_mat have samples ending in "N" (NAc); rownames(clin_Kyle)
## should match. Reorder m_mat to match clin_Kyle / expr_resid order.
stopifnot(exists("expr_resid"))
stopifnot(exists("clin_Kyle"))
stopifnot(exists("y"))
stopifnot(exists("covars"))

common_samples <- intersect(colnames(m_mat), rownames(clin_Kyle))
cat("Common samples:", length(common_samples), "\n")
stopifnot(length(common_samples) == ncol(expr_resid))

m_mat <- m_mat[, rownames(clin_Kyle)]
expr_resid <- expr_resid[, rownames(clin_Kyle)]
stopifnot(all(colnames(m_mat) == rownames(clin_Kyle)))
stopifnot(all(colnames(m_mat) == colnames(expr_resid)))
cat("Sample alignment verified.\n")

## --- 1.3 Verify no missingness ----------------------------------------------
n_na <- sum(is.na(m_mat))
cat("NAs in M-value matrix:", n_na, "\n")
if (n_na > 0) {
  warning("M-value matrix contains NAs; consider imputation before residualization.")
}

## --- 1.4 Residualize M-values against same six covariates -------------------
design_keep <- model.matrix(~ y)
design_cov  <- model.matrix(~ ., data = covars)[, -1]

cat("Residualizing M-values (this may take a few minutes)...\n")
m_resid <- removeBatchEffect(
  x          = as.matrix(m_mat),
  covariates = design_cov,
  design     = design_keep
)
rm(m_mat); gc()
cat("Residualization complete. Dimensions:", dim(m_resid), "\n")
## ============================================================================


## ===========================================================================
## STEP 2: Per-gene CpG annotation lookup
## ===========================================================================

## --- 2.1/2.2 Define target genes and look up CpGs ---------------------------
target_genes <- c("ARRDC3", "KAT5", "KDM4A", "SAT2",
                  "PRKAB2", "USE1", "SLC52A3", "ITGA1")

stopifnot(exists("anno_long"))   # parsed annotation table from earlier

## For each gene, collect annotated CpGs with collapsed region labels.
gene_cpgs <- list()
for (g in target_genes) {
  rows <- anno_long[anno_long$Gene == g, ]
  if (nrow(rows) == 0) {
    warning("No CpGs annotated to ", g, " in anno_long.")
    gene_cpgs[[g]] <- data.frame()
    next
  }
  ## Collapse multiple region annotations per CpG into one row
  collapsed <- rows %>%
    group_by(IlmnID) %>%
    summarise(
      regions = paste(sort(unique(Group)), collapse = ";"),
      island  = paste(sort(unique(Island)), collapse = ";"),
      chr     = first(CHR),
      pos     = first(MAPINFO),
      .groups = "drop"
    ) %>%
    as.data.frame()
  gene_cpgs[[g]] <- collapsed
}

## --- 2.3 Cross-reference with available CpGs in m_resid ---------------------
available_probes <- rownames(m_resid)
for (g in target_genes) {
  df <- gene_cpgs[[g]]
  if (nrow(df) == 0) next
  df <- df[df$IlmnID %in% available_probes, ]
  gene_cpgs[[g]] <- df
}

## --- 2.4 Document per-gene CpG counts ---------------------------------------
cpg_summary <- data.frame(
  gene        = target_genes,
  n_cpgs      = sapply(gene_cpgs, nrow),
  stringsAsFactors = FALSE
)
cat("\nPer-gene CpG counts:\n")
print(cpg_summary)
write.csv(cpg_summary, file.path(results_dir, "01_cpg_counts_per_gene.csv"),
          row.names = FALSE)

## Save full per-gene CpG metadata for later regional interpretation
saveRDS(gene_cpgs, file.path(results_dir, "02_gene_cpgs_metadata.rds"))
## ============================================================================


## ===========================================================================
## STEP 2.5: Map gene symbols to ENSG IDs in expr_resid
## ===========================================================================
## RNA-seq rownames are ENSG; methylation annotation uses gene symbols.
## Build a symbol -> ENSG mapping for our 8 target genes.

ensg_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = target_genes,
  columns  = c("SYMBOL", "ENSEMBL"),
  keytype  = "SYMBOL"
)
cat("\nENSG mapping for target genes:\n"); print(ensg_map)

## Some symbols may map to multiple ENSGs; keep only those present in expr_resid.
## Strip version suffixes from expr_resid rownames if present (e.g., "ENSG000001.5").
expr_ensg <- sub("\\..*$", "", rownames(expr_resid))
ensg_to_row <- setNames(seq_along(expr_ensg), expr_ensg)

gene_to_ensg <- list()
for (g in target_genes) {
  candidates <- ensg_map$ENSEMBL[ensg_map$SYMBOL == g]
  candidates <- candidates[!is.na(candidates) & candidates %in% expr_ensg]
  if (length(candidates) == 0) {
    warning("No ENSG match in expr_resid for ", g)
    gene_to_ensg[[g]] <- NA
  } else if (length(candidates) > 1) {
    ## Multiple ENSGs for one symbol: pick the one with highest mean expression
    means <- sapply(candidates, function(e) {
      mean(expr_resid[ensg_to_row[e], ], na.rm = TRUE)
    })
    chosen <- candidates[which.max(means)]
    cat("  ", g, ": multiple ENSGs (", paste(candidates, collapse = ","),
        "); chose", chosen, "by mean expression\n")
    gene_to_ensg[[g]] <- chosen
  } else {
    gene_to_ensg[[g]] <- candidates
  }
}
cat("Final gene -> ENSG mapping:\n"); print(unlist(gene_to_ensg))
## ============================================================================


## ===========================================================================
## STEP 3: Build local data matrix per gene
## ===========================================================================

build_local_matrix <- function(gene, gene_cpgs_meta, m_resid, expr_resid,
                               y, gene_to_ensg, expr_ensg) {
  cpg_ids <- gene_cpgs_meta$IlmnID
  if (length(cpg_ids) == 0) return(NULL)
  
  ensg <- gene_to_ensg[[gene]]
  if (is.na(ensg)) return(NULL)
  expr_row <- which(expr_ensg == ensg)
  if (length(expr_row) == 0) return(NULL)
  
  ## CpGs: samples x k
  cpg_block <- t(m_resid[cpg_ids, , drop = FALSE])
  ## Gene expression: samples x 1
  gene_block <- expr_resid[expr_row, , drop = FALSE]
  gene_block <- t(gene_block)
  colnames(gene_block) <- gene   # use symbol as column name
  
  ## MDD: samples x 1, integer 0/1
  mdd_block <- matrix(y, ncol = 1, dimnames = list(NULL, "MDD"))
  
  out <- cbind(cpg_block, gene_block, mdd_block)
  out
}

local_data <- list()
for (g in target_genes) {
  mat <- build_local_matrix(g, gene_cpgs[[g]], m_resid, expr_resid,
                            y, gene_to_ensg, expr_ensg)
  if (is.null(mat)) {
    cat("Skipping", g, "- no usable data\n")
    local_data[[g]] <- NULL
  } else {
    local_data[[g]] <- mat
    cat(sprintf("Local matrix for %s: %d samples x %d variables (%d CpGs + gene + MDD)\n",
                g, nrow(mat), ncol(mat), ncol(mat) - 2))
  }
}
saveRDS(local_data, file.path(results_dir, "03_local_data_matrices.rds"))
## ============================================================================


## ===========================================================================
## STEP 4: Helper for MDD-outgoing-edge constraint
## ===========================================================================
## Post-process an adjacency matrix to forbid edges from MDD to other variables.
## Removes any "MDD -> X" edge by zeroing the corresponding entry.
## In pcalg adjacency matrix convention: amat[i, j] = 1 means i -> j.
remove_outgoing_from_MDD <- function(amat, mdd_name = "MDD") {
  if (is.null(rownames(amat)) || is.null(colnames(amat))) {
    stop("Adjacency matrix must have row/column names.")
  }
  if (!mdd_name %in% rownames(amat)) return(amat)
  i <- which(rownames(amat) == mdd_name)
  ## Zero out the MDD row (no outgoing edges from MDD).
  ## Note: this preserves edges INTO MDD (which is what we want).
  amat[i, ] <- 0
  amat
}

## Same idea for FCI's PAG (which has different mark coding).
## In a PAG amat: amat[i,j] = 0 (no edge), 1 (circle), 2 (arrowhead), 3 (tail).
## Edge i -> j is encoded as amat[j,i] = 2 (arrowhead at j) and amat[i,j] = 3 (tail at i).
## We remove edges where MDD has a tail (i.e., MDD is a starting point of a directed edge).
remove_outgoing_from_MDD_pag <- function(pag, mdd_name = "MDD") {
  if (!mdd_name %in% rownames(pag)) return(pag)
  i <- which(rownames(pag) == mdd_name)
  ## For each j: if there's an edge i -> j, the marks are pag[i,j] = 3 (tail at MDD)
  ## and pag[j,i] = 2 (arrowhead at j). Remove these edges entirely.
  for (j in seq_len(ncol(pag))) {
    if (j == i) next
    if (pag[i, j] == 3 && pag[j, i] == 2) {
      pag[i, j] <- 0
      pag[j, i] <- 0
    }
  }
  pag
}
## ============================================================================


## ===========================================================================
## STEP 5: Run latent PC on each gene's local network
## ===========================================================================
alpha_local <- 0.05

run_latent_pc <- function(dat, gene_name, alpha = 0.05) {
  labels <- label_fun(dat)
  ## Sanity check: MDD should be binary
  mdd_idx <- which(colnames(dat) == "MDD")
  if (labels[mdd_idx] != "binary") {
    warning(gene_name, ": MDD not labeled binary; got '", labels[mdd_idx], "'")
  }
  
  sig <- latent_pc(dat, labels)
  sig <- as.matrix(Matrix::nearPD(sig, corr = TRUE, maxit = 10000)$mat)
  rownames(sig) <- colnames(sig) <- colnames(dat)
  
  pcFit <- pc(suffStat   = list(C = sig, n = nrow(dat)),
              indepTest  = gaussCItest,
              p          = ncol(dat),
              alpha      = alpha,
              labels     = colnames(dat),
              skel.method = "stable.fast",
              verbose    = FALSE)
  
  amat <- as(pcFit, "amat")
  amat <- remove_outgoing_from_MDD(amat, "MDD")
  
  ## IDA effects from each predictor to MDD
  effects_to_mdd <- list()
  predictor_names <- setdiff(colnames(dat), "MDD")
  predictor_idx   <- match(predictor_names, colnames(dat))
  mdd_idx_local   <- which(colnames(dat) == "MDD")
  
  effects_mat <- tryCatch(
    idaFast(mdd_idx_local, predictor_idx, sig, pcFit@graph),
    error = function(e) { warning("idaFast failed: ", e$message); NULL }
  )
  
  if (!is.null(effects_mat)) {
    for (j in seq_along(predictor_names)) {
      ms <- effects_mat[j, ]
      ms <- ms[!is.na(ms)]
      effects_to_mdd[[predictor_names[j]]] <-
        if (length(ms) > 0) min(abs(ms)) else NA_real_
    }
  }
  
  list(
    sig            = sig,
    pcFit          = pcFit,
    amat           = amat,
    effects_to_mdd = effects_to_mdd
  )
}

cat("\n=== Step 5: Latent PC ===\n")
pc_results <- list()
for (g in target_genes) {
  if (is.null(local_data[[g]])) next
  cat("Running latent PC for", g, "...\n")
  pc_results[[g]] <- tryCatch(
    run_latent_pc(local_data[[g]], g, alpha = alpha_local),
    error = function(e) { warning(g, ": PC failed: ", e$message); NULL }
  )
}
saveRDS(pc_results, file.path(results_dir, "04_latent_pc_results.rds"))
## ============================================================================


## ===========================================================================
## STEP 6: Run latent FCI on each gene's local network
## ===========================================================================

run_latent_fci <- function(dat, gene_name, sig, alpha = 0.05) {
  ## Reuse sig from PC step (latent correlation matrix)
  fciFit <- fci(suffStat   = list(C = sig, n = nrow(dat)),
                indepTest  = gaussCItest,
                p          = ncol(dat),
                alpha      = alpha,
                labels     = colnames(dat),
                skel.method = "stable.fast",
                verbose    = FALSE)
  
  pag  <- fciFit@amat
  pag  <- remove_outgoing_from_MDD_pag(pag, "MDD")
  
  list(fciFit = fciFit, pag = pag)
}

cat("\n=== Step 6: Latent FCI ===\n")
fci_results <- list()
for (g in target_genes) {
  if (is.null(pc_results[[g]])) next
  cat("Running latent FCI for", g, "...\n")
  fci_results[[g]] <- tryCatch(
    run_latent_fci(local_data[[g]], g, pc_results[[g]]$sig, alpha = alpha_local),
    error = function(e) { warning(g, ": FCI failed: ", e$message); NULL }
  )
}
saveRDS(fci_results, file.path(results_dir, "05_latent_fci_results.rds"))
## ============================================================================


## ===========================================================================
## STEP 7: Run GES on each gene's local network
## ===========================================================================

run_ges <- function(dat, gene_name) {
  ## GES uses Gaussian score; treats binary MDD as continuous (known compromise).
  score <- new("GaussL0penObsScore", as.matrix(dat))
  gesFit <- ges(score, verbose = FALSE)
  
  amat <- as(gesFit$essgraph, "matrix") * 1   # logical -> numeric
  rownames(amat) <- colnames(amat) <- colnames(dat)
  amat <- remove_outgoing_from_MDD(amat, "MDD")
  
  list(gesFit = gesFit, amat = amat)
}

cat("\n=== Step 7: GES ===\n")
ges_results <- list()
for (g in target_genes) {
  if (is.null(local_data[[g]])) next
  cat("Running GES for", g, "...\n")
  ges_results[[g]] <- tryCatch(
    run_ges(local_data[[g]], g),
    error = function(e) { warning(g, ": GES failed: ", e$message); NULL }
  )
}
saveRDS(ges_results, file.path(results_dir, "06_ges_results.rds"))
## ============================================================================


## ===========================================================================
## STEP 8: Per-gene synthesis
## ===========================================================================

## --- 8.1 Edge extraction helpers --------------------------------------------
## Convert PC/GES adjacency matrix to edge data frame.
## PC convention: amat[i,j] = 1 means i -> j; if also amat[j,i] = 1 then i - j (undirected).
amat_to_edges <- function(amat, method) {
  edges <- data.frame()
  vars <- rownames(amat)
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i >= j) next   # avoid double-counting undirected
      i_to_j <- amat[i, j] == 1
      j_to_i <- amat[j, i] == 1
      if (i_to_j && j_to_i) {
        edges <- rbind(edges, data.frame(from = vars[i], to = vars[j],
                                         type = "undirected",
                                         method = method,
                                         stringsAsFactors = FALSE))
      } else if (i_to_j) {
        edges <- rbind(edges, data.frame(from = vars[i], to = vars[j],
                                         type = "directed",
                                         method = method,
                                         stringsAsFactors = FALSE))
      } else if (j_to_i) {
        edges <- rbind(edges, data.frame(from = vars[j], to = vars[i],
                                         type = "directed",
                                         method = method,
                                         stringsAsFactors = FALSE))
      }
    }
  }
  edges
}

## Convert FCI PAG to edge data frame with mark types.
## PAG mark coding: 0 = no edge, 1 = circle, 2 = arrowhead, 3 = tail.
pag_to_edges <- function(pag, method = "FCI") {
  edges <- data.frame()
  vars <- rownames(pag)
  mark_label <- function(m) c("0" = "no", "1" = "circle",
                              "2" = "arrow", "3" = "tail")[as.character(m)]
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i >= j) next
      mij <- pag[i, j]; mji <- pag[j, i]
      if (mij == 0 && mji == 0) next
      ## Classify edge type
      type <- if (mij == 3 && mji == 2) "directed (i->j)"
      else if (mij == 2 && mji == 3) "directed (j->i)"
      else if (mij == 2 && mji == 2) "bidirected"
      else if (mij == 1 && mji == 2) "circle-arrow (j->? from i?)"
      else if (mij == 2 && mji == 1) "circle-arrow (i->? from j?)"
      else if (mij == 1 && mji == 1) "circle-circle"
      else paste0("other(", mij, ",", mji, ")")
      edges <- rbind(edges, data.frame(
        from = vars[i], to = vars[j], type = type,
        mark_at_from = mark_label(mij), mark_at_to = mark_label(mji),
        method = method, stringsAsFactors = FALSE))
    }
  }
  edges
}

## --- 8.2 Cascade classification helper --------------------------------------
classify_cascade <- function(edges_df, gene_name) {
  ## Look for any CpG -> gene directed edge AND gene -> MDD directed edge.
  ## "edges_df" can be the PC/GES edge frame or FCI edge frame (with type column).
  cpg_cols <- setdiff(unique(c(edges_df$from, edges_df$to)),
                      c(gene_name, "MDD"))
  
  cpg_to_gene_directed <- any(
    edges_df$from %in% cpg_cols & edges_df$to == gene_name &
      grepl("directed", edges_df$type) & !grepl("circle", edges_df$type)
  )
  cpg_to_gene_any <- any(
    (edges_df$from %in% cpg_cols & edges_df$to == gene_name) |
      (edges_df$to %in% cpg_cols & edges_df$from == gene_name)
  )
  cpg_to_gene_bidirected <- any(
    ((edges_df$from %in% cpg_cols & edges_df$to == gene_name) |
       (edges_df$to %in% cpg_cols & edges_df$from == gene_name)) &
      edges_df$type == "bidirected"
  )
  
  gene_to_mdd_directed <- any(
    edges_df$from == gene_name & edges_df$to == "MDD" &
      grepl("directed", edges_df$type) & !grepl("circle", edges_df$type)
  )
  gene_to_mdd_any <- any(
    (edges_df$from == gene_name & edges_df$to == "MDD") |
      (edges_df$to == gene_name & edges_df$from == "MDD")
  )
  
  if (cpg_to_gene_directed && gene_to_mdd_directed) {
    return("strong support")
  } else if (cpg_to_gene_bidirected) {
    return("inconsistent (CpG-gene confounded)")
  } else if (cpg_to_gene_any && gene_to_mdd_any) {
    return("partial support")
  } else if (gene_to_mdd_any && !cpg_to_gene_any) {
    return("inconsistent (no CpG->gene)")
  } else if (!cpg_to_gene_any && !gene_to_mdd_any) {
    return("no structure")
  } else {
    return("partial structure")
  }
}

## --- 8.3 Build per-gene synthesis -------------------------------------------
synthesis <- list()

for (g in target_genes) {
  if (is.null(local_data[[g]])) {
    synthesis[[g]] <- list(status = "skipped (no data)")
    next
  }
  
  ## Edges from each method
  pc_edges <- if (!is.null(pc_results[[g]])) {
    amat_to_edges(pc_results[[g]]$amat, method = "latent_PC")
  } else data.frame()
  
  fci_edges <- if (!is.null(fci_results[[g]])) {
    pag_to_edges(fci_results[[g]]$pag, method = "latent_FCI")
  } else data.frame()
  
  ges_edges <- if (!is.null(ges_results[[g]])) {
    amat_to_edges(ges_results[[g]]$amat, method = "GES")
  } else data.frame()
  
  ## Cascade classification per method
  pc_class  <- if (nrow(pc_edges)  > 0) classify_cascade(pc_edges, g)  else "no result"
  fci_class <- if (nrow(fci_edges) > 0) classify_cascade(fci_edges, g) else "no result"
  ges_class <- if (nrow(ges_edges) > 0) classify_cascade(ges_edges, g) else "no result"
  
  ## Specific CpGs flagged: those with directed edge to gene under any method
  cpg_meta <- gene_cpgs[[g]]
  flagged_cpgs <- character(0)
  
  for (df in list(pc_edges, fci_edges, ges_edges)) {
    if (nrow(df) == 0) next
    hits <- df$from[df$to == g & grepl("directed", df$type) & !grepl("circle", df$type)]
    flagged_cpgs <- union(flagged_cpgs, hits)
  }
  
  flagged_meta <- if (length(flagged_cpgs) > 0) {
    cpg_meta[match(flagged_cpgs, cpg_meta$IlmnID), ]
  } else data.frame()
  
  ## IDA effects: gene -> MDD and CpGs -> MDD from latent PC
  ida_effects <- if (!is.null(pc_results[[g]])) {
    pc_results[[g]]$effects_to_mdd
  } else list()
  
  ## Confidence label
  classes <- c(pc_class, fci_class, ges_class)
  strong_count <- sum(classes == "strong support")
  bidirected_count <- sum(grepl("confounded", classes))
  
  confidence <- if (strong_count == 3) {
    "High"
  } else if (strong_count >= 2 && bidirected_count == 0) {
    "Moderate"
  } else if (bidirected_count >= 1) {
    "Low (confounding flagged)"
  } else if (any(classes == "no structure")) {
    "Inconclusive"
  } else {
    "Low"
  }
  
  synthesis[[g]] <- list(
    n_cpgs           = nrow(cpg_meta),
    pc_class         = pc_class,
    fci_class        = fci_class,
    ges_class        = ges_class,
    confidence       = confidence,
    flagged_cpgs     = flagged_cpgs,
    flagged_metadata = flagged_meta,
    ida_effects_to_mdd = ida_effects,
    edges_pc         = pc_edges,
    edges_fci        = fci_edges,
    edges_ges        = ges_edges
  )
  
  cat(sprintf("\n%s: %d CpGs | PC=%s | FCI=%s | GES=%s | confidence=%s\n",
              g, nrow(cpg_meta), pc_class, fci_class, ges_class, confidence))
  if (length(flagged_cpgs) > 0) {
    cat("  Flagged CpGs (directed -> gene):", paste(flagged_cpgs, collapse = ", "), "\n")
  }
}

saveRDS(synthesis, file.path(results_dir, "07_per_gene_synthesis.rds"))


## --- 8.4 Write summary tables -----------------------------------------------
## Top-level summary across genes
summary_tbl <- data.frame(
  gene       = target_genes,
  n_cpgs     = sapply(synthesis, function(s) if (is.null(s$n_cpgs)) NA else s$n_cpgs),
  pc_class   = sapply(synthesis, function(s) if (is.null(s$pc_class))  NA else s$pc_class),
  fci_class  = sapply(synthesis, function(s) if (is.null(s$fci_class)) NA else s$fci_class),
  ges_class  = sapply(synthesis, function(s) if (is.null(s$ges_class)) NA else s$ges_class),
  confidence = sapply(synthesis, function(s) if (is.null(s$confidence)) NA else s$confidence),
  n_flagged_cpgs = sapply(synthesis, function(s) length(s$flagged_cpgs)),
  stringsAsFactors = FALSE
)
write.csv(summary_tbl, file.path(results_dir, "08_summary_per_gene.csv"),
          row.names = FALSE)

cat("\n=== Final per-gene summary ===\n")
print(summary_tbl)


## Combined edge list across all genes and methods
all_edges <- list()
for (g in target_genes) {
  s <- synthesis[[g]]
  if (is.null(s$edges_pc)  && is.null(s$edges_fci) && is.null(s$edges_ges)) next
  for (df_name in c("edges_pc", "edges_fci", "edges_ges")) {
    if (!is.null(s[[df_name]]) && nrow(s[[df_name]]) > 0) {
      df <- s[[df_name]]
      df$target_gene <- g
      all_edges[[paste(g, df_name)]] <- df
    }
  }
}
all_edges_df <- do.call(rbind, lapply(all_edges, function(x) {
  ## Pad columns so rbind works across PC/GES (no marks) and FCI (marks)
  needed <- c("from", "to", "type", "method", "target_gene",
              "mark_at_from", "mark_at_to")
  for (col in needed) if (is.null(x[[col]])) x[[col]] <- NA_character_
  x[, needed]
}))
write.csv(all_edges_df, file.path(results_dir, "09_all_edges.csv"),
          row.names = FALSE)


## Flagged CpG details across all genes
flagged_all <- list()
for (g in target_genes) {
  s <- synthesis[[g]]
  if (length(s$flagged_cpgs) == 0) next
  df <- s$flagged_metadata
  if (nrow(df) == 0) next
  df$target_gene <- g
  ## Add IDA effect of CpG on MDD if available
  df$ida_effect_to_mdd <- sapply(df$IlmnID, function(cpg) {
    eff <- s$ida_effects_to_mdd[[cpg]]
    if (is.null(eff)) NA_real_ else eff
  })
  flagged_all[[g]] <- df
}
if (length(flagged_all) > 0) {
  flagged_df <- do.call(rbind, flagged_all)
  write.csv(flagged_df, file.path(results_dir, "10_flagged_cpgs.csv"),
            row.names = FALSE)
}

cat("\nResults saved to:", results_dir, "\n")
cat("Files written:\n")
print(list.files(results_dir))
## ============================================================================



##
## Step 8.5: Visualize per-gene local networks (hybrid renderer, dot source)
## =========================================================================
## Hybrid plotting:
##   - Latent PC and GES panels: rendered via ggraph
##   - Latent FCI panel: rendered via raw Graphviz dot source (Tetrad-style)
##
## Why dot source: Rgraphviz's per-edge attribute system silently strips
## attributes when parsing fails, leading to edges being rendered without
## marks (or not at all). Writing dot source directly and shelling out to
## the system `dot` command bypasses these issues.
##
## Edge styling:
##   - CpG-gene edges:  blue
##   - gene-MDD edges:  red
##   - CpG-MDD edges:   purple
##   - CpG-CpG edges:   grey, thin (de-emphasized)
##
## Edge mark types (FCI panel only):
##   - directed (A -> B):       tail -> arrowhead
##   - bidirected (A <-> B):    arrowhead <-> arrowhead
##   - circle-arrow (A o-> B):  open circle -> arrowhead
##   - circle-circle (A o-o B): open circle <-> open circle
##
## Required packages:
##   install.packages(c("igraph", "ggraph", "ggplot2", "dplyr", "patchwork",
##                      "magick"))
##
## Required system: Graphviz `dot` command on PATH.
##   macOS:  brew install graphviz
##   Ubuntu: sudo apt-get install graphviz
##
## Verify with: system("dot -V")
###############################################################################


## Dependencies ===============================================================
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(patchwork)
library(magick)

## Make sure Homebrew Graphviz is on PATH (macOS).
## Adjust if you installed Graphviz elsewhere.
if (Sys.info()["sysname"] == "Darwin") {
  brew_paths <- c("/opt/homebrew/bin", "/usr/local/bin")
  for (bp in brew_paths) {
    if (!grepl(bp, Sys.getenv("PATH"), fixed = TRUE) && dir.exists(bp)) {
      Sys.setenv(PATH = paste(bp, Sys.getenv("PATH"), sep = ":"))
    }
  }
}

## Verify dot is accessible
dot_check <- suppressWarnings(system("dot -V", intern = TRUE,
                                     ignore.stderr = FALSE))
if (length(dot_check) == 0) {
  warning("Graphviz `dot` command not found on PATH. ",
          "FCI panels will not render. ",
          "Install with: brew install graphviz (macOS) or ",
          "sudo apt-get install graphviz (Linux).")
}
## ============================================================================


## Output directory for plots =================================================
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
## ============================================================================


## --- Helper: classify each edge by role in the cascade ---------------------
classify_edge_role <- function(from, to, gene_name) {
  is_gene <- function(x) x == gene_name
  is_mdd  <- function(x) x == "MDD"
  is_cpg  <- function(x) !is_gene(x) && !is_mdd(x)
  
  if      ((is_cpg(from)  && is_gene(to)) || (is_gene(from) && is_cpg(to)))  "CpG-gene"
  else if ((is_gene(from) && is_mdd(to))  || (is_mdd(from)  && is_gene(to))) "gene-MDD"
  else if ((is_cpg(from)  && is_mdd(to))  || (is_mdd(from)  && is_cpg(to)))  "CpG-MDD"
  else                                                                       "CpG-CpG"
}


## --- Helper: build igraph from PC/GES adjacency matrix ---------------------
amat_to_igraph <- function(amat, gene_name) {
  vars <- rownames(amat)
  edges <- list()
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i >= j) next
      i_to_j <- amat[i, j] == 1
      j_to_i <- amat[j, i] == 1
      if (i_to_j && j_to_i) {
        edges[[length(edges) + 1]] <- data.frame(
          from = vars[i], to = vars[j],
          edge_type = "undirected", arrow_style = "none",
          stringsAsFactors = FALSE)
      } else if (i_to_j) {
        edges[[length(edges) + 1]] <- data.frame(
          from = vars[i], to = vars[j],
          edge_type = "directed", arrow_style = "to",
          stringsAsFactors = FALSE)
      } else if (j_to_i) {
        edges[[length(edges) + 1]] <- data.frame(
          from = vars[j], to = vars[i],
          edge_type = "directed", arrow_style = "to",
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(edges) == 0) {
    edges_df <- data.frame(from = character(), to = character(),
                           edge_type = character(), arrow_style = character(),
                           stringsAsFactors = FALSE)
  } else {
    edges_df <- do.call(rbind, edges)
  }
  edges_df$role <- mapply(classify_edge_role, edges_df$from, edges_df$to, gene_name)
  
  g <- graph_from_data_frame(edges_df, directed = TRUE, vertices = vars)
  list(graph = g, edges_df = edges_df)
}


## --- Helper: extract edges from FCI PAG ------------------------------------
## In pcalg PAG matrix convention: pag[i, j] is the mark AT vertex j on the
## edge between i and j.
##   mark_at_i = pag[j, i]
##   mark_at_j = pag[i, j]
pag_to_edge_list <- function(pag) {
  vars <- rownames(pag)
  edges <- list()
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i >= j) next
      mij <- pag[i, j]; mji <- pag[j, i]
      if (mij == 0 && mji == 0) next
      edges[[length(edges) + 1]] <- list(
        from = vars[i], to = vars[j],
        mark_at_from = mji,
        mark_at_to   = mij
      )
    }
  }
  edges
}


## --- Helper: build igraph from FCI PAG (used for cascade subgraph filter) --
pag_to_igraph <- function(pag, gene_name) {
  vars <- rownames(pag)
  edges <- list()
  for (i in seq_along(vars)) {
    for (j in seq_along(vars)) {
      if (i >= j) next
      mij <- pag[i, j]; mji <- pag[j, i]
      if (mij == 0 && mji == 0) next
      
      mark_at_i <- mji
      mark_at_j <- mij
      
      if (mark_at_i == 3 && mark_at_j == 2) {
        from <- vars[i]; to <- vars[j]
        edge_type <- "directed"; arrow_style <- "to"
      } else if (mark_at_i == 2 && mark_at_j == 3) {
        from <- vars[j]; to <- vars[i]
        edge_type <- "directed"; arrow_style <- "to"
      } else if (mark_at_i == 2 && mark_at_j == 2) {
        from <- vars[i]; to <- vars[j]
        edge_type <- "bidirected"; arrow_style <- "both"
      } else if (mark_at_i == 1 && mark_at_j == 2) {
        from <- vars[i]; to <- vars[j]
        edge_type <- "circle-arrow"; arrow_style <- "to"
      } else if (mark_at_i == 2 && mark_at_j == 1) {
        from <- vars[j]; to <- vars[i]
        edge_type <- "circle-arrow"; arrow_style <- "to"
      } else if (mark_at_i == 1 && mark_at_j == 1) {
        from <- vars[i]; to <- vars[j]
        edge_type <- "circle-circle"; arrow_style <- "none"
      } else {
        from <- vars[i]; to <- vars[j]
        edge_type <- paste0("other(", mark_at_i, ",", mark_at_j, ")")
        arrow_style <- "none"
      }
      
      edges[[length(edges) + 1]] <- data.frame(
        from = from, to = to,
        edge_type = edge_type, arrow_style = arrow_style,
        stringsAsFactors = FALSE)
    }
  }
  if (length(edges) == 0) {
    edges_df <- data.frame(from = character(), to = character(),
                           edge_type = character(), arrow_style = character(),
                           stringsAsFactors = FALSE)
  } else {
    edges_df <- do.call(rbind, edges)
  }
  edges_df$role <- mapply(classify_edge_role, edges_df$from, edges_df$to, gene_name)
  
  g <- graph_from_data_frame(edges_df, directed = TRUE, vertices = vars)
  list(graph = g, edges_df = edges_df)
}


## --- Plot PC/GES panel via ggraph ------------------------------------------
plot_amat_ggraph <- function(graph_obj, edges_df, gene_name, method_label,
                             subgraph_only = FALSE) {
  vars <- V(graph_obj)$name
  
  if (subgraph_only) {
    relevant_cpgs <- unique(c(
      edges_df$from[edges_df$role %in% c("CpG-gene", "CpG-MDD") &
                      edges_df$from != gene_name & edges_df$from != "MDD"],
      edges_df$to[  edges_df$role %in% c("CpG-gene", "CpG-MDD") &
                      edges_df$to   != gene_name & edges_df$to   != "MDD"]
    ))
    keep <- vars %in% c(gene_name, "MDD", relevant_cpgs)
    graph_obj <- induced_subgraph(graph_obj, V(graph_obj)[keep])
    edges_df  <- edges_df[edges_df$from %in% V(graph_obj)$name &
                            edges_df$to   %in% V(graph_obj)$name, ]
    vars <- V(graph_obj)$name
  }
  
  if (length(vars) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste0(method_label, "\n(empty)"), size = 4) +
             ggtitle(method_label))
  }
  
  node_type <- ifelse(vars == gene_name, "Gene",
                      ifelse(vars == "MDD", "MDD", "CpG"))
  V(graph_obj)$node_type <- node_type
  
  display_label <- vars
  if (sum(node_type == "CpG") > 12) {
    display_label[node_type == "CpG"] <- substr(display_label[node_type == "CpG"], 1, 8)
  }
  V(graph_obj)$display_label <- display_label
  
  edge_color <- ifelse(edges_df$role == "CpG-gene", "#1f77b4",
                       ifelse(edges_df$role == "gene-MDD", "#d62728",
                              ifelse(edges_df$role == "CpG-MDD",  "#9467bd",
                                     "#cccccc")))
  edge_width <- ifelse(edges_df$role == "CpG-CpG", 0.3, 0.9)
  edge_alpha <- ifelse(edges_df$role == "CpG-CpG", 0.5, 1.0)
  edge_lty   <- ifelse(edges_df$edge_type == "bidirected",   "longdash", "solid")
  
  if (ecount(graph_obj) > 0) {
    E(graph_obj)$edge_color <- edge_color
    E(graph_obj)$edge_w     <- edge_width
    E(graph_obj)$edge_a     <- edge_alpha
    E(graph_obj)$edge_lty   <- edge_lty
    E(graph_obj)$has_arrow  <- edges_df$arrow_style == "to"
    E(graph_obj)$is_bidir   <- edges_df$arrow_style == "both"
  }
  
  set.seed(1)
  layout <- create_layout(graph_obj, layout = "fr")
  
  p <- ggraph(layout)
  
  if (ecount(graph_obj) > 0) {
    p <- p + geom_edge_link(
      aes(filter        = has_arrow,
          edge_colour   = edge_color,
          edge_width    = edge_w,
          edge_alpha    = edge_a,
          edge_linetype = edge_lty),
      arrow     = arrow(length = unit(2.5, "mm"), type = "closed"),
      end_cap   = circle(3, "mm"),
      start_cap = circle(3, "mm")
    )
    p <- p + geom_edge_link(
      aes(filter        = is_bidir,
          edge_colour   = edge_color,
          edge_width    = edge_w,
          edge_alpha    = edge_a,
          edge_linetype = edge_lty),
      arrow     = arrow(length = unit(2.5, "mm"), type = "closed", ends = "both"),
      end_cap   = circle(3, "mm"),
      start_cap = circle(3, "mm")
    )
    p <- p + geom_edge_link(
      aes(filter        = !has_arrow & !is_bidir,
          edge_colour   = edge_color,
          edge_width    = edge_w,
          edge_alpha    = edge_a,
          edge_linetype = edge_lty),
      end_cap   = circle(3, "mm"),
      start_cap = circle(3, "mm")
    )
  }
  
  p <- p +
    scale_edge_colour_identity() +
    scale_edge_width_identity() +
    scale_edge_alpha_identity() +
    scale_edge_linetype_identity() +
    geom_node_point(aes(color = node_type, size = node_type)) +
    geom_node_text(aes(label = display_label), repel = TRUE, size = 2.5) +
    scale_color_manual(values = c("Gene" = "#2ca02c",
                                  "MDD"  = "#d62728",
                                  "CpG"  = "#1f77b4"),
                       guide = "none") +
    scale_size_manual(values = c("Gene" = 6, "MDD" = 6, "CpG" = 3),
                      guide = "none") +
    theme_void() +
    theme(plot.title  = element_text(hjust = 0.5, size = 11, face = "bold"),
          plot.margin = margin(5, 5, 5, 5)) +
    ggtitle(method_label)
  
  p
}


## --- Render FCI PAG via raw Graphviz dot source ----------------------------
plot_fci_via_dot <- function(pag, gene_name, output_path,
                             format = c("pdf", "png"),
                             width_in = 6, height_in = 6,
                             subgraph_only = FALSE,
                             show_title = TRUE) {
  format <- match.arg(format)
  vars <- rownames(pag)
  edges <- pag_to_edge_list(pag)
  
  ## Optional cascade subgraph filter
  if (subgraph_only) {
    cascade_edges <- Filter(function(e) {
      role <- classify_edge_role(e$from, e$to, gene_name)
      role %in% c("CpG-gene", "gene-MDD", "CpG-MDD")
    }, edges)
    relevant <- unique(unlist(lapply(cascade_edges, function(e) c(e$from, e$to))))
    relevant <- setdiff(relevant, c(gene_name, "MDD"))
    keep_vars <- intersect(c(gene_name, "MDD", relevant), vars)
    edges <- Filter(function(e) e$from %in% keep_vars && e$to %in% keep_vars, edges)
    vars <- keep_vars
  }
  
  if (length(vars) == 0 || length(edges) == 0) {
    if (format == "pdf") {
      pdf(output_path, width = width_in, height = height_in)
    } else {
      png(output_path, width = width_in * 200, height = height_in * 200, res = 200)
    }
    plot.new()
    text(0.5, 0.5, paste0("Latent FCI: ", gene_name, "\n(no edges)"))
    dev.off()
    return(invisible(output_path))
  }
  
  ## Mark codes -> Graphviz arrow shape names
  mark_to_shape <- function(m) {
    switch(as.character(m),
           "1" = "odot",
           "2" = "normal",
           "3" = "none",
           "none")
  }
  
  ## Decide whether to shorten CpG labels
  is_cpg <- !(vars %in% c(gene_name, "MDD"))
  shorten <- sum(is_cpg) > 12
  
  ## Build node lines
  node_lines <- character(length(vars))
  for (i in seq_along(vars)) {
    v <- vars[i]
    fillcolor <- if (v == gene_name) "#2ca02c"
    else if (v == "MDD") "#d62728"
    else "#a6cee3"
    label <- if (shorten && !(v %in% c(gene_name, "MDD"))) substr(v, 1, 8) else v
    node_lines[i] <- sprintf(
      '  "%s" [label="%s", style=filled, fillcolor="%s", shape=ellipse, fontsize=11];',
      v, label, fillcolor
    )
  }
  
  ## Build edge lines
  edge_lines <- character(length(edges))
  for (i in seq_along(edges)) {
    e <- edges[[i]]
    arrowtail_shape <- mark_to_shape(e$mark_at_from)
    arrowhead_shape <- mark_to_shape(e$mark_at_to)
    
    role <- classify_edge_role(e$from, e$to, gene_name)
    color <- switch(role,
                    "CpG-gene" = "#1f77b4",
                    "gene-MDD" = "#d62728",
                    "CpG-MDD"  = "#9467bd",
                    "CpG-CpG"  = "#aaaaaa")
    penwidth <- if (role == "CpG-CpG") 1 else 2.5
    
    edge_lines[i] <- sprintf(
      '  "%s" -> "%s" [dir=both, arrowtail=%s, arrowhead=%s, color="%s", penwidth=%.1f];',
      e$from, e$to, arrowtail_shape, arrowhead_shape, color, penwidth
    )
  }
  
  ## Optional graph title
  graph_label <- if (show_title) {
    sprintf('  graph [rankdir=LR, label="Latent FCI: %s", labelloc=top, fontsize=14];',
            gene_name)
  } else {
    '  graph [rankdir=LR];'
  }
  
  ## Assemble dot source
  dot_source <- c(
    "digraph G {",
    graph_label,
    "  node [shape=ellipse, fontsize=11];",
    "  edge [arrowsize=0.8];",
    node_lines,
    edge_lines,
    "}"
  )
  
  ## Write dot file and render
  dot_file <- tempfile(fileext = ".dot")
  writeLines(dot_source, dot_file)
  
  cmd <- sprintf('dot -T%s -o "%s" "%s"', format, output_path, dot_file)
  status <- system(cmd, intern = FALSE)
  
  if (status != 0) {
    warning("dot command failed for ", gene_name,
            " (exit status ", status, "). ",
            "Verify Graphviz is installed and `dot` is on PATH. ",
            "Dot source saved at: ", dot_file)
  } else {
    file.remove(dot_file)
  }
  
  invisible(output_path)
}


## --- Convert a PNG file to a ggplot canvas for patchwork combination -------
png_to_ggplot <- function(png_path, title = NULL) {
  img <- magick::image_read(png_path)
  img_info <- magick::image_info(img)
  ratio <- img_info$height / img_info$width
  
  raster <- grid::rasterGrob(magick::image_read(png_path),
                             interpolate = TRUE,
                             width  = unit(1, "npc"),
                             height = unit(1, "npc"))
  
  p <- ggplot() +
    annotation_custom(raster, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    coord_fixed(ratio = ratio, xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
  
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  }
  p
}


## --- Build the 3-panel figure for one gene ---------------------------------
plot_gene_panel_hybrid <- function(gene_name, pc_results, fci_results, ges_results,
                                   subgraph_only = FALSE,
                                   tmp_dir = tempdir()) {
  
  pc_obj  <- pc_results[[gene_name]]
  fci_obj <- fci_results[[gene_name]]
  ges_obj <- ges_results[[gene_name]]
  
  pc_built <- if (!is.null(pc_obj))
    amat_to_igraph(pc_obj$amat, gene_name) else
      list(graph = make_empty_graph(0), edges_df = data.frame())
  p_pc <- plot_amat_ggraph(pc_built$graph, pc_built$edges_df,
                           gene_name, "Latent PC", subgraph_only)
  
  ges_built <- if (!is.null(ges_obj))
    amat_to_igraph(ges_obj$amat, gene_name) else
      list(graph = make_empty_graph(0), edges_df = data.frame())
  p_ges <- plot_amat_ggraph(ges_built$graph, ges_built$edges_df,
                            gene_name, "GES", subgraph_only)
  
  if (!is.null(fci_obj)) {
    suffix <- if (subgraph_only) "_cascade" else ""
    fci_png_path <- file.path(tmp_dir, paste0("FCI_", gene_name, suffix, ".png"))
    plot_fci_via_dot(fci_obj$pag, gene_name, fci_png_path,
                     format = "png", width_in = 6, height_in = 6,
                     subgraph_only = subgraph_only,
                     show_title = FALSE)
    p_fci <- png_to_ggplot(fci_png_path, title = "Latent FCI")
  } else {
    p_fci <- ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5,
               label = "Latent FCI\n(no result)", size = 4) +
      ggtitle("Latent FCI")
  }
  
  caption_text <- paste0(
    "Node colors: Gene (green), MDD (red), CpG (blue).  ",
    "Edge colors: CpG\u2013gene (blue), gene\u2013MDD (red), ",
    "CpG\u2013MDD (purple), CpG\u2013CpG (grey, de\u2011emphasized).\n",
    "FCI edge marks: filled arrow = directed (\u2192); ",
    "open circle = uncertain orientation; ",
    "arrows on both ends = bidirected (\u2194, latent confounding)."
  )
  
  combined <- p_pc + p_fci + p_ges +
    patchwork::plot_layout(ncol = 3) +
    patchwork::plot_annotation(
      title    = paste0("Local network for ", gene_name),
      subtitle = caption_text,
      theme    = theme(
        plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 8.5, lineheight = 1.2)
      )
    )
  combined
}


## --- Generate plots for all eight genes ------------------------------------
cat("\n=== Generating plots (hybrid renderer with dot source) ===\n")

fci_tmp_dir <- file.path(plots_dir, "fci_intermediate")
dir.create(fci_tmp_dir, showWarnings = FALSE, recursive = TRUE)

for (g in target_genes) {
  if (is.null(pc_results[[g]]) && is.null(fci_results[[g]]) && is.null(ges_results[[g]])) {
    cat("Skipping", g, "- no results\n")
    next
  }
  
  cat("Plotting", g, "...\n")
  
  ## Combined 3-panel figure (full network)
  fig <- tryCatch(
    plot_gene_panel_hybrid(g, pc_results, fci_results, ges_results,
                           subgraph_only = FALSE, tmp_dir = fci_tmp_dir),
    error = function(e) { warning("Plot failed for ", g, ": ", e$message); NULL }
  )
  if (!is.null(fig)) {
    out_path <- file.path(plots_dir, paste0(g, "_local_network.pdf"))
    ggsave(out_path, fig, width = 16, height = 6.5, device = "pdf")
  }
  
  ## Combined 3-panel figure (cascade only)
  fig_sub <- tryCatch(
    plot_gene_panel_hybrid(g, pc_results, fci_results, ges_results,
                           subgraph_only = TRUE, tmp_dir = fci_tmp_dir),
    error = function(e) NULL
  )
  if (!is.null(fig_sub)) {
    out_path_sub <- file.path(plots_dir, paste0(g, "_cascade_only.pdf"))
    ggsave(out_path_sub, fig_sub, width = 16, height = 6.5, device = "pdf")
  }
  
  ## Standalone Tetrad-style FCI panel as PDF
  if (!is.null(fci_results[[g]])) {
    fci_pdf_path <- file.path(plots_dir, paste0(g, "_FCI_only.pdf"))
    tryCatch(
      plot_fci_via_dot(fci_results[[g]]$pag, g, fci_pdf_path,
                       format = "pdf", width_in = 6, height_in = 6,
                       subgraph_only = FALSE, show_title = TRUE),
      error = function(e) warning("FCI standalone PDF failed for ", g, ": ", e$message)
    )
  }
}

cat("\nPlots saved to:", plots_dir, "\n")
print(list.files(plots_dir))

