# Causal Discovery in Multi-Omics Data for MDD

Code for the 80816 final project: applying causal discovery to 
postmortem brain RNA-seq and DNA methylation data to investigate 
molecular mechanisms of Major Depressive Disorder (MDD).

The pipeline has two steps:
1. **CStaR** on RNA-seq to prioritize genes by stable causal effect on MDD.
2. **Localized causal discovery** (latent PC, latent FCI, GES) on each 
   prioritized gene, its corresponding CpG Methylation sites, and MDD status.

## Files

| File | Purpose |
|---|---|
| `Meth_preprocess.R` | Methylation preprocessing (BMIQ, ComBat-Met, M-values, residualization). |
| `Cstar_latentPC.R` | CStaR with latent PC for binary MDD outcome. |
| `causal_discovery_meth_RNA.R` | Local CpG-gene-MDD networks per top gene. |
| `latent_pc.R` | Latent PC implementation helper functions (Cai et al. 2022). |
| `utility.R` | Helper functions Cai et al. 2022). |

Run in the order: `Meth_preprocess.R` → `Cstar_latentPC.R` → 
`causal_discovery_meth_RNA.R`.

## Dependencies

```r
install.packages(c("pcalg", "igraph", "ggraph", "ggplot2", "patchwork", 
                   "magick", "BiocManager"))
BiocManager::install(c("edgeR", "limma", "minfi", "sva"))
```

System: Graphviz (`brew install graphviz` on macOS) for FCI plot rendering.

## Data

Postmortem nucleus accumbens samples from 86 subjects (43 MDD, 43 controls), 
with 14,221 genes (RNA-seq) and 862,399 CpG sites (EPIC v2). Raw data are 
not in the repository.

## Authors

- Hao Wang (haow5@andrew.cmu.edu)
- Stephen V. Glass (svg13@pitt.edu)
