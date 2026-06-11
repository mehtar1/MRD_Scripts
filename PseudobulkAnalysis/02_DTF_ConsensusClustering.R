# =============================================================================
# Differential Transcription Factor Activity (DTF) Analysis per dASP
# =============================================================================
#
# Description:
#   Identifies transcription factors (TFs) with differential chromVAR
#   motif accessibility between MRD and IDX blast cells within each dASP.
#   Outputs a consensus matrix (DTFs × dASPs) for downstream consensus
#   clustering in 02_DTF_ConsensusCluster_Visualisation.Rmd.
#
#   Steps:
#     1. Load per-dASP chromVAR objects and ATAC Seurat objects
#     2. Build joint ATAC+chromVAR objects per dASP, filtered to blast cells
#     3. Run chromVAR FindMarkers (LR test, patient regressed) to call DTFs
#     4. Filter DTFs by gene expression (TF gene must be expressed in
#        >= 10% of cells in the relevant stage)
#     5. Aggregate DTF calls across all dASPs into a joint matrix
#     6. Build the mean z-score deviation matrix (DTFs × dASPs) used as
#        input for consensus clustering
#
# Inputs (set in USER CONFIGURATION):
#   - Per-dASP chromVAR SummarizedExperiment RDS files
#   - Per-dASP ATAC Seurat RDS files
#   - Global blast RNA Seurat object (for expression filtering)
#   - dASP cluster assignments CSV (from 01_dASP_Discovery_and_DEG_Analysis.Rmd)
#   - Blast cell metadata table
#   - Full cell metadata table (HTAN_MRD_AllCell_Metadata)
#
# Outputs:
#   - Per-dASP DTF tables (TrajCluster_<i>_PB_EnrichedTFs_MRDvsIDX.txt)
#   - DTFList.allClustersPB.rds   : all DTF calls combined (unfiltered by CC)
#   - ConsensusDTFlist.ClusterPB.rds : mean z-score matrix (DTFs × dASPs)
#   - ConsensusDTFlist.allClustersPB.rds : re-tested DTFs with scores per dASP
#
# Dependencies:
#   Seurat, Signac, chromVAR, dplyr, tidyverse, data.table, MAST,
#   RColorBrewer, BiocParallel, conflicted
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(gam)
library(Seurat)
library(Signac)
library(dplyr)
library(RColorBrewer)
library(BiocParallel)
library(ggplot2)
library(pheatmap)
library(enrichR)
library(tibble)
library(grid)
library(conflicted)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(MAST)
library(tidyr)
library(data.table)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("union",  "base")
conflicted::conflict_prefer("setdiff","base")

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Per-dASP chromVAR SummarizedExperiment RDS files (one per dASP cluster)
# Names correspond to dASP cluster indices 1–4
chromvar.paths <- list(
  "1" = "/path/to/chromVAR_Obj/ChromVAR_CellFractionCluster1_fullPeaks.rds",
  "2" = "/path/to/chromVAR_Obj/ChromVAR_CellFractionCluster2_fullPeaks.rds",
  "3" = "/path/to/chromVAR_Obj/ChromVAR_CellFractionCluster3_fullPeaks.rds",
  "4" = "/path/to/chromVAR_Obj/ChromVAR_CellFractionCluster4_fullPeaks.rds"
)

# Per-dASP ATAC Seurat RDS files
atac.paths <- list(
  "1" = "/path/to/SeuratObj/seurat_ATAC_CellFractionCluster1_FullDataset.rds",
  "2" = "/path/to/SeuratObj/seurat_ATAC_CellFractionCluster2_FullDataset.rds",
  "3" = "/path/to/SeuratObj/seurat_ATAC_CellFractionCluster3_FullDataset.rds",
  "4" = "/path/to/SeuratObj/seurat_ATAC_CellFractionCluster4_FullDataset.rds"
)

# Global blast RNA object (IDX/MRD only; used for expression filtering)
path.global.rna <- "/path/to/R_Objects/global.rna.blasts.onlyidxmrd.rds"

# Full cell metadata (all cells including non-blasts)
path.meta.all <- "/path/to/metadata/HTAN_MRD_AllCell_Metadata.txt"

# Blast metadata with dASP cluster assignments (from 01_ script)
path.blast.meta <- "/path/to/Tables/metadata_blast_aggregateCtype.csv"

# dASP cluster assignments CSV (from 01_ script)
path.dasp.clusters <- "/path/to/Tables/Ctype_aggregated_scaled_LFC_k4clustering.csv"

# Output paths
dir.dtf.tables <- "/path/to/ATAC_Processing/Tables/Chromvar_CISBP_ClusterPseudobulk/"
path.dtf.list  <- file.path(dir.dtf.tables,
                              "DTFList.allClustersPB.k4.PTregressed.rds")
path.cc.matrix <- file.path(dir.dtf.tables,
                              "ConsensusDTFlist.ClusterPB.k4.rds")
path.cc.dtflist <- file.path(dir.dtf.tables,
                               "ConsensusDTFlist.allClustersPB.k4.rds")

# DTF filtering thresholds
DTF_PADJ_CUTOFF  <- 0.01     # adjusted p-value cutoff for initial DTF calling
DTF_LFC_CUTOFF   <- 0.25     # |avg_log2FC| threshold
DTF_EXPR_FRAC    <- 0.10     # min fraction of cells expressing TF gene
DTF_CC_LFC_ABS   <- 0.5      # |avg_log2FC| threshold for CC matrix input
DTF_CC_PADJ      <- 0.01     # adjusted p-value threshold for CC matrix input

# Blast lineage assignment (maps Projected.Ctype.Wmethod → Blast_Lineage)
LINEAGE_MAP <- list(
  "Myeloid-like" = c("cDC", "DC-progenitor", "GMP", "MEP", "Monocyte", "pDC"),
  "HSC-like"     = c("HSC/MPP", "LMPP"),
  "NK_T-cell"    = c("NK", "T-cell"),
  "B-lineage"    = c("CLP", "Immature-B", "Mature-B-1", "Mature-B-2",
                      "Pre-B", "Pre-pro-B", "Pro-B")
)

assign_lineage <- function(ctype_vec) {
  result <- rep(NA_character_, length(ctype_vec))
  for (lineage in names(LINEAGE_MAP)) {
    result[ctype_vec %in% LINEAGE_MAP[[lineage]]] <- lineage
  }
  result
}


# =============================================================================
# SECTION 1: Load data
# =============================================================================

cat("Loading chromVAR and ATAC objects...\n")
chromvar.list <- lapply(chromvar.paths, readRDS)
atac.list     <- lapply(atac.paths,    readRDS)

global.rna    <- readRDS(path.global.rna)

# Load and prepare full cell metadata
meta.all <- fread(path.meta.all) %>%
  mutate(Blast_Lineage = assign_lineage(Projected.Ctype.Wmethod))

# Load blast metadata and add dASP cluster assignments
dasp.clusters <- read.csv(path.dasp.clusters, row.names = "X") %>%
  dplyr::select(Row.names, k_4_Myeloid) %>%
  dplyr::rename(PT.number = Row.names, k_cluster = k_4_Myeloid)

metadata_to_add <- read.csv(path.blast.meta, row.names = 1) %>%
  left_join(dasp.clusters, by = "PT.number") %>%
  mutate(Blast_Lineage = assign_lineage(Projected.Ctype.Wmethod)) %>%
  filter(Stage != "RLP")
rownames(metadata_to_add) <- metadata_to_add$patient_bc

global.rna <- AddMetaData(global.rna, metadata = metadata_to_add)

# Split RNA object by dASP cluster (aligned to chromVAR list order)
Idents(global.rna) <- "k_cluster"
cluster.rna.list   <- SplitObject(global.rna, split.by = "k_cluster")
cluster.rna.list   <- cluster.rna.list[order(names(cluster.rna.list))]


# =============================================================================
# SECTION 2: Build per-dASP ATAC+chromVAR objects and filter to blasts
# =============================================================================

cat("Building joint ATAC+chromVAR objects...\n")
joint.objs <- vector("list", length(chromvar.list))

for (i in seq_along(chromvar.list)) {
  z.scores <- chromvar.list[[i]]@assays@data$z

  # Strip chromVAR motif prefix to recover TF gene names
  # Format: <database>_<ID>_<TFname> → keep TFname
  rownames(z.scores) <- sapply(
    rownames(z.scores),
    function(x) tail(strsplit(x, "_")[[1]], 1)
  )

  current.atac <- atac.list[[i]]
  current.atac <- subset(current.atac,
                          cells = colnames(chromvar.list[[i]]))

  # Attach chromVAR z-scores as a new assay
  current.atac[["chromvar"]] <- CreateAssayObject(z.scores)

  # Attach blast cell metadata
  meta.ds <- meta.all[meta.all$patient_bc %in%
                         colnames(chromvar.list[[i]]), ]
  meta.ds <- data.frame(meta.ds)
  rownames(meta.ds) <- meta.ds$patient_bc
  joint.objs[[i]] <- AddMetaData(current.atac, meta.ds)

  # Filter to: blast cells, IDX/MRD stages, non-T/NK arrest states
  Idents(joint.objs[[i]]) <- "Level1.anno"
  joint.objs[[i]] <- subset(joint.objs[[i]], idents = "Blast")

  Idents(joint.objs[[i]]) <- "Stage"
  joint.objs[[i]] <- subset(joint.objs[[i]], idents = c("IDX", "MRD"))

  Idents(joint.objs[[i]]) <- "Projected.Ctype.Wmethod"
  joint.objs[[i]] <- subset(joint.objs[[i]],
                              idents = c("NK", "T-cell"), invert = TRUE)

  # Align RNA object to chromVAR-covered barcodes
  Idents(cluster.rna.list[[i]]) <- "patient_bc"
  cluster.rna.list[[i]] <- subset(
    cluster.rna.list[[i]],
    idents = joint.objs[[i]]$patient_bc
  )

  cat("dASP", i, "—", ncol(joint.objs[[i]]), "blast cells retained\n")
}


# =============================================================================
# SECTION 3: Call DTFs (chromVAR FindMarkers, LR test, patient regressed)
# =============================================================================
# For each dASP, runs a logistic regression with patient ID as a latent
# variable to identify TF motifs with differential accessibility between
# MRD and IDX. TFs are then filtered by gene expression to require the TF
# gene itself to be expressed in >= DTF_EXPR_FRAC of cells in the relevant
# stage (reduces false positives from motif-only signal).

cat("\nCalling DTFs...\n")

expression.all.DTFs <- data.frame()

for (i in seq_along(joint.objs)) {
  cat("Processing dASP", i, "\n")

  current.obj <- joint.objs[[i]]
  Idents(current.obj) <- "Stage"

  current.markers <- FindMarkers(
    current.obj,
    ident.1             = "MRD",
    ident.2             = "IDX",
    logfc.threshold     = DTF_LFC_CUTOFF,
    only.pos            = FALSE,
    test.use            = "LR",
    latent.vars         = "PTnumber",
    max.cells.per.ident = 10000,
    assay               = "chromvar"
  )

  current.markers <- current.markers[
    current.markers$p_val_adj <= DTF_PADJ_CUTOFF, ]
  current.markers <- data.table(current.markers, keep.rownames = TRUE)

  # Expression filter: IDX-enriched TFs must be expressed in IDX cells;
  #                    MRD-enriched TFs must be expressed in MRD cells
  meta.ds      <- meta.all[meta.all$patient_bc %in%
                               current.obj@meta.data$patient_bc, ]
  rna.obj      <- cluster.rna.list[[i]]
  counts.IDX   <- rna.obj@assays$RNA@counts[,
    meta.ds[meta.ds$Stage == "IDX", ]$patient_bc]
  counts.MRD   <- rna.obj@assays$RNA@counts[,
    meta.ds[meta.ds$Stage == "MRD", ]$patient_bc]

  expressed.IDX <- names(which(
    Matrix::rowSums(counts.IDX > 0) >= DTF_EXPR_FRAC * ncol(counts.IDX)
  ))
  expressed.MRD <- names(which(
    Matrix::rowSums(counts.MRD > 0) >= DTF_EXPR_FRAC * ncol(counts.MRD)
  ))

  # IDX-enriched: higher in IDX (negative LFC), must be expressed in IDX
  # MRD-enriched: higher in MRD (positive LFC), must be expressed in MRD
  dtf.IDX <- current.markers[rn %in% expressed.IDX & avg_log2FC < 0, rn]
  dtf.MRD <- current.markers[rn %in% expressed.MRD & avg_log2FC > 0, rn]
  filtered.TFs <- unique(c(dtf.IDX, dtf.MRD))

  current.filtered <- current.markers[rn %in% filtered.TFs] %>%
    dplyr::select(rn, avg_log2FC, p_val, pct.1, pct.2) %>%
    mutate(cluster = i)

  # Save per-dASP table
  write.table(current.filtered, quote = FALSE, sep = "\t", row.names = FALSE,
              file = file.path(dir.dtf.tables,
                                paste0("TrajCluster_", i,
                                       "_PB_EnrichedTFs_MRDvsIDX.txt")))

  expression.all.DTFs <- rbind(expression.all.DTFs, current.filtered)
  cat("  DTFs called:", nrow(current.filtered), "\n")
}

# Jointly adjust p-values across all dASPs
expression.all.DTFs$p_val_adj <- p.adjust(expression.all.DTFs$p_val,
                                            method = "bonferroni")
saveRDS(expression.all.DTFs, file = path.dtf.list)
cat("DTF list saved to:", path.dtf.list, "\n")


# =============================================================================
# SECTION 4: Re-test DTFs across all dASPs and build consensus matrix
# =============================================================================
# Identifies the unique set of significant DTFs across all dASPs, then
# re-runs FindMarkers for each dASP restricted to that TF set so that
# each TF has a score in every dASP (enabling cross-dASP comparison).

cat("\nBuilding consensus DTF matrix...\n")

expression.all.DTFs <- readRDS(path.dtf.list)
all.dtfs <- expression.all.DTFs %>%
  filter(p_val_adj < DTF_CC_PADJ, abs(avg_log2FC) > DTF_CC_LFC_ABS) %>%
  pull(rn) %>% unique()

cat("Unique DTFs for consensus matrix:", length(all.dtfs), "\n")

cc.matrix <- matrix(nrow = length(all.dtfs), ncol = 0)
rownames(cc.matrix) <- all.dtfs

expression.all.unique <- data.frame()

for (i in seq_along(joint.objs)) {
  cat("Scoring dASP", i, "\n")

  current.obj <- joint.objs[[i]]
  Idents(current.obj) <- "Stage"

  # Re-test: no logfc threshold, all unique DTFs as features
  current.markers <- FindMarkers(
    current.obj,
    ident.1             = "MRD",
    ident.2             = "IDX",
    only.pos            = FALSE,
    test.use            = "LR",
    latent.vars         = "PTnumber",
    max.cells.per.ident = 10000,
    min.cells.group     = 1,
    min.cells.feature   = 1,
    min.pct             = 0,
    features            = all.dtfs,
    assay               = "chromvar"
  )

  current.markers <- data.table(current.markers, keep.rownames = TRUE) %>%
    dplyr::select(rn, avg_log2FC, p_val, pct.1, pct.2) %>%
    mutate(cluster = i)

  expression.all.unique <- rbind(expression.all.unique, current.markers)

  # Compute mean z-score difference (MRD − IDX) for each TF in this dASP
  meta.ds    <- meta.all[meta.all$patient_bc %in%
                            current.obj@meta.data$patient_bc, ] %>%
    filter(Level1.anno == "Blast",
           Projected.Ctype.Wmethod %notin% c("NK", "T-cell"),
           Stage != "RLP")

  z.mat      <- chromvar.list[[i]]@assays@data$z[, meta.ds$patient_bc]
  rownames(z.mat) <- sapply(rownames(z.mat),
                              function(x) tail(strsplit(x, "_")[[1]], 1))

  mean.z <- sapply(c("MRD", "IDX"), function(stage) {
    Matrix::rowMeans(z.mat[, meta.ds[meta.ds$Stage == stage, ]$patient_bc])
  })

  delta.z <- (mean.z[, "MRD"] - mean.z[, "IDX"])[all.dtfs, drop = FALSE]

  cc.matrix <- cbind(cc.matrix, delta.z)
}

colnames(cc.matrix) <- as.character(seq_along(joint.objs))

saveRDS(cc.matrix, file = path.cc.matrix)
saveRDS(split(expression.all.unique, expression.all.unique$cluster),
        file = path.cc.dtflist)

cat("Consensus matrix saved to:", path.cc.matrix, "\n")
cat("Dimensions:", nrow(cc.matrix), "DTFs ×", ncol(cc.matrix), "dASPs\n")
