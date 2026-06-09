# =============================================================================
# Joint Assay (Multiome) Processing Pipeline
# =============================================================================
#
# Description:
#   End-to-end processing of paired scRNA-seq + scATAC-seq (10x Multiome)
#   data across multiple samples. The pipeline covers:
#     1. Building a sample-to-path metadata table from a batch manifest
#     2. Per-sample QC, doublet removal, and joint Seurat object creation
#     3. Merging all samples into a combined RNA Seurat object with
#        regression of cell-cycle, mitochondrial, and heat-shock variation
#     4. Annotating cells with ATAC QC metrics and clinical metadata
#     5. QC visualizations (FRiP and UMI boxplots per sample)
#
# Inputs:
#   - CellRanger output directories (one per sample)
#   - scATACpro output directories (one per sample)
#   - Genome annotation RDS (GRanges object for hg38)
#   - Cell cycle gene list (Regev lab format: S genes rows 1-43, G2M rows 44-97)
#   - Heat shock gene list (tab-separated, column "Approved symbol")
#   - Clinical metadata CSV with columns: USI, Subtype, MRD%, MRD_Group
#   - Batch manifest CSV (see BATCH MANIFEST section below)
#
# Outputs:
#   - metadata/AllDatapath.rds          : sample-to-path lookup table
#   - SeuratObj/BySample/JointObj/      : per-sample Seurat objects
#   - SeuratObj/seurat_combined_RNA.rds : merged RNA Seurat object
#   - metadata/combined_meta.txt        : merged cell-level metadata
#   - Figures/                          : QC and UMAP plots (.eps)
#
# Dependencies:
#   Seurat, Signac, EnsDb.Hsapiens.v86, data.table, RColorBrewer,
#   R.utils, ggplot2, tidyverse, gdata, reshape2, Matrix.utils
#
#   Local utility scripts (must be in the working directory or on R path):
#     scDataAnalysis_Utilities.R
#     dsAnalysis_utilities.R
#     DoubletsFinder.R  — provides FindDoublets()
#
# Usage:
#   1. Fill in the "USER CONFIGURATION" block below.
#   2. Provide a batch manifest CSV (see template comments).
#   3. Run sections sequentially or submit as a batch job.
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(tidyverse)
library(gdata)
library(reshape2)

source("scDataAnalysis_Utilities.R")
source("dsAnalysis_utilities.R")
source("DoubletsFinder.R")      # provides FindDoublets()

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────
# Edit all paths and parameters in this block before running.

# Root directories for upstream pipeline outputs
dir.cellranger  <- "/path/to/CellrangerResult/"   # trailing slash required
dir.scATACpro   <- "/path/to/scATACpro/"

# Project-level output directories (create these before running)
dir.output.base    <- "/path/to/project/Scripts/"
dir.seurat.samples <- file.path(dir.output.base, "SeuratObj/BySample/JointObj/")
dir.figures        <- file.path(dir.output.base, "Figures/SingleSample/")
dir.metadata       <- file.path(dir.output.base, "metadata/")

# Reference files
path.annotation    <- "/path/to/GRanges_hg38.rds"          # hg38 gene annotation
path.cell.cycle    <- "/path/to/regev_lab_cell_cycle_genes.txt"
path.heat.shock    <- "/path/to/heat_shock_geneList.txt"    # column: "Approved symbol"
path.clinical.info <- "/path/to/ClinicalInfo.csv"           # columns: USI, Subtype, MRD%, MRD_Group

# Batch manifest CSV — one row per sample, columns:
#   sample_name  : e.g. "PAPBCW_IDX"
#   batch_date   : e.g. "Feb24_2022"  (subfolder name under dir.cellranger)
# Example row: PAPBCW_IDX,Feb24_2022
path.batch.manifest <- "/path/to/batch_manifest.csv"

# Per-cell QC thresholds
QC_MAX_ATAC  <- 100000   # max ATAC UMI per cell
QC_MAX_RNA   <- 50000    # max RNA UMI per cell
QC_MIN_ATAC  <- 1000     # min ATAC UMI per cell
QC_MIN_RNA   <- 1000     # min RNA UMI per cell
QC_MAX_MT    <- 10       # max mitochondrial read % per cell

# Dimensionality reduction parameters
N_PCA_DIMS_PERSAMPLE <- 50    # PCs used in per-sample neighbor graph
N_LSI_DIMS           <- 2:40  # LSI dims for ATAC (excludes dim 1 by convention)
N_PCA_DIMS_COMBINED  <- 30    # PCs used in combined object
N_HVG_PERSAMPLE      <- 2000  # highly variable genes per sample
N_HVG_COMBINED       <- 4000  # highly variable genes in combined object
MIN_CELL_EXPR_FRAC   <- 0.01  # min fraction of cells expressing a gene to keep it
N_PCA_MAX            <- 100   # max PCs computed in combined object

# Clustering resolution
CLUSTER_RES <- 0.2


# ── SECTION 1: Build sample-to-path metadata table ───────────────────────────
#
# Reads a batch manifest CSV that maps each sample to its CellRanger batch
# subdirectory, then constructs the full input paths for CellRanger and
# scATACpro outputs.
#
# The manifest replaces the hard-coded per-batch `if` blocks in the original
# script. To add new samples, simply append rows to the CSV — no code changes
# needed.

batch.manifest <- fread(path.batch.manifest,
                        col.names = c("sample_name", "batch_date"))

data.path <- batch.manifest[, .(
  sample.name    = sample_name,
  patient.ID     = sapply(strsplit(sample_name, "_"), `[[`, 1),
  Stage          = sapply(strsplit(sample_name, "_"), `[[`, 2),
  cellranger.path = file.path(dir.cellranger, batch_date, sample_name),
  scATACpro.path  = file.path(dir.scATACpro,  sample_name)
)]

saveRDS(data.path,
        file = file.path(dir.metadata, "AllDatapath.rds"))

cat("Sample manifest saved:", nrow(data.path), "samples\n")


# ── SECTION 2: Per-sample QC and joint Seurat object creation ─────────────────
#
# For each sample:
#   1. Loads CellRanger (RNA) and scATACpro (ATAC) count matrices.
#   2. Identifies high-quality barcodes shared between both assays after
#      doublet removal (FindDoublets wraps DoubletFinder on the RNA assay).
#   3. Creates a multimodal Seurat object, computes QC metrics, and applies
#      per-cell filters.
#   4. Runs standard RNA (PCA) and ATAC (LSI) dimensionality reduction,
#      builds a weighted nearest-neighbor (WNN) graph, and generates a
#      joint UMAP.
#   5. Saves the per-sample object and two diagnostic plots.

data.path  <- readRDS(file.path(dir.metadata, "AllDatapath.rds"))
annotation <- readRDS(path.annotation)

for (sample0 in data.path$sample.name) {
  cat("\n--- Processing:", sample0, "---\n")

  sp <- data.path[data.path$sample.name == sample0, ]  # single-sample metadata row

  # -- Load count matrices --------------------------------------------------
  counts   <- Read10X(file.path(sp$cellranger.path, "outs/filtered_feature_bc_matrix/"))
  RNA.bc   <- FindDoublets(counts)          # returns barcodes after doublet removal
  atac.count <- readRDS(file.path(sp$scATACpro.path,
                                  "output/filtered_matrix/MACS2/FILTER/matrix.rds"))

  # Keep only barcodes present in both modalities
  filtered.bc          <- intersect(RNA.bc, colnames(atac.count))
  RNA.count.filter.mtx <- counts[["Gene Expression"]][, filtered.bc]
  ATAC.count.filter.mtx <- atac.count[, filtered.bc]

  # -- Build multimodal Seurat object ----------------------------------------
  seurat.obj <- CreateSeuratObject(RNA.count.filter.mtx, assay = "RNA")
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")

  seurat.obj[["ATAC"]] <- CreateChromatinAssay(
    counts    = ATAC.count.filter.mtx,
    sep       = c(":", "-"),
    fragments = file.path(sp$scATACpro.path, "output/summary",
                          paste0(sample0, ".fragments.tsv.gz")),
    annotation = annotation
  )

  seurat.obj <- RenameCells(seurat.obj, add.cell.id = sample0)
  seurat.obj$patientID <- sp$patient.ID
  seurat.obj$Stage     <- sp$Stage

  # -- ATAC QC metric --------------------------------------------------------
  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- TSSEnrichment(seurat.obj)

  # -- Per-cell filtering ----------------------------------------------------
  seurat.obj <- subset(
    seurat.obj,
    subset = nCount_ATAC  < QC_MAX_ATAC  &
             nCount_RNA   < QC_MAX_RNA   &
             nCount_ATAC  > QC_MIN_ATAC  &
             nCount_RNA   > QC_MIN_RNA   &
             percent.mt   < QC_MAX_MT
  )

  # -- RNA dimensionality reduction ------------------------------------------
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj,
                                     selection.method = "vst",
                                     nfeatures = N_HVG_PERSAMPLE)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)

  # -- ATAC dimensionality reduction (LSI) -----------------------------------
  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = 5)
  seurat.obj <- RunTFIDF(seurat.obj)
  seurat.obj <- RunSVD(seurat.obj)

  # -- Weighted nearest-neighbor joint embedding -----------------------------
  seurat.obj <- FindMultiModalNeighbors(
    object         = seurat.obj,
    reduction.list = list("pca", "lsi"),
    dims.list      = list(1:N_PCA_DIMS_PERSAMPLE, N_LSI_DIMS),
    modality.weight.name = "RNA.weight",
    verbose        = TRUE
  )
  seurat.obj <- RunUMAP(seurat.obj, nn.name = "weighted.nn",
                        assay = "RNA", verbose = TRUE)

  # -- Clustering (RNA-based) -----------------------------------------------
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:N_PCA_DIMS_PERSAMPLE)
  seurat.obj <- FindClusters(seurat.obj,  resolution = 0.5)

  # -- Diagnostic plots ------------------------------------------------------
  p_umap <- DimPlot(seurat.obj, label = TRUE, repel = TRUE, reduction = "umap") +
    NoLegend() + ggtitle(sample0)

  p_feat <- FeaturePlot(seurat.obj,
                        features = c("CD19", "MME", "CD34", "MS4A1",
                                     "GNLY", "CD3E", "LYZ", "CD36"))

  ggsave(p_umap, device = "eps", width = 9, height = 7,
         filename = file.path(dir.figures,
                              paste0(sample0, "_UMAP_doubletRemoved.eps")))

  ggsave(p_feat, device = "eps", width = 9, height = 9,
         filename = file.path(dir.figures,
                              paste0(sample0, "_FeaturePlot_doubletRemoved.eps")))

  # -- Save ------------------------------------------------------------------
  saveRDS(seurat.obj,
          file = file.path(dir.seurat.samples,
                           paste0(sample0, "_JointSeuratObj_doubletRemoved.rds")))

  cat("Finished:", sample0, "\n")
}


# ── SECTION 3: Merge per-sample objects into combined RNA Seurat object ────────
#
# Loads all per-sample joint objects, concatenates their RNA count matrices,
# then runs a full RNA analysis pipeline on the merged object including:
#   - Cell-cycle scoring and regression
#   - Heat-shock module score regression
#   - PCA, UMAP, clustering

seurat.files <- dir(dir.seurat.samples, pattern = "doubletRemoved.rds",
                    full.names = TRUE)

# -- Concatenate RNA matrices and metadata -----------------------------------
cat("Merging", length(seurat.files), "sample objects...\n")

seurat.tmp <- readRDS(seurat.files[1])
seurat.tmp$sample.name <- paste0(unique(seurat.tmp$patientID), "_",
                                 unique(seurat.tmp$Stage))
rna.mtx <- seurat.tmp[["RNA"]]@counts
meta    <- seurat.tmp@meta.data

for (f in seurat.files[-1]) {
  seurat.tmp <- readRDS(f)
  seurat.tmp$sample.name <- paste0(unique(seurat.tmp$patientID), "_",
                                   unique(seurat.tmp$Stage))
  meta    <- rbind(meta, seurat.tmp@meta.data)
  rna.mtx <- t(Matrix.utils::rBind.fill(t(rna.mtx),
                                         t(seurat.tmp[["RNA"]]@counts),
                                         fill = 0))
}

# -- Build combined Seurat object --------------------------------------------
seurat.rna <- CreateSeuratObject(rna.mtx)
seurat.rna <- NormalizeData(seurat.rna)
seurat.rna <- AddMetaData(seurat.rna, metadata = meta)

# percent.mt was stored as a fraction (0–1) in per-sample objects;
# Seurat PercentageFeatureSet returns 0–100, so divide here if needed
# (adjust this line if your per-sample objects already store 0–100)
seurat.rna$percent.mt <- seurat.rna$percent.mt / 100

# -- Highly variable gene selection with expression-frequency filter ---------
seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = N_HVG_COMBINED)

cell.expr.frac <- rowSums(seurat.rna@assays$RNA@counts >= 1) /
                  ncol(seurat.rna@assays$RNA@counts)
low.expr.genes <- names(cell.expr.frac)[cell.expr.frac < MIN_CELL_EXPR_FRAC]
VariableFeatures(seurat.rna) <- setdiff(VariableFeatures(seurat.rna),
                                         low.expr.genes)
cat("Variable features after expression filter:",
    length(VariableFeatures(seurat.rna)), "\n")

# -- Cell-cycle scoring -------------------------------------------------------
# Gene list format: rows 1–43 = S-phase genes, rows 44–97 = G2M genes
cycle.genes  <- fread(path.cell.cycle, header = FALSE)$V1
s.genes      <- cycle.genes[1:43]
g2m.genes    <- cycle.genes[44:97]

seurat.rna <- CellCycleScoring(seurat.rna,
                                s.features   = s.genes,
                                g2m.features = g2m.genes,
                                set.ident    = TRUE)

# -- Heat-shock module score --------------------------------------------------
heat.shock.genes <- fread(path.heat.shock)$`Approved symbol`
seurat.rna <- AddModuleScore(seurat.rna,
                              features = list(heat.shock.genes),
                              name     = "HeatShock.Score")

# -- Scale and regress out confounders ----------------------------------------
seurat.rna <- ScaleData(
  seurat.rna,
  features      = VariableFeatures(seurat.rna),
  vars.to.regress = c("S.Score", "G2M.Score", "percent.mt",
                      "HeatShock.Score1", "nCount_RNA")
)

# -- PCA, UMAP, clustering ----------------------------------------------------
seurat.rna <- RunPCA(seurat.rna, npcs = N_PCA_MAX, verbose = FALSE,
                     features = VariableFeatures(seurat.rna))
seurat.rna <- RunUMAP(seurat.rna, dims = 1:N_PCA_DIMS_COMBINED)
seurat.rna <- FindNeighbors(seurat.rna, dims = 1:N_PCA_DIMS_COMBINED,
                             reduction = "pca")
seurat.rna <- FindClusters(seurat.rna, resolution = CLUSTER_RES)

# -- Save metadata and object -------------------------------------------------
write.table(seurat.rna@meta.data,
            file  = file.path(dir.metadata, "combined_meta.txt"),
            quote = FALSE, sep = "\t")

saveRDS(seurat.rna,
        file = file.path(dir.output.base,
                         "SeuratObj/seurat_combined_RNA.rds"))

cat("Combined RNA object saved.\n")


# ── SECTION 4: Add ATAC QC metrics and clinical metadata to per-sample objects ─
#
# For each per-sample object, joins:
#   - Per-barcode ATAC QC metrics from scATACpro (FRiP, TSS fraction, etc.)
#   - Clinical variables from the study metadata (cancer subtype, MRD values)
# Saves updated objects and collects FRiP / UMI vectors for QC plotting.

clinical.info <- fread(path.clinical.info)

seurat.files  <- dir(dir.seurat.samples, pattern = "doubletRemoved.rds",
                     full.names = TRUE)
UMI.list  <- list()
FriP.list <- list()

for (f in seurat.files) {
  seurat.obj  <- readRDS(f)
  sample.name <- paste0(unique(seurat.obj$patientID), "_",
                        unique(seurat.obj$Stage))
  seurat.obj$sample.name <- sample.name

  meta0 <- seurat.obj@meta.data

  # Attach clinical variables via patient ID
  pid <- unique(seurat.obj$patientID)
  meta0$Subtype  <- clinical.info[clinical.info$USI == pid, ]$Subtype
  meta0$MRD      <- clinical.info[clinical.info$USI == pid, ]$`MRD%`
  meta0$MRDGroup <- clinical.info[clinical.info$USI == pid, ]$MRD_Group

  # Load scATACpro per-barcode QC table and join on barcode
  atac.qc.file <- file.path(dir.scATACpro, sample.name, "output/summary",
                             paste0(sample.name, ".MACS2.qc_per_barcode.txt"))
  scATAC.meta <- fread(atac.qc.file)
  scATAC.meta$bc <- paste0(sample.name, "_", scATAC.meta$bc)
  setkey(scATAC.meta, bc)

  meta0 <- data.table(meta0, keep.rownames = TRUE)
  setkey(meta0, rn)

  meta0[, FriP             := scATAC.meta[J(meta0$rn)]$frac_peak]
  meta0[, FractionTSS      := scATAC.meta[J(meta0$rn)]$frac_tss]
  meta0[, FractionPromoter := scATAC.meta[J(meta0$rn)]$frac_promoter]
  meta0[, FractionEnhancer := scATAC.meta[J(meta0$rn)]$frac_enhancer]

  UMI.list[[sample.name]]  <- meta0$nCount_RNA
  FriP.list[[sample.name]] <- meta0$FriP

  write.table(meta0,
              file = file.path(dir.metadata,
                               paste0(sample.name, "_metadata_withATAC.txt")))

  rownames(meta0) <- meta0$rn
  seurat.obj <- AddMetaData(seurat.obj, metadata = meta0)

  saveRDS(seurat.obj,
          file = file.path(dir.seurat.samples,
                           paste0(sample.name,
                                  "_JointSeuratObj_doubletRemoved_withMetadata.rds")))
}

cat("Metadata annotation complete for", length(seurat.files), "samples.\n")


# ── SECTION 5: QC boxplots across all samples ────────────────────────────────
#
# Generates per-sample boxplots of RNA UMI counts and ATAC FRiP scores.
# Color palette is automatically scaled to the number of samples.

n.samples  <- length(FriP.list)
myColors   <- colorRampPalette(brewer.pal(12, "Paired"))(n.samples)
names(myColors) <- names(FriP.list)

# Helper: melt a named list of numeric vectors into a two-column data.frame
melt_list <- function(lst) {
  max.len <- max(lengths(lst))
  padded  <- lapply(lst, function(x) c(x, rep(NA, max.len - length(x))))
  melt(as.data.frame(padded), na.rm = TRUE)
}

df.frip.melt <- melt_list(FriP.list)
df.umi.melt  <- melt_list(UMI.list)

p.frip <- ggplot(df.frip.melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +
  ylim(0, 1) +
  labs(title = "ATAC FRiP per sample", x = NULL, y = "FRiP") +
  theme(legend.position = "none") +
  guides(x = guide_axis(angle = 90)) +
  scale_fill_manual(values = myColors)

p.umi <- ggplot(df.umi.melt, aes(x = variable, y = value, fill = variable)) +
  geom_boxplot(alpha = 1, outlier.shape = NA) +
  ylim(0, 15000) +
  labs(title = "RNA UMI per sample", x = NULL, y = "nCount_RNA") +
  theme(legend.position = "none") +
  guides(x = guide_axis(angle = 90)) +
  scale_fill_manual(values = myColors)

ggsave(p.frip, device = "eps", width = 12, height = 14,
       filename = file.path(dir.output.base, "Figures/QC_scATAC_FRiP_allSamples.eps"))

ggsave(p.umi, device = "eps", width = 12, height = 14,
       filename = file.path(dir.output.base, "Figures/QC_scRNA_UMI_allSamples.eps"))

cat("QC plots saved.\n")

