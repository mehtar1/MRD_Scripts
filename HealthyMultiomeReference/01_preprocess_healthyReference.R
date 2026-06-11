# =============================================================================
# Healthy Donor Multiome Reference — Per-Sample Processing Pipeline
# =============================================================================
#
# Description:
#   Processes healthy donor 10x Multiome (scRNA-seq + scATAC-seq) samples
#   into per-sample Seurat objects to build a healthy reference atlas.
#   This pipeline is structurally identical to the leukemia sample preprocessing
#   in 01_jointAssay_processing.R, with three differences specific to healthy
#   donors:
#     - Metadata columns are DonorID and SortingStrategy (not patientID/Stage)
#     - Mitochondrial % threshold is relaxed to 20% (vs. 10% for leukemia)
#     - No FeaturePlot diagnostic (marker genes less informative in healthy cells)
#
#   Steps per sample:
#     1. Load CellRanger (RNA) and scATACpro (ATAC) count matrices
#     2. Intersect high-quality barcodes after doublet removal
#     3. Build multimodal Seurat object with RNA + ATAC assays
#     4. Apply per-cell QC filters
#     5. Normalize and reduce RNA (PCA) and ATAC (LSI) independently
#     6. Build weighted nearest-neighbor (WNN) joint embedding and UMAP
#     7. Cluster on RNA and save per-sample object + UMAP plot
#
# Inputs (set in USER CONFIGURATION):
#   - CellRanger output directories (one per sample)
#   - scATACpro output directories (one per sample)
#   - Genome annotation RDS (GRanges object, hg38)
#   - Batch manifest CSV (see template comment below)
#
# Outputs:
#   - metadata/HD_DataPath.rds           : sample-to-path lookup table
#   - SeuratObj/HDMultiome/JointObj/     : per-sample Seurat objects
#   - Figures/HDMultiome/                : per-sample UMAP plots (.eps)
#
# Batch manifest CSV format (no header, two columns):
#   sample_name,batch_date
#   2986_Live,Mar23_2023
#   2986_Enriched,Mar23_2023
#   ...
#
# Dependencies:
#   Seurat, Signac, EnsDb.Hsapiens.v86, data.table, RColorBrewer,
#   R.utils, ggplot2, tidyverse, gdata
#
#   Local utility scripts (must be on the R source path):
#     scDataAnalysis_Utilities.R
#     dsAnalysis_utilities.R
#     DoubletsFinder.R   — provides FindDoublets()
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

source("scDataAnalysis_Utilities.R")
source("dsAnalysis_utilities.R")
source("DoubletsFinder.R")

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Root directories for upstream pipeline outputs
dir.cellranger <- "/path/to/CellrangerResult/"
dir.scATACpro  <- "/path/to/scATACpro_HD_multiome/"

# Project-level output directories (create before running)
dir.output.base    <- "/path/to/project/Scripts/"
dir.seurat.samples <- file.path(dir.output.base, "SeuratObj/HDMultiome/JointObj/")
dir.figures        <- file.path(dir.output.base, "Figures/HDMultiome/")
dir.metadata       <- file.path(dir.output.base, "metadata/")

# Reference files
path.annotation    <- "/path/to/GRanges_hg38.rds"   # hg38 gene annotation

# Batch manifest CSV — one row per sample, columns: sample_name, batch_date
# Example row: 2986_Live,Mar23_2023
path.batch.manifest <- "/path/to/batch_manifest_HD.csv"

# Per-cell QC thresholds
# NOTE: percent.mt threshold is relaxed to 20% for healthy donors (vs. 10%
# for leukemia samples) because healthy cells have higher baseline mt expression.
QC_MAX_ATAC <- 100000
QC_MAX_RNA  <- 50000
QC_MIN_ATAC <- 1000
QC_MIN_RNA  <- 1000
QC_MAX_MT   <- 20

# Dimensionality reduction parameters
N_PCA_DIMS   <- 50      # PCs for RNA neighbor graph and WNN
N_LSI_DIMS   <- 2:40    # LSI dims for ATAC (dim 1 excluded: correlated with depth)
N_HVG        <- 2000    # highly variable genes per sample
CLUSTER_RES  <- 0.5     # clustering resolution


# =============================================================================
# SECTION 1: Build sample-to-path metadata table
# =============================================================================
# Reads a batch manifest CSV that maps each sample to its CellRanger batch
# subdirectory. The sample name encodes DonorID and SortingStrategy as
# <DonorID>_<SortingStrategy> (e.g. "2986_Live", "2986_Enriched").
#
# To add new samples, append rows to the CSV — no code changes required.

batch.manifest <- fread(path.batch.manifest,
                        col.names = c("sample_name", "batch_date"))

data.path <- batch.manifest[, .(
  sample.name      = sample_name,
  Donor.ID         = sapply(strsplit(sample_name, "_"), `[[`, 1),
  SortingStrategy  = sapply(strsplit(sample_name, "_"), `[[`, 2),
  cellranger.path  = file.path(dir.cellranger, batch_date, sample_name),
  scATACpro.path   = file.path(dir.scATACpro,  sample_name)
)]

saveRDS(data.path, file = file.path(dir.metadata, "HD_DataPath.rds"))
cat("Sample manifest saved:", nrow(data.path), "samples\n")


# =============================================================================
# SECTION 2: Per-sample QC and joint Seurat object creation
# =============================================================================
# For each healthy donor sample:
#   1. Loads CellRanger (RNA) and scATACpro (ATAC) count matrices.
#   2. Identifies high-quality barcodes shared between both assays after
#      doublet removal via FindDoublets (wraps DoubletFinder on RNA assay).
#   3. Builds a multimodal Seurat object, computes QC metrics, and filters.
#   4. Runs standard RNA (PCA) and ATAC (LSI) dimensionality reduction,
#      builds a WNN joint graph, and generates a UMAP.
#   5. Saves the per-sample object and a UMAP diagnostic plot.

data.path  <- readRDS(file.path(dir.metadata, "HD_DataPath.rds"))
annotation <- readRDS(path.annotation)

for (sample0 in data.path$sample.name) {
  cat("\n--- Processing:", sample0, "---\n")

  sp <- data.path[data.path$sample.name == sample0, ]

  # -- Load count matrices --------------------------------------------------
  counts     <- Read10X(file.path(sp$cellranger.path,
                                   "outs/filtered_feature_bc_matrix/"))
  RNA.bc     <- FindDoublets(counts)
  atac.count <- readRDS(file.path(sp$scATACpro.path,
                                   "output/filtered_matrix/MACS2/FILTER/matrix.rds"))

  # Keep only barcodes present in both modalities
  filtered.bc           <- intersect(RNA.bc, colnames(atac.count))
  RNA.count.filter.mtx  <- counts[["Gene Expression"]][, filtered.bc]
  ATAC.count.filter.mtx <- atac.count[, filtered.bc]

  # -- Build multimodal Seurat object ---------------------------------------
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
  seurat.obj$DonorID          <- sp$Donor.ID
  seurat.obj$SortingStrategy  <- sp$SortingStrategy

  # -- ATAC QC metric -------------------------------------------------------
  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- TSSEnrichment(seurat.obj)

  # -- Per-cell filtering ---------------------------------------------------
  seurat.obj <- subset(
    seurat.obj,
    subset = nCount_ATAC  < QC_MAX_ATAC &
             nCount_RNA   < QC_MAX_RNA  &
             nCount_ATAC  > QC_MIN_ATAC &
             nCount_RNA   > QC_MIN_RNA  &
             percent.mt   < QC_MAX_MT
  )

  # -- RNA dimensionality reduction -----------------------------------------
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj,
                                      selection.method = "vst",
                                      nfeatures        = N_HVG)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)

  # -- ATAC dimensionality reduction (LSI) ----------------------------------
  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = 5)
  seurat.obj <- RunTFIDF(seurat.obj)
  seurat.obj <- RunSVD(seurat.obj)

  # -- Weighted nearest-neighbor joint embedding ----------------------------
  seurat.obj <- FindMultiModalNeighbors(
    object               = seurat.obj,
    reduction.list       = list("pca", "lsi"),
    dims.list            = list(1:N_PCA_DIMS, N_LSI_DIMS),
    modality.weight.name = "RNA.weight",
    verbose              = TRUE
  )
  seurat.obj <- RunUMAP(seurat.obj,
                         nn.name = "weighted.nn",
                         assay   = "RNA",
                         verbose = TRUE)

  # -- Clustering (RNA-based) -----------------------------------------------
  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:N_PCA_DIMS)
  seurat.obj <- FindClusters(seurat.obj,  resolution = CLUSTER_RES)

  # -- Diagnostic UMAP plot -------------------------------------------------
  p_umap <- DimPlot(seurat.obj, label = TRUE, repel = TRUE,
                     reduction = "umap") +
    NoLegend() +
    ggtitle(sample0)

  ggsave(p_umap, device = "eps", width = 9, height = 7,
         filename = file.path(dir.figures,
                               paste0(sample0, "_UMAP_doubletRemoved.eps")))

  # -- Save -----------------------------------------------------------------
  saveRDS(seurat.obj,
          file = file.path(dir.seurat.samples,
                            paste0(sample0,
                                   "_JointSeuratObj_doubletRemoved.rds")))

  cat("Finished:", sample0, "\n")
}

cat("\nAll samples processed.\n")
