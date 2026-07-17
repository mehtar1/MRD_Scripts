# =============================================================================
# Healthy Donor Multiome Reference — Full Construction Pipeline
# =============================================================================
#
# Description:
#   Processes healthy donor 10x Multiome (scRNA-seq + scATAC-seq) samples
#   into a merged, jointly analyzed healthy reference atlas.
#
#   Pipeline overview:
#     SECTION 1 : Build sample-to-path metadata table
#     SECTION 2 : Per-sample QC and joint Seurat object creation
#     SECTION 3 : Merge all per-sample objects into one combined object
#     SECTION 4 : RNA processing — Gini-iterative feature selection + PCA
#     SECTION 5 : ATAC processing — DAP-based feature selection + LSI
#     SECTION 6 : WNN joint embedding (multi-k neighbors + multi-UMAP)
#     SECTION 7 : Label transfer from old RNA reference
#     SECTION 8 : Cell-cycle & heat-shock score regression
#     SECTION 9 : B-trajectory downsampling and final joint analysis
#     SECTION 10: Manual cell-type annotation
#     SECTION 11: Pseudotime trajectory inference with Slingshot
#
# Inputs (set in USER CONFIGURATION):
#   - CellRanger output directories (one per sample)
#   - scATACpro output directories (one per sample)
#   - Combined ATAC fragment file (all donors merged)
#   - Genome annotation RDS (GRanges object, hg38)
#   - Batch manifest CSV (sample_name, batch_date)
#   - Old scRNA-seq reference RDS (for label transfer)
#   - Cell-cycle gene list (Regev lab format, 97 genes total)
#   - Heat-shock gene list
#
# Outputs:
#   - metadata/HD_DataPath.rds
#   - SeuratObj/HDMultiome/JointObj/          per-sample objects
#   - SeuratObj/HDMultiome/Combined/          merged + jointly analyzed objects
#   - SeuratObj/HDMultiome/Downsample_Btraj/  downsampled final reference
#   - Trajectory_Objects/                     Slingshot trajectory object
#   - Input4VisCello/                         VisCello session (eset + clist)
#   - Figures/HDMultiome/                     per-sample and combined UMAP plots
#
# Dependencies:
#   Seurat, Signac, Harmony, EnsDb.Hsapiens.v86, data.table, RColorBrewer,
#   R.utils, ggplot2, tidyverse, gdata, Matrix.utils, slingshot,
#   SingleCellExperiment, Biobase
#
#   Local utility scripts:
#     scDataAnalysis_Utilities.R  — provides ifg_select() (Gini-based feature selection)
#     dsAnalysis_utilities.R
#     DoubletsFinder.R            — provides FindDoublets()
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
library(Matrix.utils)
library(slingshot)
library(SingleCellExperiment)
library(Biobase)
source("scDataAnalysis_Utilities.R")
source("dsAnalysis_utilities.R")
source("DoubletsFinder.R")
`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────
dir.cellranger <- "/path/to/CellrangerResult/"
dir.scATACpro  <- "/path/to/scATACpro_HD_multiome/"

dir.output.base    <- "/path/to/project/Scripts/"
dir.seurat.samples <- file.path(dir.output.base, "SeuratObj/HDMultiome/JointObj/")
dir.seurat.combined <- file.path(dir.output.base, "SeuratObj/HDMultiome/Combined/")
dir.seurat.downsample <- file.path(dir.output.base, "SeuratObj/HDMultiome/Downsample_Btraj/")
dir.trajectory     <- file.path(dir.output.base, "Trajectory_Objects/")
dir.viscello       <- file.path(dir.output.base, "Input4VisCello/")
dir.figures        <- file.path(dir.output.base, "Figures/HDMultiome/")
dir.metadata       <- file.path(dir.output.base, "metadata/")

path.annotation    <- "/path/to/GRanges_hg38.rds"
path.batch.manifest <- "/path/to/batch_manifest_HD.csv"

# Combined fragment file for all donors (produced by scATACpro merge step)
path.fragments.combined <- file.path(dir.scATACpro,
                                     "Fragments/Combine_7Donors/7Donors_combined.fragments.filtered.sorted.bed.gz")

# Old scRNA-seq reference (for label transfer)
path.old.rna.ref <- "/path/to/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC.rds"

# Signature gene lists
path.cell.cycle.genes <- "/path/to/regev_lab_cell_cycle_genes.txt"  # 97 genes: 1:43 = S, 44:97 = G2M
path.heat.shock.genes <- "/path/to/heat_shock_geneList.txt"         # column: "Approved symbol"

# Per-cell QC thresholds (Sections 1–2 only)
QC_MAX_ATAC <- 100000
QC_MAX_RNA  <- 50000
QC_MIN_ATAC <- 1000
QC_MIN_RNA  <- 1000
QC_MAX_MT   <- 20          # relaxed to 20% for healthy donors

# Dimensionality reduction parameters
N_PCA_DIMS   <- 50
N_LSI_DIMS   <- 2:40
N_HVG        <- 2000
CLUSTER_RES  <- 0.5

# Cell types in the reference atlas
CELLTYPE_LEVELS <- c("HSC/MPP", "LMPP", "CLP", "Pre-pro-B", "Pro-B", "Pre-B",
                     "Immature-B", "Mature-B-1", "Mature-B-2",
                     "GMP", "MEP", "DC-progenitor", "Monocyte", "cDC", "pDC",
                     "T-cell", "NK")
B_LINEAGE_TYPES <- c("HSC/MPP", "LMPP", "CLP", "Pre-pro-B", "Pro-B", "Pre-B",
                     "Immature-B", "Mature-B-1", "Mature-B-2")

# Color palette (17 cell types)
CELLTYPE_COLORS <- c("#DB8312", "#36A438", "#39A6D6", "#AA7AAD", "#B41273",
                     "#22ACA3", "#CC79AE", "#DA5D8D", "#27A241", "#58C66A",
                     "#1DA774", "#E26E63", "#B0A21B", "#4791C4", "#86AA35",
                     "#CA9611", "#22ACA3")


# =============================================================================
# SECTION 1: Build sample-to-path metadata table
# =============================================================================
# Sample name encodes <DonorID>_<SortingStrategy> (e.g. "2986_Live").
# Append rows to the CSV to add new samples without code changes.

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
# For each sample: load matrices → doublet removal → build multimodal Seurat
# object → QC filtering → RNA PCA → ATAC LSI → WNN UMAP → save.

data.path  <- readRDS(file.path(dir.metadata, "HD_DataPath.rds"))
annotation <- readRDS(path.annotation)

for (sample0 in data.path$sample.name) {
  cat("\n--- Processing:", sample0, "---\n")
  sp <- data.path[data.path$sample.name == sample0, ]

  counts     <- Read10X(file.path(sp$cellranger.path,
                                   "outs/filtered_feature_bc_matrix/"))
  RNA.bc     <- FindDoublets(counts)
  atac.count <- readRDS(file.path(sp$scATACpro.path,
                                   "output/filtered_matrix/MACS2/FILTER/matrix.rds"))

  filtered.bc           <- intersect(RNA.bc, colnames(atac.count))
  RNA.count.filter.mtx  <- counts[["Gene Expression"]][, filtered.bc]
  ATAC.count.filter.mtx <- atac.count[, filtered.bc]

  seurat.obj <- CreateSeuratObject(RNA.count.filter.mtx, assay = "RNA")
  seurat.obj[["percent.mt"]] <- PercentageFeatureSet(seurat.obj, pattern = "^MT-")
  seurat.obj[["ATAC"]] <- CreateChromatinAssay(
    counts     = ATAC.count.filter.mtx,
    sep        = c(":", "-"),
    fragments  = file.path(sp$scATACpro.path, "output/summary",
                            paste0(sample0, ".fragments.tsv.gz")),
    annotation = annotation
  )
  seurat.obj <- RenameCells(seurat.obj, add.cell.id = sample0)
  seurat.obj$DonorID         <- sp$Donor.ID
  seurat.obj$SortingStrategy <- sp$SortingStrategy

  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- TSSEnrichment(seurat.obj)

  seurat.obj <- subset(
    seurat.obj,
    subset = nCount_ATAC  < QC_MAX_ATAC &
             nCount_RNA   < QC_MAX_RNA  &
             nCount_ATAC  > QC_MIN_ATAC &
             nCount_RNA   > QC_MIN_RNA  &
             percent.mt   < QC_MAX_MT
  )

  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst",
                                      nfeatures = N_HVG)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)

  DefaultAssay(seurat.obj) <- "ATAC"
  seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = 5)
  seurat.obj <- RunTFIDF(seurat.obj)
  seurat.obj <- RunSVD(seurat.obj)

  seurat.obj <- FindMultiModalNeighbors(
    object               = seurat.obj,
    reduction.list       = list("pca", "lsi"),
    dims.list            = list(1:N_PCA_DIMS, N_LSI_DIMS),
    modality.weight.name = "RNA.weight",
    verbose              = TRUE
  )
  seurat.obj <- RunUMAP(seurat.obj, nn.name = "weighted.nn",
                         assay = "RNA", verbose = TRUE)

  DefaultAssay(seurat.obj) <- "RNA"
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:N_PCA_DIMS)
  seurat.obj <- FindClusters(seurat.obj,  resolution = CLUSTER_RES)

  p_umap <- DimPlot(seurat.obj, label = TRUE, repel = TRUE,
                     reduction = "umap") +
    NoLegend() + ggtitle(sample0)
  ggsave(p_umap, device = "eps", width = 9, height = 7,
         filename = file.path(dir.figures,
                               paste0(sample0, "_UMAP_doubletRemoved.eps")))

  saveRDS(seurat.obj,
          file = file.path(dir.seurat.samples,
                            paste0(sample0, "_JointSeuratObj_doubletRemoved.rds")))
  cat("Finished:", sample0, "\n")
}
cat("\nAll samples processed.\n")


# =============================================================================
# SECTION 3: Merge all per-sample objects into one combined object
# =============================================================================
# Reads back each per-sample RDS and assembles:
#   - rna.mtx   : column-merged RNA count matrix (sparse, genes padded with 0)
#   - atac.mtx  : column-merged ATAC count matrix (re-loaded from scATACpro
#                 reConstruct_matrix to ensure a common peak set)
#   - meta      : row-bound metadata from each sample
#   - scATAC.meta : per-barcode ATAC QC metrics (frac_peak, etc.)
#
# NOTE: rBind.fill pads missing features with 0 so matrices from samples with
# different peak sets can be safely cbind-ed.

data.path  <- readRDS(file.path(dir.metadata, "HD_DataPath.rds"))
annotation <- readRDS(path.annotation)
sample.name <- data.path$sample.name

cat("Merging sample:", sample.name[1], "\n")
seurat.tmp <- readRDS(file.path(dir.seurat.samples,
                                 paste0(sample.name[1], "_JointSeuratObj_doubletRemoved.rds")))
seurat.tmp$sample.name <- sample.name[1]
meta      <- seurat.tmp@meta.data

scATAC.meta <- fread(file.path(dir.scATACpro, sample.name[1],
                                "output/summary",
                                paste0(sample.name[1], ".MACS2.qc_per_barcode.txt")))
scATAC.meta$bc <- paste0(sample.name[1], "_", scATAC.meta$bc)
scATAC.meta    <- scATAC.meta[scATAC.meta$bc %in% rownames(meta), ]
scATAC.meta    <- data.frame(scATAC.meta); rownames(scATAC.meta) <- scATAC.meta$bc

rna.mtx  <- seurat.tmp[["RNA"]]@counts

atac.mtx <- readRDS(file.path(dir.scATACpro, sample.name[1],
                               "output/reConstruct_matrix_HD_7Donors/matrix.rds"))
colnames(atac.mtx) <- paste0(sample.name[1], "_", colnames(atac.mtx))
atac.mtx <- atac.mtx[, colnames(rna.mtx)]

for (sample0 in sample.name[-1]) {
  cat("Merging sample:", sample0, "\n")
  seurat.tmp <- readRDS(file.path(dir.seurat.samples,
                                   paste0(sample0, "_JointSeuratObj_doubletRemoved.rds")))
  seurat.tmp$sample.name <- sample0
  meta0 <- seurat.tmp@meta.data

  scATAC.meta0 <- fread(file.path(dir.scATACpro, sample0, "output/summary",
                                   paste0(sample0, ".MACS2.qc_per_barcode.txt")))
  scATAC.meta0$bc <- paste0(sample0, "_", scATAC.meta0$bc)
  scATAC.meta0    <- scATAC.meta0[scATAC.meta0$bc %in% rownames(meta0), ]
  scATAC.meta0    <- data.frame(scATAC.meta0); rownames(scATAC.meta0) <- scATAC.meta0$bc

  rna.mtx0  <- seurat.tmp[["RNA"]]@counts
  atac.mtx0 <- readRDS(file.path(dir.scATACpro, sample0,
                                  "output/reConstruct_matrix_HD_7Donors/matrix.rds"))
  colnames(atac.mtx0) <- paste0(sample0, "_", colnames(atac.mtx0))
  atac.mtx0 <- atac.mtx0[, colnames(rna.mtx0)]

  meta       <- rbind(meta, meta0)
  scATAC.meta <- rbind(scATAC.meta, scATAC.meta0)
  rna.mtx    <- t(rBind.fill(t(rna.mtx),  t(rna.mtx0),  fill = 0))
  atac.mtx   <- t(rBind.fill(t(atac.mtx), t(atac.mtx0), fill = 0))
}

# Drop columns recalculated below; drop irrelevant ATAC QC fields
meta$nCount_RNA <- meta$nCount_ATAC <- meta$nFeature_ATAC <-
  meta$nFeature_RNA <- meta$RNA.weight <- meta$ATAC.weight <-
  meta$RNA_snn_res.0.5 <- meta$seurat_clusters <- NULL

scATAC.meta$bc           <- NULL
scATAC.meta$tss_enrich_score <- NULL
scATAC.meta$frac_mito    <- NULL

# Build combined multimodal Seurat object
seurat.HDRef <- CreateSeuratObject(counts = rna.mtx, assay = "RNA")
seurat.HDRef[["ATAC"]] <- CreateChromatinAssay(
  counts     = atac.mtx,
  sep        = c(":", "-"),
  fragments  = path.fragments.combined,
  annotation = annotation
)
seurat.HDRef <- AddMetaData(seurat.HDRef, metadata = meta)
seurat.HDRef <- AddMetaData(seurat.HDRef, metadata = scATAC.meta)

# Label batch origin (PseudoNormal = enriched/sorted; TrueNormal = fresh)
seurat.HDRef$Batch <- seurat.HDRef$sample.name
seurat.HDRef$Batch[seurat.HDRef$Batch %in%
                     c("2986_Live", "2986_Enriched", "2932_Live", "2932_Enriched",
                       "3757_Live", "3757_Enriched", "3852_Enriched", "3852_Live")] <- "PseudoNormal"
seurat.HDRef$Batch[seurat.HDRef$Batch %in%
                     c("4970_Enriched", "4970_Live", "5123_Combine")] <- "TrueNormal"

VlnPlot(seurat.HDRef, features = c("nCount_RNA", "nCount_ATAC", "frac_peak", "percent.mt"),
        group.by = "sample.name", pt.size = 0)

rm(rna.mtx, atac.mtx, rna.mtx0, atac.mtx0, seurat.tmp)


# =============================================================================
# SECTION 4: RNA processing — Gini-iterative feature selection + PCA
# =============================================================================
# Strategy:
#   1. Standard HVG selection (2 000 features), remove lowly expressed genes.
#   2. Scale regressing out percent.mt and nCount_RNA; run PCA + initial UMAP.
#   3. Repeat (niter times): replace HVGs with Gini-informative genes per
#      current clusters, then re-scale, re-PCA, re-cluster.
#   Gini selection captures genes with highly unequal expression across
#   clusters, which outperforms variance-based selection for rare cell types.

DefaultAssay(seurat.HDRef) <- "RNA"
seurat.HDRef$percent.mt <- seurat.HDRef$percent.mt / 100  # convert to fraction if stored as %

seurat.HDRef <- NormalizeData(seurat.HDRef)
seurat.HDRef <- FindVariableFeatures(seurat.HDRef, nfeatures = N_HVG)

# Expression prevalence filter: remove genes expressed in < 1% of cells
cell.express.perc <- rowSums(seurat.HDRef@assays$RNA@counts >= 1) /
                       ncol(seurat.HDRef@assays$RNA@counts)
sel.features <- setdiff(VariableFeatures(seurat.HDRef),
                         names(cell.express.perc)[cell.express.perc < 0.01])
VariableFeatures(seurat.HDRef) <- sel.features

seurat.HDRef <- ScaleData(seurat.HDRef, features = VariableFeatures(seurat.HDRef),
                           vars.to.regress = c("percent.mt", "nCount_RNA"))
seurat.HDRef <- RunPCA(seurat.HDRef, npc = 50, verbose = FALSE,
                        features = VariableFeatures(seurat.HDRef))
seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:20)
seurat.HDRef <- FindNeighbors(seurat.HDRef, dims = 1:20, reduction = "pca")
seurat.HDRef <- FindClusters(seurat.HDRef, res = 0.1)

# Gini iteration
niter <- 2
for (k in seq_len(niter)) {
  gini_genes  <- ifg_select(seurat.HDRef@assays$RNA@counts,
                             seurat.HDRef$seurat_clusters,
                             gini_cut_qt = 0.9)$include_g
  sel.features <- setdiff(gini_genes,
                           names(cell.express.perc)[cell.express.perc < 0.01])
  VariableFeatures(seurat.HDRef) <- sel.features
  cat("Iteration", k, "— Gini feature count:", length(sel.features), "\n")

  seurat.HDRef <- ScaleData(seurat.HDRef, features = VariableFeatures(seurat.HDRef),
                             vars.to.regress = c("percent.mt", "nCount_RNA"))
  seurat.HDRef <- RunPCA(seurat.HDRef, npc = 50, verbose = FALSE,
                          features = VariableFeatures(seurat.HDRef))
  seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:20)
  seurat.HDRef <- FindNeighbors(seurat.HDRef, dims = 1:20, reduction = "pca")
  seurat.HDRef <- FindClusters(seurat.HDRef, res = 0.1)
}

# Generate multi-dimension UMAPs for comparison
for (d in c(10, 20, 30, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:d,
                           reduction.name = paste0("UMAP_pca", d),
                           reduction.key  = paste0("UMAP_pca", d, "_"))
}

# Diagnostic marker check on RNA UMAP
FeaturePlot(seurat.HDRef,
            features  = c("CD34", "CD38", "IGLL1", "DNTT", "MME", "CD19",
                          "MS4A1", "CD14", "CD3E", "GNLY", "GATA1", "MPO"),
            reduction = "UMAP_pca20") &
  scale_colour_gradientn(colours = brewer.pal(9, "Reds"))


# =============================================================================
# SECTION 5: ATAC processing — DAP-based feature selection + LSI
# =============================================================================
# Strategy:
#   1. Initial LSI on top features (q90 cutoff); cluster; find DAPs per cluster.
#   2. Filter low-quality peaks (max frequency across samples < 0.02, or
#      overall accessibility < 1%).
#   3. Re-run LSI on the DAP-selected features only; explore multiple LSI dims.
#
# A second DAP round (using the DAP-LSI-derived clusters) refines the feature
# set further, which improves separation of rare progenitor populations.

DefaultAssay(seurat.HDRef) <- "ATAC"

# Round 1: initial LSI for clustering
seurat.HDRef <- FindTopFeatures(seurat.HDRef, min.cutoff = "q90")
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
DepthCor(seurat.HDRef)

seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:20,
                         reduction.name = "UMAP_lsi20", reduction.key = "UMAP_lsi20_")
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:20)
seurat.HDRef <- FindClusters(seurat.HDRef, verbose = FALSE, algorithm = 1, resolution = 0.2)

# Cluster-level DAP discovery
all.daps <- FindAllMarkers(seurat.HDRef, test.use = "wilcox",
                            logfc.threshold = 0, slot = "data",
                            max.cells.per.ident = 200,
                            min.pct = 0.1, only.pos = TRUE,
                            min.diff.pct = 0.05)
all.daps <- data.table(all.daps, keep.rownames = TRUE)

# Peak quality filter: remove peaks absent in all samples or globally rare
mtx <- seurat.HDRef@assays$ATAC@counts
freq_by_batch <- sapply(unique(seurat.HDRef$sample.name), function(x)
  Matrix::rowMeans(mtx[, seurat.HDRef$sample.name == x] > 0))
freq.max     <- apply(freq_by_batch, 1, max)
names(freq.max) <- rownames(freq_by_batch)
filter.pks1  <- names(which(freq.max < 0.02))
rs           <- Matrix::rowSums(mtx > 0)
filter.pks2  <- names(which(rs <= 0.01 * ncol(mtx)))
filter.pks   <- unique(c(filter.pks1, filter.pks2))

sele.daps    <- all.daps[p_val_adj < 0.01 & avg_log2FC > 1]
sele.features <- unique(setdiff(sele.daps$gene, filter.pks))
cat("DAP round 1 — selected features:", length(sele.features), "\n")

# Round 1 DAP-based LSI
VariableFeatures(seurat.HDRef) <- sele.features
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
DepthCor(seurat.HDRef)

for (d in c(10, 20, 30, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:d,
                           reduction.name = paste0("UMAP_lsi", d),
                           reduction.key  = paste0("UMAP_lsi", d, "_"))
}
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:20)
seurat.HDRef <- FindClusters(seurat.HDRef, verbose = FALSE, algorithm = 1, resolution = 0.2)

# Round 2 DAP: refine with cluster from DAP-LSI space
all.daps.2 <- FindAllMarkers(seurat.HDRef, test.use = "wilcox",
                              logfc.threshold = 0, slot = "data",
                              max.cells.per.ident = 200,
                              min.pct = 0.05, only.pos = TRUE,
                              min.diff.pct = 0.05)
all.daps.2   <- data.table(all.daps.2, keep.rownames = TRUE)
sele.daps.2  <- all.daps.2[p_val_adj < 0.01 & avg_log2FC > 1]
sele.features.2 <- unique(setdiff(sele.daps.2$gene, filter.pks))
cat("DAP round 2 — selected features:", length(sele.features.2), "\n")

VariableFeatures(seurat.HDRef) <- sele.features.2
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
DepthCor(seurat.HDRef)

for (d in c(10, 15, 20, 25, 30, 35, 40, 45, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:d,
                           reduction.name = paste0("UMAP_lsi", d),
                           reduction.key  = paste0("UMAP_lsi", d, "_"))
}
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:20)
seurat.HDRef <- FindClusters(seurat.HDRef, verbose = FALSE, algorithm = 1, resolution = 0.5)


# =============================================================================
# SECTION 6: WNN joint embedding — multi-k neighbors + multi-resolution UMAPs
# =============================================================================
# Build four WNN graphs with different k values (10, 20, 30, 50) and generate
# corresponding UMAPs. Use dims PCA 1:20 and LSI 2:25 (DAP-refined LSI space).
# Explore RNA modality weight to confirm balanced contribution of both assays.

dims.pca.wnn <- 1:20
dims.lsi.wnn <- 2:25

for (k_nn in c(20, 30, 10, 50)) {
  nn.name <- paste0("weighted.nn.", k_nn)
  seurat.HDRef <- FindMultiModalNeighbors(
    object               = seurat.HDRef,
    reduction.list       = list("pca", "lsi"),
    dims.list            = list(dims.pca.wnn, dims.lsi.wnn),
    weighted.nn.name     = nn.name,
    k.nn                 = k_nn,
    modality.weight.name = "RNA.weight",
    verbose              = TRUE
  )
}

for (k_nn in c(20, 30, 10, 50)) {
  nn.name   <- paste0("weighted.nn.", k_nn)
  red.name  <- paste0("UMAP_wnn", k_nn)
  red.key   <- paste0("UMAP_wnn", k_nn, "_")
  seurat.HDRef <- RunUMAP(seurat.HDRef, nn.name = nn.name, assay = "RNA",
                           reduction.name = red.name, reduction.key = red.key,
                           verbose = TRUE)
}

# Compare k values and batch structure
p.wnn <- lapply(c(10, 20, 30, 50), function(k) {
  DimPlot(seurat.HDRef, label = TRUE, group.by = "sample.name", repel = TRUE,
          reduction = paste0("UMAP_wnn", k)) + NoLegend() + ggtitle(paste0("k=", k))
})
p.wnn[[1]] + p.wnn[[2]] + p.wnn[[3]] + p.wnn[[4]]

FeaturePlot(seurat.HDRef, features = c("RNA.weight", "ATAC.weight"),
            reduction = "UMAP_wnn20") &
  scale_colour_gradientn(colours = brewer.pal(9, "Reds"))

saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.combined,
                          "Seurat_AllDonors_JointAnalysis_RNAGini_ATACDAP.rds"))


# =============================================================================
# SECTION 7: Label transfer from old RNA reference
# =============================================================================
# Anchors between the old 5-donor scRNA reference and the new multiome object
# allow probabilistic cell-type transfer and pseudotime projection (MapQuery).
# The transferred labels (predicted.celltype) guide manual re-annotation in
# Section 10

DefaultAssay(seurat.HDRef) <- "RNA"
seurat.rna.old.ref <- readRDS(path.old.rna.ref)

HD.rna.anchors <- FindTransferAnchors(
  reference          = seurat.rna.old.ref,
  query              = seurat.HDRef,
  features           = VariableFeatures(seurat.rna.old.ref),
  dims               = 1:20,
  reference.reduction = "pca"
)
predictions <- TransferData(anchorset = HD.rna.anchors,
                             refdata   = seurat.rna.old.ref$Ctype,
                             dims      = 1:20)
predictions <- predictions[, c("predicted.id", "prediction.score.max")]
names(predictions) <- c("predicted.celltype", "labeltransfer.score")
seurat.HDRef <- AddMetaData(seurat.HDRef, metadata = predictions)

# Project query onto reference UMAP and transfer pseudotime
seurat.rna.old.ref <- RunUMAP(seurat.rna.old.ref, dims = 1:20, reduction = "pca",
                               return.model = TRUE,
                               reduction.name = "UMAP_ref", reduction.key = "UMAP_ref_")
seurat.HDRef <- MapQuery(
  anchorset         = HD.rna.anchors,
  reference         = seurat.rna.old.ref,
  query             = seurat.HDRef,
  refdata           = list(celltype = "Ctype", BTraj.pseudotime = "pseudotime"),
  reference.reduction = "pca",
  reduction.model   = "UMAP_ref"
)

myColors <- colorRampPalette(brewer.pal(12, "Paired"))(16)
DimPlot(seurat.HDRef, label = TRUE, group.by = "predicted.celltype",
        repel = TRUE, reduction = "UMAP_wnn30", cols = myColors)

saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.combined,
                          "Seurat_AllDonors_JointAnalysis_RNAGini_ATACDAP.rds"))


# =============================================================================
# SECTION 8: Cell-cycle & heat-shock score regression (optional re-run)
# =============================================================================
# Re-build the Seurat object from scratch using stored matrices to add
# cell-cycle and heat-shock module scores as regression covariates.
# This removes confounding transcriptional programs that mask cell identity.


seurat.HDRef <- readRDS(file.path(dir.seurat.combined,
                                   "Seurat_AllDonors_JointAnalysis_RNAGini_ATACDAP.rds"))

# Extract matrices and metadata, then rebuild object to avoid memory bloat
rna.mtx  <- seurat.HDRef@assays$RNA@counts
atac.mtx <- seurat.HDRef@assays$ATAC@counts
meta     <- seurat.HDRef@meta.data
rm(seurat.HDRef)

seurat.HDRef <- CreateSeuratObject(counts = rna.mtx, assay = "RNA")
seurat.HDRef[["ATAC"]] <- CreateChromatinAssay(
  counts = atac.mtx, sep = c(":", "-"),
  fragments = path.fragments.combined, annotation = annotation)
seurat.HDRef <- AddMetaData(seurat.HDRef, metadata = meta)
rm(rna.mtx, atac.mtx)

# Cell cycle scoring
cycle3   <- fread(path.cell.cycle.genes, header = FALSE)$V1
s.genes  <- cycle3[1:43]
g2m.genes <- cycle3[44:97]
seurat.HDRef <- CellCycleScoring(seurat.HDRef, s.features = s.genes,
                                  g2m.features = g2m.genes, set.ident = TRUE)

# Heat-shock module score
heat.shock.gene <- fread(path.heat.shock.genes)$`Approved symbol`
seurat.HDRef    <- AddModuleScore(seurat.HDRef, features = list(heat.shock.gene),
                                   name = "HeatShock.Score")

# RNA: Gini iteration with expanded regression
DefaultAssay(seurat.HDRef) <- "RNA"
seurat.HDRef <- NormalizeData(seurat.HDRef)
seurat.HDRef <- FindVariableFeatures(seurat.HDRef, nfeatures = N_HVG)
cell.express.perc <- rowSums(seurat.HDRef@assays$RNA@counts >= 1) /
                       ncol(seurat.HDRef@assays$RNA@counts)
sel.features <- setdiff(VariableFeatures(seurat.HDRef),
                         names(cell.express.perc)[cell.express.perc < 0.01])
VariableFeatures(seurat.HDRef) <- sel.features

seurat.HDRef <- ScaleData(seurat.HDRef, features = VariableFeatures(seurat.HDRef),
                           vars.to.regress = c("S.Score", "G2M.Score",
                                               "percent.mt", "HeatShock.Score1",
                                               "nCount_RNA"))
seurat.HDRef <- RunPCA(seurat.HDRef, npc = 50, verbose = FALSE,
                        features = VariableFeatures(seurat.HDRef))
seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:20)
seurat.HDRef <- FindNeighbors(seurat.HDRef, dims = 1:20, reduction = "pca")
seurat.HDRef <- FindClusters(seurat.HDRef, res = 0.1)

for (k in seq_len(niter)) {
  gini_genes  <- ifg_select(seurat.HDRef@assays$RNA@counts,
                             seurat.HDRef$seurat_clusters,
                             gini_cut_qt = 0.9)$include_g
  sel.features <- setdiff(gini_genes,
                           names(cell.express.perc)[cell.express.perc < 0.01])
  VariableFeatures(seurat.HDRef) <- sel.features
  cat("CC/HS iteration", k, "— features:", length(sel.features), "\n")
  seurat.HDRef <- ScaleData(seurat.HDRef, features = VariableFeatures(seurat.HDRef),
                             vars.to.regress = c("S.Score", "G2M.Score",
                                                 "percent.mt", "HeatShock.Score1",
                                                 "nCount_RNA"))
  seurat.HDRef <- RunPCA(seurat.HDRef, npc = 50, verbose = FALSE,
                          features = VariableFeatures(seurat.HDRef))
  seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:20)
  seurat.HDRef <- FindNeighbors(seurat.HDRef, dims = 1:20, reduction = "pca")
  seurat.HDRef <- FindClusters(seurat.HDRef, res = 0.1)
}

for (d in c(10, 20, 30, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:d,
                           reduction.name = paste0("UMAP_pca", d),
                           reduction.key  = paste0("UMAP_pca", d, "_"))
}

# ATAC: re-run DAP pipeline
DefaultAssay(seurat.HDRef) <- "ATAC"
seurat.HDRef <- FindTopFeatures(seurat.HDRef, min.cutoff = "q90")
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:20,
                         reduction.name = "UMAP_lsi20", reduction.key = "UMAP_lsi20_")
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:20)
seurat.HDRef <- FindClusters(seurat.HDRef, verbose = FALSE, algorithm = 1, resolution = 0.2)

all.daps.cc <- FindAllMarkers(seurat.HDRef, test.use = "wilcox",
                               logfc.threshold = 0, slot = "data",
                               max.cells.per.ident = 200,
                               min.pct = 0.1, only.pos = TRUE,
                               min.diff.pct = 0.05)
all.daps.cc <- data.table(all.daps.cc, keep.rownames = TRUE)
sele.daps.cc <- all.daps.cc[p_val_adj < 0.01 & avg_log2FC > 1]
sele.features.cc <- unique(setdiff(sele.daps.cc$gene, filter.pks))

VariableFeatures(seurat.HDRef) <- sele.features.cc
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
for (d in c(10, 20, 30, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:d,
                           reduction.name = paste0("UMAP_lsi", d),
                           reduction.key  = paste0("UMAP_lsi", d, "_"))
}

# WNN
for (k_nn in c(20, 30, 10, 50)) {
  seurat.HDRef <- FindMultiModalNeighbors(
    seurat.HDRef, reduction.list = list("pca", "lsi"),
    dims.list = list(dims.pca.wnn, dims.lsi.wnn),
    weighted.nn.name = paste0("weighted.nn.", k_nn),
    k.nn = k_nn, modality.weight.name = "RNA.weight", verbose = TRUE)
  seurat.HDRef <- RunUMAP(seurat.HDRef,
                           nn.name = paste0("weighted.nn.", k_nn), assay = "RNA",
                           reduction.name = paste0("UMAP_wnn", k_nn),
                           reduction.key  = paste0("UMAP_wnn", k_nn, "_"),
                           verbose = TRUE)
}

saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.combined,
                          "Seurat_AllDonors_JointAnalysis_RNAGini_regressCCandHS_ATACDAP.rds"))


# =============================================================================
# SECTION 9: B-trajectory downsampling and final joint analysis
# =============================================================================
# B-lineage cells dominate the dataset. Downsample B cells + CLP to 15 000 to
# balance cell-type representation before building the final reference.
# After downsampling, reprocess RNA (Gini) and ATAC (DAP) from scratch so that
# feature selection reflects the downsampled composition.

seurat.HDRef <- readRDS(file.path(dir.seurat.combined,
                                   "Seurat_AllDonors_JointAnalysis_RNAGini_ATACDAP.rds"))

# Bin cell types for proportional downsampling
seurat.HDRef$binned.celltype <- as.character(seurat.HDRef$predicted.celltype)
seurat.HDRef$binned.celltype[seurat.HDRef$binned.celltype %in%
                               c("Pre-pro-B", "Pro-B", "Pre-B", "Immature-B", "Mature-B")] <- "B"
seurat.HDRef$binned.celltype[seurat.HDRef$binned.celltype %in%
                               c("DC-Progenitor", "GMP", "MEP")] <- "Prog.Mye"
seurat.HDRef$binned.celltype[seurat.HDRef$binned.celltype %in% c("CLP")] <- "Prog.Lym"
seurat.HDRef$binned.celltype[seurat.HDRef$binned.celltype %in%
                               c("HSPC", "LMPP")] <- "Prog.Multi"
seurat.HDRef$binned.celltype[seurat.HDRef$binned.celltype %in%
                               c("cDC", "pDC", "Mono")] <- "Mye"
seurat.HDRef$binned.celltype[seurat.HDRef$binned.celltype %in%
                               c("T", "NK")] <- "T/NK"

other.cell   <- rownames(seurat.HDRef@meta.data[
  seurat.HDRef@meta.data$binned.celltype %notin% c("B", "Prog.Lym"), ])
B.downsampled <- sample(rownames(seurat.HDRef@meta.data[
  seurat.HDRef@meta.data$binned.celltype %in% c("B", "Prog.Lym"), ]), 15000)
selected.cell <- unique(c(other.cell, B.downsampled))

cat("Cells before downsampling:", ncol(seurat.HDRef),
    "\nCells after:", length(selected.cell), "\n")

# Rebuild from subsetted matrices
rna.mtx  <- seurat.HDRef@assays$RNA@counts[, selected.cell]
atac.mtx <- seurat.HDRef@assays$ATAC@counts[, selected.cell]
meta     <- seurat.HDRef@meta.data[selected.cell, ]
rm(seurat.HDRef)

seurat.HDRef <- CreateSeuratObject(counts = rna.mtx, assay = "RNA")
seurat.HDRef[["ATAC"]] <- CreateChromatinAssay(
  counts = atac.mtx, sep = c(":", "-"),
  fragments = path.fragments.combined, annotation = annotation)
seurat.HDRef <- AddMetaData(seurat.HDRef, metadata = meta)
rm(rna.mtx, atac.mtx)

# RNA Gini processing (same as Section 4, without CC/HS regression)
DefaultAssay(seurat.HDRef) <- "RNA"
seurat.HDRef <- NormalizeData(seurat.HDRef)
seurat.HDRef <- FindVariableFeatures(seurat.HDRef, nfeatures = N_HVG)
cell.express.perc <- rowSums(seurat.HDRef@assays$RNA@counts >= 1) /
                       ncol(seurat.HDRef@assays$RNA@counts)
sel.features <- setdiff(VariableFeatures(seurat.HDRef),
                         names(cell.express.perc)[cell.express.perc < 0.01])
VariableFeatures(seurat.HDRef) <- sel.features
seurat.HDRef <- ScaleData(seurat.HDRef, features = VariableFeatures(seurat.HDRef),
                           vars.to.regress = c("percent.mt", "nCount_RNA"))
seurat.HDRef <- RunPCA(seurat.HDRef, npc = 50, verbose = FALSE,
                        features = VariableFeatures(seurat.HDRef))
seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:20)
seurat.HDRef <- FindNeighbors(seurat.HDRef, dims = 1:20, reduction = "pca")
seurat.HDRef <- FindClusters(seurat.HDRef, res = 0.1)

for (k in seq_len(niter)) {
  gini_genes  <- ifg_select(seurat.HDRef@assays$RNA@counts,
                             seurat.HDRef$seurat_clusters,
                             gini_cut_qt = 0.9)$include_g
  sel.features <- setdiff(gini_genes,
                           names(cell.express.perc)[cell.express.perc < 0.01])
  VariableFeatures(seurat.HDRef) <- sel.features
  cat("Downsample Gini iteration", k, "— features:", length(sel.features), "\n")
  seurat.HDRef <- ScaleData(seurat.HDRef, features = VariableFeatures(seurat.HDRef),
                             vars.to.regress = c("percent.mt", "nCount_RNA"))
  seurat.HDRef <- RunPCA(seurat.HDRef, npc = 50, verbose = FALSE,
                          features = VariableFeatures(seurat.HDRef))
  seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:20)
  seurat.HDRef <- FindNeighbors(seurat.HDRef, dims = 1:20, reduction = "pca")
  seurat.HDRef <- FindClusters(seurat.HDRef, res = 0.1)
}

for (d in c(10, 20, 30, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, dims = 1:d,
                           reduction.name = paste0("UMAP_pca", d),
                           reduction.key  = paste0("UMAP_pca", d, "_"))
}

saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.downsample,
                          "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini.rds"))

# ATAC DAP processing on downsampled object (2-round DAP, same as Section 5)
seurat.HDRef <- readRDS(file.path(dir.seurat.downsample,
                                   "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini.rds"))
DefaultAssay(seurat.HDRef) <- "ATAC"
seurat.HDRef <- FindTopFeatures(seurat.HDRef, min.cutoff = "q90")
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:20,
                         reduction.name = "UMAP_lsi20", reduction.key = "UMAP_lsi20_")
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:20)
seurat.HDRef <- FindClusters(seurat.HDRef, verbose = FALSE, algorithm = 1, resolution = 0.2)

all.daps.ds <- FindAllMarkers(seurat.HDRef, test.use = "wilcox",
                               logfc.threshold = 0, slot = "data",
                               max.cells.per.ident = 200,
                               min.pct = 0.05, only.pos = TRUE,
                               min.diff.pct = 0.05)
all.daps.ds <- data.table(all.daps.ds, keep.rownames = TRUE)
saveRDS(all.daps.ds, file.path(dir.seurat.downsample, "DAP_round1_DownsampleBtraj15000.rds"))

mtx <- seurat.HDRef@assays$ATAC@counts
freq_by_batch <- sapply(unique(seurat.HDRef$sample.name), function(x)
  Matrix::rowMeans(mtx[, seurat.HDRef$sample.name == x] > 0))
freq.max   <- apply(freq_by_batch, 1, max); names(freq.max) <- rownames(freq_by_batch)
filter.pks1 <- names(which(freq.max < 0.02))
rs <- Matrix::rowSums(mtx > 0)
filter.pks2 <- names(which(rs <= 0.01 * ncol(mtx)))
filter.pks  <- unique(c(filter.pks1, filter.pks2))

sele.daps.ds <- all.daps.ds[p_val_adj < 0.01 & avg_log2FC > 1]
sele.features.ds <- unique(setdiff(sele.daps.ds$gene, filter.pks))

VariableFeatures(seurat.HDRef) <- sele.features.ds
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:20)
seurat.HDRef <- FindClusters(seurat.HDRef, verbose = FALSE, algorithm = 1, resolution = 0.2)

# DAP round 2
all.daps.ds.2 <- FindAllMarkers(seurat.HDRef, test.use = "wilcox",
                                 logfc.threshold = 0, slot = "data",
                                 max.cells.per.ident = 200,
                                 min.pct = 0.05, only.pos = TRUE,
                                 min.diff.pct = 0.05)
all.daps.ds.2 <- data.table(all.daps.ds.2, keep.rownames = TRUE)
saveRDS(all.daps.ds.2, file.path(dir.seurat.downsample, "DAP_round2_DownsampleBtraj15000.rds"))

sele.daps.ds.2  <- all.daps.ds.2[p_val_adj < 0.01 & avg_log2FC > 1]
sele.features.ds.2 <- unique(setdiff(sele.daps.ds.2$gene, filter.pks))

VariableFeatures(seurat.HDRef) <- sele.features.ds.2
seurat.HDRef <- RunTFIDF(seurat.HDRef)
seurat.HDRef <- RunSVD(seurat.HDRef)
for (d in c(10, 15, 20, 25, 30, 35, 40, 45, 50)) {
  seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:d,
                           reduction.name = paste0("UMAP_lsi", d),
                           reduction.key  = paste0("UMAP_lsi", d, "_"))
}

# WNN using PCA 1:20 and LSI 2:25
for (k_nn in c(20, 30, 10, 50)) {
  seurat.HDRef <- FindMultiModalNeighbors(
    seurat.HDRef, reduction.list = list("pca", "lsi"),
    dims.list = list(1:20, 2:25),
    weighted.nn.name = paste0("weighted.nn.", k_nn),
    k.nn = k_nn, modality.weight.name = "RNA.weight", verbose = TRUE)
  seurat.HDRef <- RunUMAP(seurat.HDRef,
                           nn.name = paste0("weighted.nn.", k_nn), assay = "RNA",
                           reduction.name = paste0("UMAP_wnn", k_nn),
                           reduction.key  = paste0("UMAP_wnn", k_nn, "_"),
                           verbose = TRUE)
}

saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.downsample,
                          "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP.rds"))


# =============================================================================
# SECTION 10: Manual cell-type annotation and quality filtering
# =============================================================================
# Combines:
#   (a) wsnn-based clustering at res=3
#   (b) label-transferred cell types from Section 7
#   (c) manual inspection via VlnPlot / FeaturePlot of canonical markers to assign AdjustLabel to each cluster.

seurat.HDRef <- readRDS(file.path(dir.seurat.downsample,
                                   "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP.rds"))

DefaultAssay(seurat.HDRef) <- "ATAC"
seurat.HDRef <- FindNeighbors(seurat.HDRef, reduction = "lsi", dims = 2:25)
seurat.HDRef <- FindClusters(seurat.HDRef, graph.name = "wsnn",
                              algorithm = 3, resolution = 3, verbose = FALSE)

# Inspect cluster vs. label-transfer concordance
table(seurat.HDRef$predicted.celltype, seurat.HDRef$seurat_clusters)

myColors <- colorRampPalette(brewer.pal(12, "Paired"))(17)
DimPlot(seurat.HDRef, group.by = "predicted.celltype",
        label = TRUE, reduction = "UMAP_wnn20", cols = myColors)

DefaultAssay(seurat.HDRef) <- "RNA"
VlnPlot(seurat.HDRef, features = c("MPO", "IGLL1", "ELANE", "LYZ",
                                    "MME",  "MS4A1", "DNTT", "CD19"))

# ── Manual cluster-to-celltype mapping ──────────────────────────────────────
# annotated

seurat.HDRef$AdjustLabel <- as.character(seurat.HDRef$seurat_clusters)
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("0","2","38","43")]              <- "HSC/MPP"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("4","35")]                       <- "LMPP"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("8","30")]                       <- "GMP"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("11","19","37","29","26","42")]  <- "MEP"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("15")]                           <- "CLP"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("3","4")]                        <- "Pro-B"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("14","20","32")]                 <- "Pre-B"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("1")]                            <- "Immature-B"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("22","36")]                      <- "Mature-B-1"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("23")]                           <- "Mature-B-2"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("6","27","25","21","33","39","41","12")] <- "T-cell"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("16","31")]                      <- "NK"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("5")]                            <- "DC-progenitor"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("7","28")]                       <- "pDC"
seurat.HDRef$AdjustLabel[seurat.HDRef$AdjustLabel %in% c("10","18","13","9","34","17","40")] <- "Monocyte"

seurat.HDRef$NewCellType <- as.character(seurat.HDRef$AdjustLabel)
seurat.HDRef$NewCellType[seurat.HDRef$seurat_clusters == "19" & seurat.HDRef$predicted.celltype == "HSC/MPP"] <- "HSC/MPP"
seurat.HDRef$NewCellType[seurat.HDRef$seurat_clusters == "15" & seurat.HDRef$predicted.celltype == "Pre-pro-B"] <- "Pre-pro-B"
seurat.HDRef$NewCellType[seurat.HDRef$seurat_clusters == "5" & seurat.HDRef$predicted.celltype == "cDC"] <- "cDC"

seurat.HDRef$CellType <- seurat.HDRef$NewCellType

# Rebuild UMAP models with return.model = TRUE
DefaultAssay(seurat.HDRef) <- "RNA"
seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "pca", dims = 1:20,
                         reduction.name = "UMAP_pca20", reduction.key = "UMAP_pca20_",
                         return.model = TRUE)
DefaultAssay(seurat.HDRef) <- "ATAC"
seurat.HDRef <- RunUMAP(seurat.HDRef, reduction = "lsi", dims = 2:25,
                         reduction.name = "UMAP_lsi25", reduction.key = "UMAP_lsi25_",
                         return.model = TRUE)
seurat.HDRef <- RunUMAP(seurat.HDRef, nn.name = "weighted.nn.20", assay = "RNA",
                         reduction.name = "UMAP_wnn20", reduction.key = "UMAP_wnn20_",
                         return.model = TRUE)

saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.downsample,
                          "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels.rds"))



# =============================================================================
# SECTION 11: Pseudotime trajectory inference with Slingshot
# =============================================================================
# Fits a principal curve through the UMAP_pca20 embedding, starting from
# HSC/MPP and ending at Mature-B, Monocyte, and MEP.
# Using the RNA-only PCA UMAP (rather than WNN) is intentional: it produces
# a smoother B-lineage curve without ATAC discontinuities.
# Only Lineage1 (the B-differentiation trajectory) is retained as Btraj.PseudoTime.

seurat.HDRef <- readRDS(file.path(dir.seurat.downsample,
                                   "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels.rds"))

# Merge Mature-B subtypes for trajectory endpoint
seurat.HDRef$CelltypeforPseudo <- as.character(seurat.HDRef$NewCellTypeUpdate)
seurat.HDRef$CelltypeforPseudo[seurat.HDRef$CelltypeforPseudo %in%
                                  c("Mature-B-1", "Mature-B-2")] <- "Mature-B"

HD_Trajectory <- slingshot(
  data          = Embeddings(seurat.HDRef, "UMAP_pca20"),
  clusterLabels = seurat.HDRef$CelltypeforPseudo,
  start.clus    = "HSC/MPP",
  end.clus      = c("Mature-B", "Monocyte", "MEP"),
  approx_points = 150,
  omega         = TRUE,
  omega_scale   = 2
)

HD_Trajectory@metadata$lineages

pseudo.time <- data.frame(HD_Trajectory@assays@data@listData$pseudotime)
pseudo.time <- pseudo.time["Lineage1"]          # B-differentiation trajectory
names(pseudo.time) <- "Btraj.PseudoTime"
seurat.HDRef <- AddMetaData(seurat.HDRef, metadata = pseudo.time)

FeaturePlot(seurat.HDRef, features = "Btraj.PseudoTime", reduction = "UMAP_wnn20") &
  scale_colour_gradientn(colours = brewer.pal(9, "Reds"))

saveRDS(HD_Trajectory,
        file = file.path(dir.trajectory,
                          "Slingshot_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels.rds"))
saveRDS(seurat.HDRef,
        file = file.path(dir.seurat.downsample,
                          "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels.rds"))


cat("\nHealthy Donor Multiome Reference pipeline complete.\n")