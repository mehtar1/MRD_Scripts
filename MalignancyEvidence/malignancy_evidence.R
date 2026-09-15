# =============================================================================
# Malignancy Evidence Pipeline
# =============================================================================
#
# Description:
#   Identifies and validates leukemic blast populations within a subtype-specific
#   scRNA-seq Seurat object. The same pipeline is applied to each B-ALL subtype
#   (KMT2Ar, iAMP21, Ph-like, Ph-positive) by changing the parameters in the
#   USER CONFIGURATION block below.
#
#   The pipeline covers five sequential analyses:
#     1. Data processing  — normalization, HVG selection, PCA, UMAP, clustering
#     2. Cell annotation  — Level 1 labels (Blast / T/NK / B / Erythroid / Mono)
#                           based on cluster IDs identified during visual inspection
#     3. Shannon entropy  — chromatin/transcriptional heterogeneity score per cell
#                           (requires entropycalc.R utility)
#     4. Healthy similarity — CCA-based label transfer from a healthy donor
#                              reference to score how "normal-like" each blast is
#     5. Blast percentage validation — correlates annotated blast % with
#                                       clinical blast count; Pearson r + t-test
#     6. InferCNV integration — loads pre-run inferCNV output and adds CNV
#                                scores per cell back into the Seurat object
#     7. Per-patient re-clustering — isolates ambiguous patients for fine-grained
#                                     cluster inspection and re-labeling
#
# Inputs (all set in USER CONFIGURATION):
#   - Subtype Seurat RDS (pre-merged, pre-QC'd RNA object)
#   - Healthy donor reference Seurat RDS
#   - Clinical blast percentage CSV (columns: USI, Stage, leukemia..)
#   - InferCNV output directory
#   - Per-patient barcode CSVs for cluster re-labeling (subtype-specific)
#   - Output directories for plots and R objects
#
# Outputs:
#   - Updated Seurat RDS with Level1.anno, entropy scores, and similarity scores
#   - PDF figures: UMAPs, FeaturePlots, VlnPlots, similarity boxplots,
#                  blast % scatter, entropy feature plots
#
# Usage:
#   Edit the USER CONFIGURATION block for your subtype, then run the script
#   end-to-end, OR source individual sections interactively.
#
#   Example (command line):
#     Rscript malignancy_evidence.R   # after editing config below
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(infercnv)

source("entropycalc.R")   # provides get.entropy(); included in this folder

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────
# Edit this block to switch between subtypes. All subtype-specific parameters
# are collected here; nothing else in the script needs to change.
#
# Supported subtypes (set SUBTYPE to one of these strings):
#   "KMT2Ar" | "iAMP21" | "Ph-like" | "Ph-positive"

SUBTYPE <- "KMT2Ar"   # <-- change this to run a different subtype

# Input files
path.seurat.input   <- "/path/to/RObjects/<subtype>.rna.1.rds"
path.healthy.ref    <- "/path/to/healthy_reference_seurat.rds"
path.blast.pct.csv  <- "/path/to/LeukemicBlastPercentage.csv"   # cols: USI, Stage, leukemia..
path.infercnv.dir   <- "/path/to/InferCNV_Results/<subtype>/"   # directory with inferCNV output
path.infercnv.rds   <- "/path/to/InferCNV_Results/<subtype>/<subtype>_downsampledSeuratforInferCNV.rds"

# Output directories (create before running)
dir.plots   <- "/path/to/Plots/<subtype>/"
dir.robjs   <- "/path/to/RObjects/"

# Dimensionality reduction parameters
# These differ slightly per subtype; set to match your original analysis.
N_PCA_DIMS   <- 75     # dims for UMAP and neighbor graph  (KMT2Ar=75, others=50 or 30)
N_PCA_MAX    <- 100    # total PCs computed
UMAP_NAME    <- paste0("umapPCA", N_PCA_DIMS)   # e.g. "umapPCA75"
CLUSTER_RES  <- 0.50   # initial clustering resolution

# Level 1 annotation cluster IDs
# Replace these vectors with the cluster numbers from YOUR visual inspection.
# All clusters not listed below will be labeled "Blast" by default.
# These are the subtype-specific values from the original scripts:
#
#   Subtype       T/NK clusters       B clusters   Erythroid clusters   Mono clusters
#   KMT2Ar        15,12,22,43,38,     28,19        16,34,49             31,32,51
#                 29,44,48,67,27
#   iAMP21        1,3,42              10           18,22,30             24
#   Ph-like       0,13                15           17,21                22
#   Ph-positive   13,7                17           14                   16

CLUSTERS_TNK        <- c(15, 12, 22, 43, 38, 29, 44, 48, 67, 27)
CLUSTERS_B          <- c(28, 19)
CLUSTERS_ERYTHROID  <- c(16, 34, 49)
CLUSTERS_MONO       <- c(31, 32, 51)

# Annotation resolution: which clustering column to use for Level1 labels
# (set to whichever resolution you inspected visually, e.g. "RNA_snn_res.2")
ANNO_RES_COL <- "RNA_snn_res.2"

# Healthy reference subsampling size (reduce if object is smaller than 50k cells)
HEALTHY_SUBSAMPLE_N <- 50000

# Per-patient cluster re-labeling
# Each entry is a named list with:
#   csv_path  : path to barcode CSV (rownames = cell barcodes)
#   patientID : patient ID string to filter on
#   anno      : original Level1.anno label to filter on (the wrongly labeled cells)
#   new_label : what to relabel them as ("Blast" or a normal cell type)
#
# Add one list entry per patient/cluster requiring correction.
# Leave as an empty list() if no re-labeling is needed.
recluster.corrections <- list(
  # Example (KMT2Ar has no per-patient re-labeling; iAMP21 example shown):
  # list(csv_path  = "/path/to/iamp21_Bcluster_1.csv",
  #      patientID = "PAVALZ",
  #      anno      = "B",
  #      new_label = "Blast"),
  # list(csv_path  = "/path/to/iamp21_Tcluster_1.csv",
  #      patientID = "PATSBP",
  #      anno      = "T/NK",
  #      new_label = "Blast")
)

# Patients that require individual sub-clustering for re-labeling inspection.
# Each entry: list(patientID, plot_title)
# Leave as empty list() if not needed.
patients.to.subcluster <- list(
  # list(patientID = "PAVALZ", plot_title = "PAVALZ_Bclust"),
  # list(patientID = "PATSBP", plot_title = "PATSBP_Tclust")
)


# ── Derived paths (auto-built from SUBTYPE; no edits needed below) ───────────

path.seurat.output         <- file.path(dir.robjs, paste0(tolower(SUBTYPE), ".rna.annotated.rds"))
path.anchors.rds           <- file.path(dir.robjs, paste0(tolower(SUBTYPE), ".healthy.donor.anchors.rds"))
path.similarity.rds        <- file.path(dir.robjs, paste0(tolower(SUBTYPE), "_downsample", HEALTHY_SUBSAMPLE_N, "K_withHealthySimilarity.rds"))

pdf.umap                   <- file.path(dir.plots, paste0(SUBTYPE, "_RNA_UMAPs.pdf"))
pdf.allmarkers             <- file.path(dir.plots, paste0(SUBTYPE, "_RNA_AllMarkers.pdf"))
pdf.entropy                <- file.path(dir.plots, paste0(SUBTYPE, "_ShannonEntropy.pdf"))
pdf.similarity             <- file.path(dir.plots, paste0(SUBTYPE, "_SimilarityScore.pdf"))
pdf.blast.pct              <- file.path(dir.plots, paste0(SUBTYPE, "_BlastPercentageValidation.pdf"))
pdf.subcluster             <- file.path(dir.plots, paste0(SUBTYPE, "_PerPatientSubclustering.pdf"))
pdf.umap.relabeled         <- file.path(dir.plots, paste0(SUBTYPE, "_RNA_UMAPs_relabeled.pdf"))

# Standard marker genes used for all subtypes
MARKER_GENES <- c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ")
myColors     <- colorRampPalette(brewer.pal(6, "Set1"))(6)


# =============================================================================
# SECTION 1: Data Processing
# =============================================================================
# Loads the subtype Seurat object and runs normalization, variable feature
# selection, PCA, UMAP, and Louvain clustering.

cat("Loading Seurat object for subtype:", SUBTYPE, "\n")
seurat.rna <- readRDS(path.seurat.input)

seurat.rna <- NormalizeData(seurat.rna)
seurat.rna <- FindVariableFeatures(seurat.rna, nfeatures = 4000)
seurat.rna <- RunPCA(seurat.rna, npcs = N_PCA_MAX, verbose = TRUE,
                     features = VariableFeatures(seurat.rna))
seurat.rna <- RunUMAP(seurat.rna, dims = 1:N_PCA_DIMS,
                      reduction.name = UMAP_NAME)
seurat.rna <- FindNeighbors(seurat.rna, dims = 1:N_PCA_DIMS)
seurat.rna <- FindClusters(seurat.rna, resolution = CLUSTER_RES)

# Quick diagnostic: marker gene expression on the initial clustering
FeaturePlot(seurat.rna, features = MARKER_GENES,
            reduction = UMAP_NAME, ncol = 3,
            label = TRUE, repel = TRUE, order = TRUE)

DotPlot(seurat.rna, features = MARKER_GENES)


# =============================================================================
# SECTION 2: Level 1 Cell Type Annotation
# =============================================================================
# Assigns broad cell type labels based on cluster IDs identified during visual
# inspection. Clusters not assigned to a normal population are labeled "Blast".
#
# NOTE: Review ANNO_RES_COL and the CLUSTERS_* vectors in USER CONFIGURATION
# before running this section. These are the only values that differ between
# subtypes.

cat("Annotating cell types...\n")
seurat.rna@meta.data$Level1.anno <- as.character(
  seurat.rna@meta.data[[ANNO_RES_COL]]
)

# All unassigned clusters default to Blast
all.normal.clusters <- c(CLUSTERS_TNK, CLUSTERS_B,
                          CLUSTERS_ERYTHROID, CLUSTERS_MONO)
seurat.rna@meta.data$Level1.anno[
  seurat.rna@meta.data$Level1.anno %notin% as.character(all.normal.clusters)
] <- "Blast"

seurat.rna@meta.data$Level1.anno[
  seurat.rna@meta.data$Level1.anno %in% as.character(CLUSTERS_TNK)
] <- "T/NK"

seurat.rna@meta.data$Level1.anno[
  seurat.rna@meta.data$Level1.anno %in% as.character(CLUSTERS_B)
] <- "B"

seurat.rna@meta.data$Level1.anno[
  seurat.rna@meta.data$Level1.anno %in% as.character(CLUSTERS_ERYTHROID)
] <- "Erythroid Progenitor"

seurat.rna@meta.data$Level1.anno[
  seurat.rna@meta.data$Level1.anno %in% as.character(CLUSTERS_MONO)
] <- "Mono"

cat("Level1.anno distribution:\n")
print(table(seurat.rna@meta.data$Level1.anno))

# -- Annotation plots --------------------------------------------------------
pdf(pdf.umap, height = 10, width = 10)
DimPlot(seurat.rna, group.by = "Level1.anno",    label = FALSE, cols = myColors,    raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
DimPlot(seurat.rna, group.by = "Stage",          cols = myColors[1:3],              raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
DimPlot(seurat.rna, group.by = "Subtype",        cols = myColors[1:4],              raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
DimPlot(seurat.rna, group.by = "MRDGroup",       cols = myColors[1:3],              raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
DimPlot(seurat.rna, group.by = "patientID",      label = TRUE,  repel = TRUE,       raster.dpi = c(1024, 1024),                reduction = UMAP_NAME) + NoLegend()
for (gene in MARKER_GENES) {
  print(FeaturePlot(seurat.rna, features = gene, order = TRUE, label = TRUE,
                    raster = TRUE, raster.dpi = c(1024, 1024),
                    label.size = 3, reduction = UMAP_NAME))
}
VlnPlot(seurat.rna, features = MARKER_GENES, stack = TRUE,
        group.by = "Level1.anno", pt.size = 0, flip = TRUE)
dev.off()

pdf(pdf.allmarkers, height = 11.35, width = 20)
FeaturePlot(seurat.rna, features = MARKER_GENES, order = TRUE, label = TRUE,
            raster = TRUE, raster.dpi = c(1024, 1024), label.size = 3, ncol = 3)
dev.off()


# =============================================================================
# SECTION 3: Shannon Entropy Score
# =============================================================================
# Computes per-cell transcriptional entropy across multiple clustering
# resolutions using get.entropy() from entropycalc.R. Higher entropy indicates
# greater transcriptional heterogeneity (associated with blast plasticity).

cat("Computing Shannon entropy scores...\n")

entropy.resolutions <- c(0.25, 0.5, 0.75, 1, 2)
for (res in entropy.resolutions) {
  seurat.rna <- FindClusters(seurat.rna, resolution = res)
}

entropy.cols <- paste0("RNA_snn_res.", entropy.resolutions)
seurat.rna   <- get.entropy(seurat.rna, entropy.cols, join.to.seurat = TRUE)

pdf(pdf.entropy, width = 10, height = 10)
FeaturePlot(seurat.rna,
            features = c(paste0(entropy.cols, "_entropy"), "average_entropy"),
            raster = TRUE, raster.dpi = c(1024, 1024),
            combine = FALSE, reduction = UMAP_NAME)
dev.off()


# =============================================================================
# SECTION 4: Similarity to Healthy Reference
# =============================================================================
# Subsamples the leukemia object, then uses Seurat CCA anchor transfer to
# project leukemic cells onto a healthy donor reference. The max prediction
# score per cell ("prediction.score.max") quantifies how normal-like each cell
# is; blast cells typically score lower than contaminating normal cells.

cat("Running healthy reference label transfer...\n")

# Load and prepare healthy reference
seurat.hd <- readRDS(path.healthy.ref)
seurat.hd <- NormalizeData(seurat.hd)
seurat.hd <- FindVariableFeatures(seurat.hd)
seurat.hd <- ScaleData(seurat.hd)

# Subsample leukemia object (CCA is memory-intensive on large objects)
n.cells <- ncol(seurat.rna)
subsample.size <- min(HEALTHY_SUBSAMPLE_N, n.cells)
if (subsample.size < n.cells) {
  cat("Subsampling to", subsample.size, "cells for anchor finding.\n")
}
seurat.sub <- seurat.rna[, sample(colnames(seurat.rna), subsample.size, replace = FALSE)]
seurat.sub <- NormalizeData(seurat.sub)
seurat.sub <- FindVariableFeatures(seurat.sub)
seurat.sub <- ScaleData(seurat.sub)

# Find transfer anchors and transfer cell type labels
anchors <- FindTransferAnchors(reference = seurat.hd, query = seurat.sub,
                                reduction = "cca")
saveRDS(anchors, file = path.anchors.rds)

predictions <- TransferData(anchorset = anchors,
                             refdata   = seurat.hd$Ctype,
                             k.weight  = 50,
                             weight.reduction = "cca")
seurat.sub <- AddMetaData(seurat.sub, metadata = predictions)
saveRDS(seurat.sub, file = path.similarity.rds)

# -- Similarity plots --------------------------------------------------------
pdf(pdf.similarity, width = 7, height = 7)
seurat.sub@meta.data %>%
  ggplot(aes(x = Level1.anno, y = prediction.score.max, fill = Level1.anno)) +
  geom_boxplot(outlier.size = 1) +
  labs(title  = paste(SUBTYPE, "scRNA-Seq Similarity to Healthy Reference"),
       x      = NULL,
       y      = "Similarity Score (max prediction score)") +
  theme(legend.position = "none")
dev.off()

VlnPlot(seurat.sub, group.by = "Level1.anno", features = "prediction.score.max")

FeaturePlot(seurat.sub, features = "prediction.score.max", order = TRUE) +
  DimPlot(seurat.sub, group.by = "predicted.id") +
  DimPlot(seurat.sub, group.by = "Stage") +
  DimPlot(seurat.sub, group.by = "patientID")


# =============================================================================
# SECTION 5: Blast Percentage Validation
# =============================================================================
# Computes the fraction of cells annotated as "Blast" per patient/stage from
# the Seurat object, then compares to clinical flow cytometry blast counts.
# Reports Pearson correlation and a two-sample t-test.

cat("Validating blast percentages against clinical counts...\n")

# Normal cell type labels (everything not listed is counted as Blast)
normal.labels <- c("Mono", "Erythroid Progenitor", "T/NK", "B")

blasts.meta <- seurat.rna@meta.data %>%
  mutate(NormalStatus = case_when(
    Level1.anno %in% normal.labels ~ "Normal",
    Level1.anno == "Blast"         ~ "Blast",
    TRUE                           ~ NA_character_
  )) %>%
  group_by(patientID, NormalStatus, Stage) %>%
  summarize(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = NormalStatus, values_from = count) %>%
  mutate(BlastLabeledPercent = Blast / (Blast + Normal))

# Join clinical blast percentages
blast.pct.clinical <- read.csv(path.blast.pct.csv)
blasts.meta <- blasts.meta %>%
  left_join(blast.pct.clinical, by = c("patientID" = "USI", "Stage")) %>%
  rename(countedLeukemia = "leukemia..") %>%
  mutate(diff_in_percents = countedLeukemia - BlastLabeledPercent) %>%
  select(-any_of("Subtype")) %>%
  unite(sampleID, c(patientID, Stage), remove = TRUE)

cor_r   <- round(cor(blasts.meta$BlastLabeledPercent,
                     blasts.meta$countedLeukemia,
                     method = "pearson"), digits = 3)
cor_lab <- paste0("Pearson r = ", cor_r)

pdf(pdf.blast.pct, width = 5, height = 5)
ggplot(blasts.meta, aes(x = BlastLabeledPercent, y = countedLeukemia)) +
  geom_point() +
  geom_smooth(method = lm, se = FALSE) +
  geom_text_repel(label = blasts.meta$sampleID) +
  labs(title   = SUBTYPE,
       x       = "Annotated Blast Percentage (scRNA)",
       y       = "Clinical Blast Percentage (flow)",
       caption = cor_lab)
dev.off()

t.test(x = blasts.meta$BlastLabeledPercent,
       y = blasts.meta$countedLeukemia,
       var.equal = TRUE)


# =============================================================================
# SECTION 6: InferCNV Integration
# =============================================================================
# Loads the pre-computed inferCNV output for this subtype and adds per-cell
# CNV scores (proportion_scaled_cnv_chrN, proportion_scaled_dupli_chrN) back
# into a downsampled Seurat object for visualization.
#
# NOTE: inferCNV must be run separately (see infercnv_rna_bySample.R).
# The chromosomes of interest vary by subtype:
#   KMT2Ar      → chr6 (duplication)
#   iAMP21      → chr21 (amplification), chr6
#   Ph-like     → chr13, chr19, chr22
#   Ph-positive → chr1, chr3

cat("Integrating inferCNV results...\n")

seurat.infercnv <- readRDS(path.infercnv.rds)
seurat.infercnv <- NormalizeData(seurat.infercnv)
seurat.infercnv <- add_to_seurat(seurat_obj        = seurat.infercnv,
                                  infercnv_output_path = path.infercnv.dir)

# Visualize CNV scores — replace chromosome names with those relevant to your subtype
DimPlot(seurat.infercnv, group.by = "Level1.anno")
DimPlot(seurat.infercnv, group.by = "top_dupli_4")

# Example chromosome of interest — update per subtype:
# FeaturePlot(seurat.infercnv, features = "proportion_scaled_dupli_chr6")
# FeaturePlot(seurat.infercnv, features = "proportion_scaled_cnv_chr21")


# =============================================================================
# SECTION 7: Per-patient Sub-clustering and Cluster Re-labeling
# =============================================================================
# Some patients have small clusters of cells initially mis-labeled as normal
# (e.g., a T/NK cluster that is actually malignant based on CNV or marker
# evidence). This section:
#   (a) Sub-clusters individual patients for visual review.
#   (b) Applies re-labeling decisions recorded in recluster.corrections.
#
# Both steps are driven entirely by USER CONFIGURATION — no code changes needed.

# -- (a) Per-patient sub-clustering for visual inspection --------------------
if (length(patients.to.subcluster) > 0) {
  cat("Sub-clustering individual patients...\n")
  pdf(pdf.subcluster, height = 11.35, width = 20)

  for (pt in patients.to.subcluster) {
    cells.pt  <- seurat.rna@meta.data %>%
      dplyr::filter(patientID == pt$patientID) %>% rownames()
    rna.pt    <- subset(seurat.rna, cells = cells.pt)
    rna.pt    <- NormalizeData(rna.pt)
    rna.pt    <- FindVariableFeatures(rna.pt, nfeatures = 4000)
    rna.pt    <- RunPCA(rna.pt, npcs = 30, verbose = FALSE,
                        features = VariableFeatures(rna.pt))
    rna.pt    <- RunUMAP(rna.pt, dims = 1:30)
    rna.pt    <- FindNeighbors(rna.pt, dims = 1:30)
    rna.pt    <- FindClusters(rna.pt, resolution = 0.25)

    print(
      DimPlot(rna.pt, group.by = "RNA_snn_res.0.25") + ggtitle(pt$plot_title) +
      DimPlot(rna.pt, group.by = "Level1.anno") +
      FeaturePlot(rna.pt, features = "MME",  order = TRUE) +
      FeaturePlot(rna.pt, features = "CD3D", order = TRUE) +
      FeaturePlot(rna.pt, features = "MS4A1", order = TRUE)
    )
  }
  dev.off()
}

# -- (b) Apply re-labeling corrections ---------------------------------------
# Loads per-patient barcode CSVs and reassigns Level1.anno for any cells
# whose marker expression or CNV evidence contradicts the initial label.

if (length(recluster.corrections) > 0) {
  cat("Applying cluster re-labeling corrections...\n")

  for (corr in recluster.corrections) {
    barcodes <- read.csv(corr$csv_path, row.names = 1) %>%
      dplyr::filter(Level1.anno == corr$anno,
                    patientID   == corr$patientID) %>%
      rownames()

    seurat.rna@meta.data$Level1.anno[
      rownames(seurat.rna@meta.data) %in% barcodes
    ] <- corr$new_label

    cat("  Re-labeled", length(barcodes), "cells from",
        corr$patientID, corr$anno, "→", corr$new_label, "\n")
  }

  cat("Level1.anno after re-labeling:\n")
  print(table(seurat.rna@meta.data$Level1.anno))

  # -- Final annotated UMAPs after re-labeling --------------------------------
  pdf(pdf.umap.relabeled, height = 10, width = 10)
  DimPlot(seurat.rna, group.by = "Level1.anno",    label = FALSE, cols = myColors,    raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
  DimPlot(seurat.rna, group.by = "Stage",          cols = myColors[1:3],              raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
  DimPlot(seurat.rna, group.by = "Subtype",        cols = myColors[1:4],              raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
  DimPlot(seurat.rna, group.by = "MRDGroup",       cols = myColors[1:3],              raster = TRUE, raster.dpi = c(1024, 1024), reduction = UMAP_NAME)
  DimPlot(seurat.rna, group.by = "patientID",      label = TRUE, repel = TRUE,        raster.dpi = c(1024, 1024),                reduction = UMAP_NAME) + NoLegend()
  for (gene in MARKER_GENES) {
    print(FeaturePlot(seurat.rna, features = gene, order = TRUE, label = TRUE,
                      raster = TRUE, raster.dpi = c(1024, 1024),
                      label.size = 3, reduction = UMAP_NAME))
  }
  VlnPlot(seurat.rna, features = MARKER_GENES, stack = TRUE,
          group.by = "Level1.anno", pt.size = 0, flip = TRUE)
  dev.off()
}


# =============================================================================
# SECTION 8: Save Final Annotated Object
# =============================================================================

saveRDS(seurat.rna, file = path.seurat.output)
cat("Saved annotated Seurat object to:", path.seurat.output, "\n")
