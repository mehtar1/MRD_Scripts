# =============================================================================
# Project Patient Cells onto Healthy Multiome Reference
# =============================================================================
#
# Description:
#   Projects leukemic patient scRNA-seq cells onto a pre-built healthy donor
#   UMAP reference using a direct PCA projection approach. For each sample:
#     1. The healthy donor PCA scaling coefficients (learned from regression
#        of mitochondrial % and UMI count) are recovered and applied to
#        patient data so both datasets occupy the same feature space.
#     2. Patient cells are projected into the healthy donor PCA space using
#        the shared gene loadings, then embedded into the healthy UMAP using
#        uwot::umap_transform().
#     3. Each patient cell is assigned the cell type of its nearest healthy
#        donor neighbor in UMAP space (Euclidean distance).
#     4. A per-sample projection plot and a cell-to-healthy-neighbor mapping
#        table are saved.
#
#   This script generalizes the original Ph-positive-only analysis to work
#   on any annotated subtype Seurat object by reading sample identities
#   directly from the object's metadata.
#
# Inputs (set in USER CONFIGURATION):
#   - Healthy donor Seurat RDS (joint RNA+ATAC, with Ctype annotation)
#   - Annotated patient Seurat RDS (output of malignancy_evidence.R;
#     must have $sample.name and $Level1.anno in metadata)
#   - Per-sample raw RNA count matrices (saved during preprocessing as
#     Matrix_forProjection/Matrix_related2_<sample>_JointSeuratObj_doubletRemoved.rds)
#
# Outputs:
#   - Per-sample projection plots (.eps) overlaying patient cells on healthy UMAP
#   - Per-sample cell-to-neighbor mapping tables (.csv)
#   - Combined mapping table across all samples (.csv)
#
# Dependencies:
#   Seurat, uwot, ggplot2, RColorBrewer, data.table
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(uwot)
library(ggplot2)
library(RColorBrewer)
library(data.table)

source("scDataAnalysis_Utilities.R")


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Healthy donor reference Seurat object (must have $predicted.celltype,
# $percent.mt, $nCount_RNA, and a PCA reduction)
path.healthy.ref    <- "/path/to/SeuratObj/HDMultiome/Combined/Seurat_AllDonors_JointAnalysis.rds"

# Annotated patient Seurat object (any subtype; must have $sample.name,
# $Level1.anno in metadata). To project a specific subtype only, filter
# the samples vector in the loop below.
path.patient.rds    <- "/path/to/RObjects/<subtype>.rna.annotated.rds"

# Directory containing per-sample raw RNA count matrix RDS files, named:
#   Matrix_related2_<sample.name>_JointSeuratObj_doubletRemoved.rds
dir.matrices        <- "/path/to/Matrix_forProjection/"

# Output directories
dir.plots           <- "/path/to/Figures/scRNA/Projection_Multiome/"
dir.metadata.sample <- "/path/to/metadata/Projection_MultiomeHD/BySample/"
path.output.combined <- "/path/to/metadata/Projection_MultiomeHD/AllSamples_proj2healthyCells.csv"

# UMAP model parameters — must match those used to build the healthy reference
UMAP_N_NEIGHBORS <- 30
UMAP_MIN_DIST    <- 0.3
UMAP_METRIC      <- "cosine"
N_PCA_DIMS       <- 20      # number of PCA dimensions used in projection

# Set to TRUE to project all cells; FALSE to project Blast cells only
PROJECT_ALL_CELLS <- TRUE


# =============================================================================
# SECTION 1: Recover healthy donor scaling model
# =============================================================================
# Extracts the PCA feature loadings and re-fits per-gene linear regression
# models (expr ~ percent.mt + nCount_RNA) on the healthy donor data. The
# resulting coefficients, gene means, and SDs are used in Section 3 to
# apply identical scaling to patient data before projection.
#
# NOTE: This recovers the transformation post-hoc from the stored Seurat
# object rather than requiring the original ScaleData output to be saved
# separately. This ensures reproducibility without extra storage.

cat("Loading healthy donor reference...\n")
seurat.hd <- readRDS(path.healthy.ref)
DefaultAssay(seurat.hd) <- "RNA"

vegs         <- VariableFeatures(seurat.hd)
feature.load <- seurat.hd@reductions$pca@feature.loadings[, 1:N_PCA_DIMS]

cat("Recovering scaling coefficients for", length(vegs), "variable genes...\n")
data.hd   <- seurat.hd@assays$RNA@data
meta.hd   <- seurat.hd@meta.data

coef.lm      <- list()
centers.vegs <- list()
sd.vegs      <- list()

for (gene in vegs) {
  expr       <- data.hd[gene, ]
  perc.mito  <- meta.hd$percent.mt
  numi       <- meta.hd$nCount_RNA
  model      <- lm(expr ~ perc.mito + numi)
  residuals  <- residuals(model)

  coef.lm[[gene]]      <- coef(model)
  centers.vegs[[gene]] <- mean(residuals)
  sd.vegs[[gene]]      <- sd(residuals)
}

coef.lm      <- do.call("rbind", coef.lm)
centers.vegs <- do.call("c",     centers.vegs)
sd.vegs      <- do.call("c",     sd.vegs)


# =============================================================================
# SECTION 2: Build (or recover) the healthy donor UMAP model
# =============================================================================
# Fits a UMAP model on the healthy donor PCA embeddings using uwot, retaining
# the model object so patient cells can be transformed into the same space.
# The random seed must match the original analysis for reproducibility.
#
# TIP: Save umap_model to disk after the first run to avoid re-fitting:
#   saveRDS(umap_model, file = file.path(dir.output, "umap_model_healthy.rds"))
# Then load it on subsequent runs:
#   umap_model <- readRDS(file.path(dir.output, "umap_model_healthy.rds"))

cat("Fitting healthy donor UMAP model...\n")
set.seed(42)
umap_model <- uwot::umap(
  as.matrix(seurat.hd@reductions$pca@cell.embeddings[, 1:N_PCA_DIMS]),
  n_neighbors = UMAP_N_NEIGHBORS,
  min_dist    = UMAP_MIN_DIST,
  ret_model   = TRUE,
  metric      = UMAP_METRIC
)
rownames(umap_model$embedding) <- rownames(seurat.hd@meta.data)

# Build healthy donor annotation and UMAP coordinate table for plotting
health_ann       <- data.table(seurat.hd@meta.data, keep.rownames = TRUE)
setkey(health_ann, rn)

umap.healthy     <- data.table(umap_model$embedding, keep.rownames = TRUE)
names(umap.healthy) <- c("rn", "refUMAP_1", "refUMAP_2")


# =============================================================================
# SECTION 3: Project patient samples onto healthy reference
# =============================================================================
# For each sample in the patient object:
#   1. Load the raw count matrix saved during preprocessing.
#   2. Apply healthy-donor-derived scaling (subtract regression fit, center,
#      scale) to the shared variable genes.
#   3. Project into PCA space via matrix multiplication with healthy loadings.
#   4. Embed into UMAP using umap_transform() on the saved model.
#   5. Assign each patient cell the cell type of its nearest healthy neighbor.
#   6. Save per-sample plot and cell-neighbor mapping table.

cat("Loading patient annotation object...\n")
patient.obj <- readRDS(path.patient.rds)
mllr_ann    <- data.table(patient.obj@meta.data, keep.rownames = TRUE)
setnames(mllr_ann, "rn", "V1")
setkey(mllr_ann, V1)

sample.ids <- unique(mllr_ann$sample.name)
cat("Samples to project:", length(sample.ids), "\n")

# Free healthy reference object memory — no longer needed as Seurat object
rm(seurat.hd)
gc()

cell.map <- NULL

for (sample0 in sample.ids) {
  cat("\n-- Projecting:", sample0, "--\n")

  # -- Load raw count matrix for this sample ---------------------------------
  mtx.path <- file.path(dir.matrices,
                        paste0("Matrix_related2_", sample0,
                               "_JointSeuratObj_doubletRemoved.rds"))
  if (!file.exists(mtx.path)) {
    warning("Matrix not found, skipping: ", mtx.path)
    next
  }
  mtx0 <- readRDS(mtx.path)

  # Filter annotation to this sample (all cells, or Blast only)
  if (PROJECT_ALL_CELLS) {
    mllr_ann1 <- mllr_ann[mllr_ann$sample.name == sample0, ]
  } else {
    mllr_ann1 <- mllr_ann[mllr_ann$sample.name == sample0 &
                           mllr_ann$Level1.anno == "Blast", ]
  }
  setkey(mllr_ann1, V1)

  # -- Compute mitochondrial fraction ----------------------------------------
  mito.genes <- grep("^MT-", rownames(mtx0), value = TRUE)
  perc.mito  <- Matrix::colSums(mtx0[mito.genes, ]) / Matrix::colSums(mtx0)

  seurat.tmp <- CreateSeuratObject(mtx0)
  seurat.tmp@meta.data$perc.mito <- perc.mito[colnames(seurat.tmp)]
  seurat.tmp <- NormalizeData(seurat.tmp)

  # -- Apply healthy-donor scaling to shared variable genes ------------------
  shared.genes  <- intersect(rownames(seurat.tmp), vegs)
  coef.lm1      <- coef.lm[shared.genes, ]
  centers.vegs1 <- centers.vegs[shared.genes]
  sd.vegs1      <- sd.vegs[shared.genes]

  mtx1       <- seurat.tmp@assays$RNA@data[shared.genes, ]
  meta1      <- seurat.tmp@meta.data
  cov.mat    <- cbind(1, meta1$perc.mito, meta1$nCount_RNA)

  # Scale: subtract regression fit, center on healthy mean, divide by healthy SD
  mtx1.scale <- (mtx1 - coef.lm1 %*% t(cov.mat) - centers.vegs1) / sd.vegs1

  # -- Project into PCA space ------------------------------------------------
  U    <- feature.load[shared.genes, 1:N_PCA_DIMS]
  pca1 <- t(mtx1.scale) %*% U

  # -- Embed into healthy UMAP -----------------------------------------------
  umap1             <- uwot::umap_transform(as.matrix(pca1), umap_model)
  rownames(umap1)   <- rownames(pca1)

  # Keep only barcodes present in the annotation table
  umap1 <- umap1[rownames(umap1) %in% mllr_ann1$V1, , drop = FALSE]

  # Attach UMAP coordinates to sample annotation
  umap1.dt       <- data.table(umap1, keep.rownames = TRUE)
  names(umap1.dt) <- c("rn", "refUMAP_1", "refUMAP_2")
  setkey(umap1.dt, rn)
  setkey(mllr_ann1, V1)
  mllr_ann1 <- cbind(mllr_ann1, umap1.dt[mllr_ann1$V1])

  # -- Projection plot: patient cells overlaid on healthy UMAP ---------------
  ctype1  <- mllr_ann1$Level1.anno
  n.ctype <- length(unique(ctype1))

  pdata <- rbind(
    data.table(umap_model$embedding, ctype = "Healthy"),
    data.table(umap1,                ctype = ctype1)
  )
  names(pdata)[1:2] <- c("UMAP_1", "UMAP_2")
  pdata$ctype <- factor(pdata$ctype,
                        levels = c("Healthy", sort(unique(ctype1))))

  plot.colors               <- colorRampPalette(brewer.pal(9, "Set1"))(n.ctype + 1)
  names(plot.colors)        <- c(sort(unique(ctype1)), "Healthy")
  plot.colors["Healthy"]    <- "#cccccc"

  pp <- ggplot(pdata, aes(x = UMAP_1, y = UMAP_2, col = ctype)) +
    geom_point(size = 0.2) +
    scale_colour_manual(name = "Cell type", values = plot.colors) +
    ggtitle(paste("Projection:", sample0)) +
    theme_classic()

  ggsave(pp, device = "eps", width = 9, height = 7,
         filename = file.path(dir.plots,
                              paste0("proj_direct_", sample0,
                                     "_toHealthyRef.eps")))

  # -- Assign nearest healthy neighbor cell type -----------------------------
  # Computes pairwise Euclidean distances between patient and healthy UMAP
  # coordinates, then assigns each patient cell the label of its closest
  # healthy neighbor.
  umap.combined <- rbind(umap_model$embedding, umap1)
  dist.mat      <- as.matrix(dist(umap.combined))

  dist.sub      <- dist.mat[rownames(umap1), rownames(umap_model$embedding)]
  nearest.hd    <- sapply(seq_len(nrow(umap1)),
                          function(i) names(which.min(dist.sub[i, ])))

  cell.map.tmp <- data.table(
    patient_bc  = rownames(umap1),
    healthy_bc  = nearest.hd
  )
  setkey(cell.map.tmp, healthy_bc)
  cell.map.tmp[, ctype_healthy  := health_ann[J(cell.map.tmp$healthy_bc)]$predicted.celltype]
  cell.map.tmp[, sample_healthy := health_ann[J(cell.map.tmp$healthy_bc)]$sample.name]

  setkey(cell.map.tmp, patient_bc)
  setkey(mllr_ann1, V1)
  cell.map.tmp[, ctype_patient  := mllr_ann1[J(cell.map.tmp$patient_bc)]$Level1.anno]
  cell.map.tmp[, sample_patient := mllr_ann1[J(cell.map.tmp$patient_bc)]$sample.name]
  cell.map.tmp[, refUMAP_1      := mllr_ann1[J(cell.map.tmp$patient_bc)]$refUMAP_1]
  cell.map.tmp[, refUMAP_2      := mllr_ann1[J(cell.map.tmp$patient_bc)]$refUMAP_2]

  # Save per-sample mapping table
  write.table(cell.map.tmp,
              file      = file.path(dir.metadata.sample,
                                    paste0(sample0, "_proj2healthyCells.csv")),
              sep       = "\t", row.names = FALSE, quote = FALSE)

  cell.map <- rbind(cell.map, cell.map.tmp)
  cat("  Projected", nrow(umap1), "cells.\n")
}

# Save combined mapping table across all samples
write.table(cell.map, file = path.output.combined,
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nProjection complete.", nrow(cell.map), "cells across",
    length(sample.ids), "samples.\n")
cat("Combined table saved to:", path.output.combined, "\n")
