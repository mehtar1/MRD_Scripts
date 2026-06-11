# =============================================================================
# CITE-Seq Preprocessing and Integration Pipeline
# =============================================================================
#
# Description:
#   Processes raw CITE-Seq CellRanger output (RNA + ADT) from multiple samples
#   into a single, batch-corrected Seurat object ready for annotation.
#   Downstream annotation and inferCNV integration are handled in
#   02_citeseq_annotate.R.
#
#   Pipeline steps:
#     1. Per-sample loading  — reads h5 matrices, creates RNA + ADT assays,
#                              computes QC metrics (mt%, ribo%, hb%)
#     2. Per-sample QC plots — violin and UMAP plots for threshold inspection
#                              (uncomment pdf block to generate)
#     3. Per-sample filtering — applies sample-specific thresholds defined
#                               in SAMPLE_FILTERS; produces pre-processed objects
#     4. Cell cycle + heat shock scoring — for regression in ScaleData
#     5. Per-sample healthy reference projection — assigns celltype.short and
#                                                   celltype.binned labels
#     6. Doublet detection — flags doublets for removal at merge step
#     7. Merge + 2nd-pass filter — removes doublets and low-quality clusters
#                                   identified visually after initial merge
#     8. Cohort merge — combines new samples with pilot rounds 1 & 2
#     9. RPCA integration — corrects batch effects across patient samples
#    10. Post-integration projection — refreshes healthy reference labels on
#                                       the integrated embedding
#    11. Barcode export — writes per-sample barcode TSV files
#
# Inputs (set in USER CONFIGURATION):
#   - CellRanger h5 files (one per sample)
#   - Pilot rounds 1 & 2 combined Seurat RDS
#   - Healthy donor reference RDS
#   - Cell cycle gene list (Regev lab format)
#   - Heat shock gene list
#
# Outputs:
#   - filtered per-sample objects RData (intermediate)
#   - Round3 merged + filtered RDS (intermediate)
#   - Final merged + integrated Seurat RDS → input for 02_citeseq_annotate.R
#   - Per-sample barcode TSV files
#
# Dependencies:
#   Seurat, dplyr, ggplot2, tidyverse, patchwork, dsb, ComplexHeatmap,
#   data.table, RColorBrewer, RSpectra, uwot
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(dsb)
library(data.table)
library(RColorBrewer)
library(RSpectra)
library(uwot)

options(future.globals.maxSize = 100 * 1000 * 1024^2)
`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# CellRanger h5 output files — one per sample, in the same order as SAMPLE_NAMES
CELLRANGER_FILES <- c(
  "/path/to/4526R_CITESeq/outs/filtered_feature_bc_matrix.h5",
  "/path/to/5481_CITESeq/outs/filtered_feature_bc_matrix.h5",
  "/path/to/5850_CITESeq/outs/filtered_feature_bc_matrix.h5",
  "/path/to/5907_CITESeq/outs/filtered_feature_bc_matrix.h5",
  "/path/to/6024_CITESeq/outs/filtered_feature_bc_matrix.h5",
  "/path/to/6146_CITESeq/outs/filtered_feature_bc_matrix.h5"
)
SAMPLE_NAMES <- c("4526R", "5481", "5850", "5907", "6024", "6146")

# Per-sample QC filtering thresholds.
# Each entry must match a name in SAMPLE_NAMES.
# Adjust min/max values based on per-sample VlnPlot inspection (Section 2).
SAMPLE_FILTERS <- list(
  "4526R" = list(nFeature_RNA_min = 500,  nFeature_RNA_max = 3000,
                 nCount_ADT_max   = 10000, nCount_RNA_max   = 30000,
                 percent_mt_max   = 20,    percent_hb_max   = 1),
  "5481"  = list(nFeature_RNA_min = 500,  nFeature_RNA_max = 4500,
                 nCount_ADT_max   = 10000, nCount_RNA_max   = 30000,
                 percent_mt_max   = 20,    percent_hb_max   = 1),
  "5850"  = list(nFeature_RNA_min = 300,  nFeature_RNA_max = 4500,
                 nCount_ADT_max   = 10000, nCount_RNA_max   = 30000,
                 percent_mt_max   = 20,    percent_hb_max   = 1),
  "5907"  = list(nFeature_RNA_min = 500,  nFeature_RNA_max = 10000,
                 nCount_ADT_max   = 15000, nCount_RNA_max   = 30000,
                 percent_mt_max   = 20,    percent_hb_max   = 1),
  "6024"  = list(nFeature_RNA_min = 500,  nFeature_RNA_max = 4500,
                 nCount_ADT_max   = 15000, nCount_RNA_max   = 30000,
                 percent_mt_max   = 20,    percent_hb_max   = 1),
  "6146"  = list(nFeature_RNA_min = 500,  nFeature_RNA_max = 4500,
                 nCount_ADT_max   = 10000, nCount_RNA_max   = 30000,
                 percent_mt_max   = 20,    percent_hb_max   = 1)
)

# Cluster IDs to remove after initial unintegrated merge (low-quality clusters
# identified by high mt% or hb% on VlnPlot — update after visual inspection)
ROUND3_REMOVE_CLUSTERS  <- c(0, 1, 2, 21, 25)   # from RNA_snn_res default

# Cluster IDs to remove after cohort merge (update after visual inspection)
COHORT_REMOVE_CLUSTERS  <- c(11, 17)              # from RNA_snn_res.0.4

# Samples to exclude from RPCA integration (failed QC or redundant runs)
EXCLUDE_FROM_INTEGRATION <- c("6099", "4526")

# Annotation cluster IDs (RPCA clusters at res 1.4) assigned before inferCNV
# Update these after inspecting the integrated UMAP
LEVEL1_CLUSTERS <- list(
  Myeloid     = c(28, 32, 15, 33, 29, 19),
  T_NK        = c(12, 13, 20, 27, 17, 30, 21),
  B           = c(11, 26, 14),
  Ambiguous_1 = c(16),
  Ambiguous_2 = c(24)
  # All other clusters → "Blast" by default
)

# Patient-to-sample mapping (Sample.Names → patient ID)
# "R" suffix samples are renamed before matching (e.g. "6099R" → "6099")
PATIENT_MAP <- list(
  PT1 = c("4526", "5481"),
  PT2 = c("5995", "6024"),
  PT3 = c("6092", "6061"),
  PT4 = c("6099", "6146"),
  PT5 = c("5850", "5907")
)

# IDX vs MRD assignment by sample name
MRD_SAMPLES <- c("5481", "6024", "6061", "6146", "5907")   # all others → IDX

# Reference files
path.healthy.ref   <- "/path/to/HealthyReference/full.ref.gini.genes.cc.removed.rds"
path.cell.cycle    <- "/path/to/References/regev_lab_cell_cycle_genes.txt"
path.heat.shock    <- "/path/to/References/heat_shock_geneList.txt"
path.pilot.rds     <- "/path/to/Objects/pilotRound1And2Combined_mergedFiltered.rds"

# Output paths
dir.plots          <- "/path/to/Plots/Sample_Plots/"
dir.barcodes       <- "/path/to/Tables/cell_barcodes/"
path.filtered.rdata <- "/path/to/Objects/sequencing_r3_filteredObjs.Rdata"
path.round3.rds    <- "/path/to/Objects/Round3_mergedFiltered.rds"
path.output.rds    <- "/path/to/Objects/5Pairs_IntegratedAnnotated.rds"

# PCA dimensions
N_PCA <- 30
N_PCA_INTEGRATION <- 40


# ── Shared reference gene lists (loaded once) ─────────────────────────────────

hd.ref         <- readRDS(path.healthy.ref)
cycle.genes    <- fread(path.cell.cycle, header = FALSE)$V1
s.genes        <- cycle.genes[1:43]
g2m.genes      <- cycle.genes[44:97]
heat.shock.genes <- fread(path.heat.shock)$`Approved symbol`


# =============================================================================
# SECTION 1: Per-sample loading and QC metrics
# =============================================================================
# Reads each CellRanger h5 file, creates a Seurat object with RNA and ADT
# assays, and computes mitochondrial, ribosomal, and hemoglobin percentages.

cat("Loading", length(SAMPLE_NAMES), "samples...\n")

sample.list <- lapply(seq_along(CELLRANGER_FILES), function(i) {
  mtx       <- Read10X_h5(CELLRANGER_FILES[i])
  s         <- CreateSeuratObject(counts = mtx$`Gene Expression`)
  s[["ADT"]] <- CreateAssay5Object(counts = mtx$`Antibody Capture`)

  s <- NormalizeData(s, assay = "RNA")
  s <- NormalizeData(s, assay = "ADT", normalization.method = "CLR")

  s[["percent.mt"]]   <- PercentageFeatureSet(s, pattern = "^MT-")
  s <- PercentageFeatureSet(s, "^RP[SL]",        col.name = "percent.ribo")
  s <- PercentageFeatureSet(s, "^HB[^(P|E|S)]",  col.name = "percent.hb")
  s$Sample <- SAMPLE_NAMES[i]
  s
})
names(sample.list) <- SAMPLE_NAMES


# =============================================================================
# SECTION 2: Per-sample QC plots (optional — uncomment to generate)
# =============================================================================
# Generates unfiltered VlnPlots and UMAPs for setting per-sample thresholds.
# Run this section first, then adjust SAMPLE_FILTERS in USER CONFIGURATION.

sample.list <- lapply(sample.list, function(s) {
  NormalizeData(s, assay = "RNA", normalization.method = "LogNormalize",
                scale.factor = 10000) %>%
    FindVariableFeatures(assay = "RNA", selection.method = "vst",
                         nfeatures = 2000) %>%
    ScaleData(assay = "RNA", verbose = FALSE) %>%
    RunPCA(npcs = N_PCA, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
    FindClusters(resolution = c(0.4, 0.8), verbose = FALSE)
})

# -- Uncomment to generate per-sample unfiltered QC PDFs ---------------------
# for (sname in SAMPLE_NAMES) {
#   s <- sample.list[[sname]]
#   pdf(file.path(dir.plots, paste0(sname, "_QC_unfiltered.pdf")),
#       width = 8, height = 6)
#   print(VlnPlot(s, group.by = "orig.ident", pt.size = 0, log = TRUE,
#                 features = c("nCount_ADT", "nFeature_ADT", "percent.ribo",
#                              "nCount_RNA", "nFeature_RNA", "percent.mt",
#                              "percent.hb"), ncol = 4))
#   print(DimPlot(s, reduction = "umap", group.by = "seurat_clusters",
#                 label = TRUE) + ggtitle(paste(sname, "Unfiltered")))
#   print(FeaturePlot(s, features = c("nFeature_RNA", "nFeature_ADT",
#                                     "nCount_RNA", "nCount_ADT"),
#                     order = TRUE, min.cutoff = "q1",
#                     max.cutoff = "q99") & NoAxes())
#   dev.off()
# }


# =============================================================================
# SECTION 3: Per-sample filtering
# =============================================================================
# Applies the thresholds from SAMPLE_FILTERS to each sample. Each filter is
# expressed separately so the thresholds are visible and auditable.

cat("Filtering samples...\n")

pre.processed.objs <- lapply(SAMPLE_NAMES, function(sname) {
  s   <- sample.list[[sname]]
  flt <- SAMPLE_FILTERS[[sname]]

  # Diagnostic: highlight cells near thresholds before filtering
  print(VlnPlot(s, group.by = "orig.ident", pt.size = 0, log = TRUE,
                flip = TRUE, ncol = 4,
                features = c("nCount_ADT", "nFeature_ADT", "percent.ribo",
                             "nCount_RNA",  "nFeature_RNA", "percent.mt",
                             "percent.hb")) +
          ggtitle(sname))

  subset(s,
         subset = nFeature_RNA > flt$nFeature_RNA_min &
                  nFeature_RNA < flt$nFeature_RNA_max &
                  percent.mt   < flt$percent_mt_max   &
                  nCount_ADT   < flt$nCount_ADT_max   &
                  nCount_RNA   < flt$nCount_RNA_max   &
                  percent.hb   < flt$percent_hb_max)
})
names(pre.processed.objs) <- SAMPLE_NAMES


# =============================================================================
# SECTION 4: Cell cycle + heat shock scoring, normalization, and projection
# =============================================================================
# For each filtered sample:
#   - Scores cell cycle phase (S, G2M) and heat shock module
#   - Runs PCA and UMAP with confound regression
#   - Projects onto healthy donor reference to assign celltype.short labels
#   - Detects doublets with DoubletFinder

cat("Scoring, normalizing, and projecting", length(SAMPLE_NAMES), "samples...\n")

filtered.r1.objs <- lapply(SAMPLE_NAMES, function(sname) {
  s <- pre.processed.objs[[sname]]

  # Cell cycle and heat shock scores
  s <- CellCycleScoring(s, s.features = s.genes,
                        g2m.features = g2m.genes, set.ident = FALSE)
  s <- AddModuleScore(s, features = list(heat.shock.genes),
                      name = "HeatShock.Score")

  # Normalize, scale (regressing confounders), reduce
  s <- NormalizeData(s, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 10000) %>%
    FindVariableFeatures(assay = "RNA", selection.method = "vst",
                         nfeatures = 2000) %>%
    ScaleData(assay = "RNA",
              vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"),
              verbose = FALSE) %>%
    RunPCA(npcs = N_PCA, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
    FindNeighbors(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
    FindClusters(resolution = c(0.2, 0.4, 0.6, 0.7, 0.8, 1), verbose = FALSE)

  # Project onto healthy donor reference
  anchors     <- FindTransferAnchors(reference = hd.ref, query = s,
                                      dims = 1:N_PCA,
                                      reference.reduction = "pca")
  predictions <- TransferData(
    anchorset = anchors,
    list(celltype = hd.ref$cell.type.short,
         celltype.bin = hd.ref$cell.type.binned),
    dims = 1:N_PCA
  )
  predictions <- lapply(predictions, function(x) {
    dplyr::select(x, predicted.id, prediction.score.max)
  })
  predictions         <- cbind(predictions[[1]], predictions[[2]])
  colnames(predictions) <- c("celltype.short", "max.score.celltype.short",
                              "celltype.binned", "max.score.celltype.binned")
  s <- AddMetaData(s, metadata = predictions)

  # Doublet detection (FindDoublets wraps DoubletFinder)
  s <- FindDoublets(s)
  s
})
names(filtered.r1.objs) <- SAMPLE_NAMES

# Rename cells to avoid barcode collisions on merge
for (sname in SAMPLE_NAMES) {
  filtered.r1.objs[[sname]] <- RenameCells(filtered.r1.objs[[sname]],
                                            add.cell.id = sname)
}

save(filtered.r1.objs, file = path.filtered.rdata)
cat("Per-sample filtered objects saved.\n")


# =============================================================================
# SECTION 5: Merge round-3 samples and second-pass filtering
# =============================================================================
# Merges all round-3 filtered samples, removes doublets, then re-clusters to
# identify any remaining low-quality clusters for removal (set in
# ROUND3_REMOVE_CLUSTERS after visual inspection).

load(path.filtered.rdata)

merged.r3 <- merge(filtered.r1.objs[[1]], filtered.r1.objs[2:6],
                   add.cell.ids = SAMPLE_NAMES)
merged.r3 <- subset(merged.r3, subset = Doublet_Singlet != "Doublet")

merged.r3 <- NormalizeData(merged.r3, assay = "RNA",
                            normalization.method = "LogNormalize",
                            scale.factor = 10000) %>%
  FindVariableFeatures(assay = "RNA", selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(assay = "RNA",
            vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"),
            verbose = FALSE) %>%
  RunPCA(npcs = N_PCA, verbose = FALSE, assay = "RNA") %>%
  RunUMAP(assay = "RNA", reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
  FindClusters(resolution = c(0.2, 0.4, 0.6, 0.7, 0.8, 1), verbose = FALSE)

# Inspect these plots to confirm ROUND3_REMOVE_CLUSTERS before filtering
VlnPlot(merged.r3, group.by = "seurat_clusters",
        features = c("nCount_RNA", "percent.mt", "percent.hb"),
        stack = TRUE, flip = TRUE) +
  DimPlot(merged.r3, group.by = "seurat_clusters", label = TRUE,
          repel = TRUE) & NoLegend()

merged.r3 <- subset(merged.r3,
                    subset = seurat_clusters %notin% ROUND3_REMOVE_CLUSTERS)

saveRDS(merged.r3, path.round3.rds)
cat("Round 3 merged object saved.\n")


# =============================================================================
# SECTION 6: Cohort merge with pilot rounds 1 & 2
# =============================================================================
# Combines round-3 samples with the previously processed pilot cohort, then
# applies a second round of cluster-based QC filtering.

pilot.objs <- readRDS(path.pilot.rds)
merged.all <- merge(merged.r3, pilot.objs)

merged.all <- NormalizeData(merged.all, assay = "RNA",
                             normalization.method = "LogNormalize",
                             scale.factor = 10000) %>%
  FindVariableFeatures(assay = "RNA", selection.method = "vst",
                       nfeatures = 2000) %>%
  ScaleData(assay = "RNA",
            vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"),
            verbose = FALSE) %>%
  RunPCA(npcs = N_PCA, verbose = FALSE, assay = "RNA") %>%
  RunUMAP(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:N_PCA, verbose = FALSE) %>%
  FindClusters(resolution = c(0.2, 0.4, 0.6, 0.7, 0.8, 1), verbose = FALSE)

# Inspect to confirm COHORT_REMOVE_CLUSTERS
VlnPlot(merged.all, group.by = "RNA_snn_res.0.6",
        features = c("nCount_RNA", "percent.mt", "percent.hb"),
        stack = TRUE, flip = TRUE) +
  DimPlot(merged.all, group.by = "RNA_snn_res.0.6", label = TRUE,
          repel = TRUE) & NoLegend()

merged.all <- subset(merged.all,
                     subset = RNA_snn_res.0.4 %notin% COHORT_REMOVE_CLUSTERS)


# =============================================================================
# SECTION 7: Assign clinical metadata
# =============================================================================
# Adds Stage (IDX/MRD), patient ID, and composite grouping columns.
# All mappings are driven by USER CONFIGURATION — no hardcoded sample lists.

# Build reverse lookup: sample name → patient ID
sample_to_patient <- setNames(
  rep(names(PATIENT_MAP), lengths(PATIENT_MAP)),
  unlist(PATIENT_MAP)
)

merged.all@meta.data <- mutate(
  merged.all@meta.data,
  # Normalize "R"-suffix sample names for matching (e.g. "6099R" → "6099")
  Sample.Names = case_when(
    Sample == "6099R" ~ "6099",
    Sample == "4526R" ~ "4526",
    .default = as.character(Sample)
  ),
  Stage = case_when(
    Sample.Names %in% MRD_SAMPLES ~ "MRD",
    .default = "IDX"
  )
)

merged.all$patient.sample <- sample_to_patient[merged.all$Sample.Names]
merged.all@meta.data <- mutate(
  merged.all@meta.data,
  Patient.Stage.All       = paste0(patient.sample, "_", Stage),
  patient.sample.Stage    = paste0(patient.sample, "_", Stage)
)


# =============================================================================
# SECTION 8: RPCA integration
# =============================================================================
# Corrects batch effects across patient samples using reciprocal PCA (RPCA).
# Samples in EXCLUDE_FROM_INTEGRATION are dropped before integration
# (typically early pilot runs superseded by repeat samples).

merged.all <- subset(merged.all,
                     subset = Sample %notin% EXCLUDE_FROM_INTEGRATION)
merged.all <- JoinLayers(merged.all)
merged.all[["RNA"]] <- split(merged.all[["RNA"]],
                              f = merged.all$patient.sample)

merged.all <- FindVariableFeatures(merged.all) %>%
  ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt",
                                 "percent.ribo", "HeatShock.Score1")) %>%
  RunPCA()

merged.all <- IntegrateLayers(
  object          = merged.all,
  method          = RPCAIntegration,
  orig.reduction  = "pca",
  new.reduction   = "integrated.rpca",
  verbose         = FALSE
)

merged.all <- FindNeighbors(merged.all, reduction = "integrated.rpca",
                             dims = 1:N_PCA_INTEGRATION) %>%
  FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4),
               cluster.name = c("rpca_clusters_0.2", "rpca_clusters_0.4",
                                 "rpca_clusters_0.6", "rpca_clusters_0.8",
                                 "rpca_clusters_1.0", "rpca_clusters_1.2",
                                 "rpca_clusters_1.4")) %>%
  RunUMAP(reduction = "integrated.rpca", dims = 1:N_PCA_INTEGRATION,
          reduction.name = "umap.rpca.pc40")

DimPlot(merged.all, group.by = "Sample",
        reduction = "umap.rpca.pc40", shuffle = TRUE)


# =============================================================================
# SECTION 9: Post-integration healthy reference projection
# =============================================================================
# Re-projects each sample onto the healthy donor reference using the
# integrated embedding. Overwrites the per-sample celltype labels from
# Section 4 with ones computed on the batch-corrected data.

cat("Re-projecting onto healthy reference post-integration...\n")

sample.objs      <- SplitObject(merged.all, split.by = "Sample")
projections.list <- lapply(sample.objs, function(s) {
  VariableFeatures(s) <- VariableFeatures(hd.ref)
  s <- NormalizeData(s, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 10000) %>%
    ScaleData(assay = "RNA",
              vars.to.regress = c("nCount_RNA", "percent.mt", "percent.ribo"),
              verbose = FALSE) %>%
    RunPCA(npcs = N_PCA, verbose = FALSE) %>%
    RunUMAP(reduction = "pca", dims = 1:N_PCA, verbose = FALSE)

  anchors     <- FindTransferAnchors(reference = hd.ref, query = s,
                                      dims = 1:N_PCA,
                                      reference.reduction = "pca")
  predictions <- TransferData(
    anchorset = anchors,
    list(celltype = hd.ref$cell.type.short,
         celltype.bin = hd.ref$cell.type.binned),
    dims = 1:N_PCA
  )
  predictions <- lapply(predictions, function(x) {
    dplyr::select(x, predicted.id, prediction.score.max)
  })
  predictions         <- cbind(predictions[[1]], predictions[[2]])
  colnames(predictions) <- c("celltype.short", "max.score.celltype.short",
                              "celltype.binned", "max.score.celltype.binned")
  predictions
})

projections.merged <- do.call(rbind, projections.list)

# Remove old celltype columns before adding updated ones
merged.all@meta.data <- merged.all@meta.data[,
  !grepl("celltype\\.", colnames(merged.all@meta.data))]
merged.all <- AddMetaData(merged.all, projections.merged)

saveRDS(merged.all, path.output.rds)
cat("Integrated object saved to:", path.output.rds, "\n")
cat("Hand off to 02_citeseq_annotate.R for Level1 annotation.\n")


# =============================================================================
# SECTION 10: Barcode export (utility)
# =============================================================================
# Writes per-sample barcode TSV files (e.g. for use with velocyto or pySCENIC).
# Strips the sample-prefix added by RenameCells to recover original barcodes.

cat("Exporting barcodes...\n")
dir.create(dir.barcodes, recursive = TRUE, showWarnings = FALSE)

merged.split <- SplitObject(merged.all, split.by = "Sample")
for (sname in names(merged.split)) {
  barcodes <- colnames(merged.split[[sname]])
  # Strip "sampleName_" prefix added during merge
  barcodes <- sapply(barcodes,
                     function(x) strsplit(x, "_")[[1]][2], USE.NAMES = FALSE)
  cat(barcodes,
      file = file.path(dir.barcodes, paste0(sname, "_barcodes.tsv")),
      sep = "\n")
}
cat("Barcodes exported to:", dir.barcodes, "\n")
