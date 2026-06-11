# =============================================================================
# CITE-Seq Blast Annotation, inferCNV Integration, and ADT Visualization
# =============================================================================
#
# Description:
#   Annotates leukemic blast populations in the integrated CITE-Seq cohort
#   using cluster-based Level1 annotation, inferCNV copy-number evidence,
#   and healthy donor projection labels. Generates all downstream ADT protein
#   heatmaps, RNA-ADT correlation plots, supplementary UMAPs, and volcano
#   plots for marker comparisons.
#
#   Prerequisite: run 01_citeseq_preprocess.R first to generate the
#   integrated Seurat object at path.integrated.rds.
#
#   Pipeline steps:
#     0. Initial Level1 annotation — assigns Blast / T_NK / B / Myeloid from
#                                     RPCA cluster IDs; handles ambiguous clusters
#     1. ADT heatmaps and RNA-ADT correlation (pre-inferCNV)
#     2. Per-pair inferCNV re-clustering — merges IDX+MRD per patient pair,
#                                          re-clusters, re-assigns Level1.anno
#                                          using CNV evidence
#     3. CNV metadata transfer — adds CNV scores to merged integrated object
#                                 and rebuilds Blast_Lineage labels
#     4. ADT heatmaps post-inferCNV (by lineage, stage, sample)
#     5. Supplementary UMAP plots
#     6. ADT differential expression and volcano plots
#
# Inputs (set in USER CONFIGURATION):
#   - Integrated Seurat RDS (output of 01_citeseq_preprocess.R)
#   - Per-sample inferCNV Seurat RDS files (output of 03_citeseq_infercnv.R)
#
# Outputs:
#   - Annotated Seurat RDS (with Level1.anno, Blast_Lineage, CNV metadata)
#   - ADT heatmap PDFs
#   - RNA-ADT correlation PDFs
#   - Supplementary UMAP PDF
#   - Volcano plot PDF
#
# Dependencies:
#   Seurat, dplyr, ggplot2, tidyverse, patchwork, dsb, ComplexHeatmap,
#   data.table, RColorBrewer, RSpectra, uwot, ggpubr, EnhancedVolcano
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(dsb)
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(RSpectra)
library(uwot)
library(ggpubr)
library(EnhancedVolcano)

options(future.globals.maxSize = 100 * 1000 * 1024^2)
`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Input: integrated object from 01_citeseq_preprocess.R
path.integrated.rds <- "/path/to/Objects/5Pairs_IntegratedAnnotated.rds"

# Per-sample inferCNV Seurat RDS files (output of 03_citeseq_infercnv.R).
# Keys are sample IDs; values are file paths.
infercnv.paths <- list(
  "6099R" = "/path/to/Objects/6099R_InferCNV.rds",
  "6061"  = "/path/to/Objects/6061_InferCNV.rds",
  "6146"  = "/path/to/Objects/6146_InferCNV.rds",
  "6092"  = "/path/to/Objects/6092_InferCNV.rds"
)

# RPCA cluster IDs for initial Level1 annotation (res 1.4).
# Update after inspecting the integrated UMAP — these are the cluster numbers
# that correspond to each normal cell type; all others become "Blast".
LEVEL1_CLUSTERS <- list(
  Myeloid     = c(28, 32, 15, 33, 29, 19),
  T_NK        = c(12, 13, 20, 27, 17, 30, 21),
  B           = c(11, 26, 14),
  Ambiguous_1 = c(16),    # inspect further; likely remove or relabel
  Ambiguous_2 = c(24)     # inspect further; likely remove or relabel
)

# RNA_snn_res.1 cluster IDs for final Level1.anno after ambiguity review.
# Set after inspecting VlnPlot and DotPlot in Section 0.
LEVEL1_FINAL_CLUSTERS <- list(
  Myeloid = c(17, 21, 29, 27),
  T_NK    = c(11, 12, 18, 24, 25),
  B       = c(10, 13)
)

# Patient pair config for per-pair inferCNV re-clustering.
# clusters.myeloid / .tnk / .b: cluster IDs from merged_res_1.0 to label as
# normal; all others → Blast. CNV features: chromosomes with recurrent events.
pair.config <- list(
  pair3 = list(
    idx.sample       = "6092",
    mrd.sample       = "6061",
    clusters.myeloid = c(10, 16, 17),
    clusters.tnk     = c(11, 15),
    clusters.b       = c(18),
    cnv.features     = c("proportion_dupli_chr4",  "proportion_dupli_chr8",
                          "proportion_dupli_chr21", "proportion_dupli_chr22")
  ),
  pair4 = list(
    idx.sample       = "6099R",
    mrd.sample       = "6146",
    clusters.myeloid = c(13),
    clusters.tnk     = c(6, 9, 14, 15),
    clusters.b       = c(7, 8),
    cnv.features     = c("proportion_dupli_chr21", "proportion_dupli_chrX",
                          "proportion_loss_chrX")
  )
  # Add further pairs here following the same structure
)

# CNV metadata columns to transfer from per-pair objects to the merged object
METADATA_CNV_COLS <- c(
  "has_dupli_chr21",  "proportion_cnv_chr21",  "proportion_dupli_chr21",
  "has_dupli_chr22",  "proportion_cnv_chr22",  "proportion_dupli_chr22",
  "has_dupli_chrX",   "proportion_cnv_chrX",   "proportion_dupli_chrX",
  "has_loss_chrX",    "proportion_loss_chrX",
  "has_dupli_chr4",   "proportion_cnv_chr4",   "proportion_dupli_chr4",
  "Level1.anno"
)

# Cell type groupings for blast lineage assignment
HSPC_CTYPES  <- c("HSPC", "LMPP")
MY_CTYPES    <- c("GMP", "MEP", "DC-Progenitor", "pDC", "cDC", "Mono")
B_CTYPES     <- c("CLP", "Pre-pro-B", "Pro-B", "Pre-B", "Immature-B", "Mature-B")
T_CTYPES     <- c("Pro-T", "Pre-T", "CyclingDP", "Naive_T", "alpha-beta",
                   "alpha-beta(mature)", "EffectorT", "gamma-delta", "NK")

BLAST_LINEAGE_LEVELS <- c("B-like", "Myeloid-like", "HSPC-like",
                           "T-cell-like", "B", "T_NK", "Myeloid")
LINEAGE_COLS         <- c("#377eb8", "#e41a1c", "#4daf4a", "#984ea3",
                           "#cbd5e8", "#fdcdac", "#f4cae4")

# ADT panel for heatmaps
CORE_ADT_PANEL <- c("CD33.1", "CD371", "CD36.1", "CD305",
                     "CD3",    "CD20",  "CD10",   "CD24.1", "CD19.1", "CD45",
                     "CD72.1", "CD123", "CD64",   "CD83.1", "CD13",   "CD117",
                     "CD38.1", "CD34.1")

# ADT-RNA gene pairs for correlation plots (ADT name → RNA gene symbol)
# Remove any pairs where the gene is absent from your dataset
ADT_RNA_PAIRS <- data.frame(
  ADT = c("CD3",    "CD56.1", "CD19.1", "CD33.1", "CD11c.1","CD123.1","CD4.1",
           "CD16.1", "CD25.1", "CD20",   "HLA-DR", "CD11b",  "CD38.1", "CD127",
           "CD10",   "CD90",   "CD117",  "CD34.1", "CD22.1", "CD24.1",
           "CD36.1", "CD64",   "CD13",   "CD14",   "CD45",   "CD305",  "CD99.1",
           "CD163",  "CD105",  "CD98",   "ATP1B3", "CD119",  "CD83.1", "CD9.1",
           "CD72.1", "CD371",  "FLT3.1"),
  RNA = c("CD3E",   "NCAM1",  "CD19",   "CD33",   "ITGAX",  "IL3RA",  "CD4",
           "FCGR3A", "IL2RA",  "MS4A1",  "HLA-DRA","ITGAM",  "CD38",   "IL7R",
           "MME",    "THY1",   "KIT",    "CD34",   "CD22",   "CD24",
           "CD36",   "FCGR1A", "ANPEP",  "CD14",   "PTPRC",  "LAIR1",  "CD99",
           "MSR1",   "ENG",    "SLC7A5", "ATP1B3", "IFNGR1", "CD83",   "CD9",
           "CD72",   "CLEC12A","FLT3")
)

# Samples to exclude from ADT correlation batch annotation
ADT_BATCH_MAP <- list(
  ADT_1 = c("6099", "5995"),
  ADT_2 = c("6099R", "4526", "6092", "6061")
  # all others → ADT_3
)

# Subtype assignment by patient sample ID
SUBTYPE_MAP <- list(KMT2Ar = c("PT3", "PT5"))

# DE cutoffs for volcano plots
DE_P_CUT  <- 0.05
DE_FC_CUT <- 0.25

# Output paths
dir.plots      <- "/path/to/Plots/Sample_Plots/"
path.output.rds <- "/path/to/Objects/5Pairs_Integrated_inferCNVAnno.rds"


# ── Helper functions ──────────────────────────────────────────────────────────

# Assign blast lineage labels using healthy donor projection cell types.
assign_blast_lineage <- function(seurat_obj) {
  seurat_obj@meta.data <- mutate(
    seurat_obj@meta.data,
    Blast_Lineage = case_when(
      Level1.anno == "Blast" & celltype.short %in% HSPC_CTYPES ~ "HSPC-like",
      Level1.anno == "Blast" & celltype.short %in% MY_CTYPES   ~ "Myeloid-like",
      Level1.anno == "Blast" & celltype.short %in% B_CTYPES    ~ "B-like",
      Level1.anno == "Blast" & celltype.short %in% T_CTYPES    ~ "T-cell-like",
      .default = Level1.anno
    )
  )
  seurat_obj$Blast_Lineage <- factor(seurat_obj$Blast_Lineage,
                                      levels = BLAST_LINEAGE_LEVELS)
  seurat_obj
}

# ADT median expression heatmap grouped by a metadata variable.
adt_heatmap <- function(seurat_obj,
                        prots       = rownames(seurat_obj@assays$ADT@features)[1:37],
                        groupByVar,
                        clusterIdents = TRUE,
                        title         = NULL) {
  require(ComplexHeatmap)
  s         <- JoinLayers(seurat_obj, assay = "ADT")
  s         <- NormalizeData(s, assay = "ADT", normalization.method = "CLR")
  cite_data <- t(GetAssayData(s, slot = "data", assay = "ADT"))
  d         <- cbind(s@meta.data, cite_data)

  dat_plot <- d %>%
    dplyr::group_by(.data[[groupByVar]]) %>%
    dplyr::summarize_at(.vars = prots, .funs = median) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(groupByVar)

  prot_col <- circlize::colorRamp2(
    breaks = seq(0.25, 2, by = 0.0625),
    colors = viridis::viridis(n = 29, option = "B")
  )
  p <- Heatmap(t(dat_plot)[prots, ],
               name            = "protein",
               col             = prot_col,
               use_raster      = TRUE,
               cluster_columns = clusterIdents,
               row_names_gp    = gpar(color = "black", fontsize = 10),
               column_title    = title,
               column_title_gp = gpar(fontsize = 12, col = "black"))
  print(p)
  invisible(p)
}

# Stacked bar chart of cell type composition per sample/group.
cluster.plot <- function(table, x.var, y.var,
                          return.table = FALSE, color.scale = NULL) {
  freq <- table %>%
    dplyr::select({{ x.var }}, {{ y.var }}) %>%
    group_by({{ x.var }}) %>% mutate(total = n()) %>%
    group_by({{ x.var }}, {{ y.var }}) %>%
    mutate(freq = n() / total) %>%
    dplyr::select({{ x.var }}, {{ y.var }}, freq) %>% distinct()

  totals <- table %>%
    dplyr::select({{ x.var }}) %>%
    group_by({{ x.var }}) %>% mutate(total_n = n()) %>% unique()

  p <- ggplot(freq) +
    geom_col(aes(x = {{ x.var }}, y = freq, fill = {{ y.var }}),
             color = "gray") +
    geom_text(data = totals,
              aes(x = {{ x.var }}, y = 1, label = total_n), vjust = -0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  if (!is.null(color.scale)) p <- p + scale_fill_manual(values = color.scale)
  if (return.table) return(list(p, freq)) else return(p)
}


# =============================================================================
# SECTION 0: Initial Level1 annotation on integrated object
# =============================================================================
# Assigns broad cell type labels from RPCA cluster IDs, adds composite
# metadata columns, and generates diagnostic plots to confirm annotation
# before running inferCNV.

cat("Loading integrated object...\n")
merged.all <- readRDS(path.integrated.rds)
merged.all <- RunUMAP(merged.all, dims = 1:40)

# -- First-pass annotation from RPCA clusters (res 1.4) ----------------------
merged.all@meta.data <- mutate(
  merged.all@meta.data,
  Level1.anno = case_when(
    rpca_clusters_1.4 %in% as.character(LEVEL1_CLUSTERS$Myeloid)     ~ "Myeloid",
    rpca_clusters_1.4 %in% as.character(LEVEL1_CLUSTERS$T_NK)        ~ "T_NK",
    rpca_clusters_1.4 %in% as.character(LEVEL1_CLUSTERS$B)           ~ "B",
    rpca_clusters_1.4 %in% as.character(LEVEL1_CLUSTERS$Ambiguous_1) ~ "Ambiguous_1",
    rpca_clusters_1.4 %in% as.character(LEVEL1_CLUSTERS$Ambiguous_2) ~ "Ambiguous_2",
    .default = "Blast"
  )
)

# Inspect ambiguous clusters: review VlnPlot + DotPlot to decide whether to
# relabel as Blast or a normal cell type, then update LEVEL1_FINAL_CLUSTERS
pdf(file.path(dir.plots, "InitialAnnotation_AmbiguousClusterReview.pdf"),
    height = 8, width = 7)

DimPlot(merged.all, reduction = "umap.rpca.pc40",
        group.by = "rpca_clusters_1.4", label = TRUE, repel = TRUE) +
  NoLegend() + ggtitle("RPCA clusters (res 1.4)")

VlnPlot(merged.all,
        features = c("CD19", "MS4A1", "CD22", "CD24", "CD79A", "MME",
                     "CD34", "CD38", "PTPRC", "nCount_RNA", "nFeature_RNA",
                     "max.score.celltype.short"),
        group.by = "Level1.anno", stack = TRUE, flip = TRUE) + NoLegend()

DimPlot(merged.all, group.by = "Level1.anno",
        reduction = "umap.rpca.pc40") +
  DimPlot(merged.all, group.by = "celltype.short",
          reduction = "umap.rpca.pc40")

DotPlot(merged.all,
        features = c("CD19", "MME", "PTPRC", "CD34", "MS4A1",
                     "CEBPA", "CD14", "ITGAM", "CD3E"),
        group.by = "RNA_snn_res.1", scale = TRUE,
        col.min = -1, col.max = 1, cluster.idents = TRUE) +
  RotatedAxis() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

dev.off()

# -- Final annotation using RNA_snn_res.1 after ambiguity review --------------
merged.all@meta.data <- mutate(
  merged.all@meta.data,
  Level1.anno = case_when(
    RNA_snn_res.1 %in% as.character(LEVEL1_FINAL_CLUSTERS$Myeloid) ~ "Myeloid",
    RNA_snn_res.1 %in% as.character(LEVEL1_FINAL_CLUSTERS$T_NK)    ~ "T_NK",
    RNA_snn_res.1 %in% as.character(LEVEL1_FINAL_CLUSTERS$B)       ~ "B",
    .default = "Blast"
  ),
  Level1.anno.Sample = case_when(
    Level1.anno == "Blast" ~ paste0("Blast_", Sample),
    .default = as.character(Level1.anno)
  )
)

# Factor levels for consistent ordering in plots
merged.all$Level1.anno.Sample <- factor(
  merged.all$Level1.anno.Sample,
  levels = c("Blast_4526R", "Blast_5481",  "Blast_5995",  "Blast_6024",
             "Blast_6092",  "Blast_6061",  "Blast_6099R", "Blast_6146",
             "Blast_5850",  "Blast_5907",  "B", "Myeloid", "T_NK")
)

merged.all$celltype.short <- factor(
  merged.all$celltype.short,
  levels = c(B_CTYPES, HSPC_CTYPES,
             "GMP", "MEP", "DC-Progenitor", "pDC", "cDC", "Mono", T_CTYPES)
)

merged.all@meta.data <- mutate(
  merged.all@meta.data,
  patient.sample.Stage = case_when(
    Level1.anno == "Blast" ~ paste0(patient.sample, "_", Stage),
    .default = Level1.anno
  )
)

merged.all$patient.sample.Stage <- factor(
  merged.all$patient.sample.Stage,
  levels = c("T_NK", "B", "Myeloid",
             "PT1_IDX", "PT1_MRD", "PT2_IDX", "PT2_MRD",
             "PT3_IDX", "PT3_MRD", "PT4_IDX", "PT4_MRD",
             "PT5_IDX", "PT5_MRD")
)

# Assign blast lineage using projection labels
merged.all <- assign_blast_lineage(merged.all)

# Highlight individual samples to verify annotation
pdf(file.path(dir.plots, "InitialAnnotation_SampleHighlights.pdf"),
    height = 6, width = 8)
DimPlot(merged.all, group.by = "Level1.anno", reduction = "umap.rpca.pc40") +
  DimPlot(merged.all, group.by = "celltype.short",  reduction = "umap.rpca.pc40")
DimPlot(merged.all, group.by = "Stage")
for (sname in unique(merged.all$Sample)) {
  print(
    DimPlot(merged.all,
            cells.highlight = colnames(subset(merged.all,
                                               subset = Sample == sname)),
            sizes.highlight = 0.5) + ggtitle(sname)
  )
}
dev.off()

saveRDS(merged.all,
        gsub("\\.rds", "_preInferCNV.rds", path.output.rds))
cat("Pre-inferCNV annotated object saved.\n")


# =============================================================================
# SECTION 1: ADT heatmaps and RNA-ADT correlation (pre-inferCNV)
# =============================================================================
# Generates ADT heatmaps grouped by patient.sample.Stage and blast celltype,
# plus per-marker RNA–ADT scatter plots and a per-sample correlation matrix.
# These are run on the initial annotation before inferCNV refines blast labels.

# -- ADT heatmaps by patient and blast celltype --------------------------------
blasts.all <- subset(merged.all, subset = Level1.anno == "Blast")
blasts.all$celltype.short <- factor(
  blasts.all$celltype.short,
  levels = c(HSPC_CTYPES, B_CTYPES, "MEP", "GMP", "Mono",
             "DC-Progenitor", "pDC", "cDC")
)
blasts.subset <- subset(blasts.all,
                         subset = celltype.short %in%
                           c(HSPC_CTYPES, B_CTYPES,
                             "MEP", "GMP", "Mono", "DC-Progenitor",
                             "pDC", "cDC"))

pdf(file.path(dir.plots, "ADTHeatmaps_PreInferCNV.pdf"),
    height = 8, width = 5)
adt_heatmap(merged.all, groupByVar = "patient.sample.Stage")
adt_heatmap(merged.all, groupByVar = "patient.sample.Stage",
            clusterIdents = FALSE)
adt_heatmap(merged.all, groupByVar = "patient.sample.Stage",
            prots = c("CD305", "CD204", "CD105", "CD98", "CD119", "CD83.1",
                       "CD371", "CD9.1", "CD298", "CD13",  "CD33.1",
                       "CD19.1", "CD10"))
adt_heatmap(blasts.subset, groupByVar = "celltype.short")
adt_heatmap(blasts.subset, groupByVar = "celltype.short", clusterIdents = FALSE)
dev.off()

# Highlight celltype-like blast populations on UMAP
celltypes.highlight <- c(HSPC_CTYPES, B_CTYPES,
                          "MEP", "GMP", "Mono", "DC-Progenitor", "pDC", "cDC")
pdf(file.path(dir.plots, "ADTHeatmaps_BlastHighlights.pdf"),
    height = 3, width = 5)
for (ct in celltypes.highlight) {
  Idents(blasts.all) <- "celltype.short"
  print(
    DimPlot(merged.all,
            cells.highlight = WhichCells(blasts.all, idents = ct),
            sizes.highlight = 0.5) +
      ggtitle(paste0(ct, "-like Blasts"))
  )
}
dev.off()

# -- RNA-ADT scatter and correlation plots ------------------------------------
# Assign ADT concentration batch (affects CLR normalization comparability)
merged.all@meta.data <- mutate(
  merged.all@meta.data,
  Concentration_Batch = case_when(
    Sample %in% ADT_BATCH_MAP$ADT_1 ~ "ADT_1",
    Sample %in% ADT_BATCH_MAP$ADT_2 ~ "ADT_2",
    .default = "ADT_3"
  )
)

merged.all <- JoinLayers(merged.all, assay = "RNA")
merged.all <- JoinLayers(merged.all, assay = "ADT")

# Scatter plots (downsampled for speed)
merged.ds <- subset(merged.all, downsample = 500)
scatter.list <- lapply(seq_len(nrow(ADT_RNA_PAIRS)), function(j) {
  FeatureScatter(merged.ds,
                 feature1   = ADT_RNA_PAIRS$RNA[j],
                 feature2   = ADT_RNA_PAIRS$ADT[j],
                 pt.size    = 0.5,
                 group.by   = "Concentration_Batch") +
    labs(title = ADT_RNA_PAIRS$ADT[j]) +
    theme(plot.title  = element_text(size = 12),
          axis.text   = element_text(size = 7),
          axis.title  = element_text(size = 7))
})

pdf(file.path(dir.plots, "ADT_RNA_CorrelationScatters.pdf"),
    height = 12, width = 10)
print(ggarrange(plotlist = scatter.list, ncol = 6, nrow = 7,
                common.legend = TRUE, legend = "bottom"))
dev.off()

# Per-sample Pearson correlation matrix (RNA vs ADT, CLR normalized)
merged.all <- NormalizeData(merged.all, assay = "ADT",
                             normalization.method = "CLR")
sample.objs    <- SplitObject(merged.all, split.by = "Sample")
pairs.no.cd90  <- ADT_RNA_PAIRS[ADT_RNA_PAIRS$ADT != "CD90", ]  # CD90 excluded

cor.mtx <- sapply(names(sample.objs), function(sname) {
  obj     <- sample.objs[[sname]]
  rna.mtx <- GetAssayData(obj, layer = "data", assay = "RNA")
  adt.mtx <- GetAssayData(obj, layer = "data", assay = "ADT")
  sapply(seq_len(nrow(pairs.no.cd90)), function(j) {
    cor(rna.mtx[pairs.no.cd90$RNA[j], ],
        adt.mtx[pairs.no.cd90$ADT[j], ])
  })
})
rownames(cor.mtx) <- pairs.no.cd90$ADT

prot_col <- circlize::colorRamp2(
  breaks = c(min(cor.mtx, na.rm = TRUE), 0, max(cor.mtx, na.rm = TRUE)),
  colors = c("blue", "white", "red")
)

pdf(file.path(dir.plots, "ADT_RNA_CorrelationMatrix_BySample.pdf"),
    height = 7, width = 6)
print(Heatmap(cor.mtx,
              name             = "CLR Normalization",
              col              = prot_col,
              use_raster       = TRUE,
              na_col           = "grey",
              cluster_rows     = FALSE,
              cluster_columns  = FALSE,
              row_names_gp     = gpar(color = "black", fontsize = 10),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.2f", cor.mtx[i, j]), x, y,
                          gp = gpar(fontsize = 6, col = "black"))
              }))
dev.off()


# =============================================================================
# SECTION 2: Per-pair inferCNV re-clustering and Level1 re-annotation
# =============================================================================
# For each patient pair, merges IDX + MRD inferCNV objects, re-clusters,
# and refines Level1.anno using CNV scores and marker genes.
# All logic is driven by pair.config — no per-pair code duplication.

infercnv.objects <- lapply(infercnv.paths, readRDS)

pdf(file.path(dir.plots, "inferCNVResults_PairwiseClustering.pdf"),
    height = 7, width = 8)

pair.objects <- list()

for (pair_name in names(pair.config)) {
  cfg <- pair.config[[pair_name]]
  cat("\nProcessing", pair_name, ":", cfg$idx.sample, "vs", cfg$mrd.sample, "\n")

  pair.s <- merge(infercnv.objects[[cfg$idx.sample]],
                  infercnv.objects[[cfg$mrd.sample]])
  pair.s <- pair.s %>%
    ScaleData() %>% FindVariableFeatures() %>%
    RunPCA(ndims = 40) %>% FindNeighbors() %>%
    FindClusters(resolution = 1.0, cluster.name = "merged_res_1.0") %>%
    RunUMAP(dims = 1:40)

  print(DimPlot(pair.s, group.by = "Sample") +
        DimPlot(pair.s, group.by = "celltype.binned"))
  print(FeaturePlot(pair.s, features = cfg$cnv.features,
                    ncol = 2, order = TRUE) & NoAxes())
  print(VlnPlot(pair.s, features = cfg$cnv.features,
                stack = TRUE, flip = TRUE,
                group.by = "merged_res_1.0") + NoLegend())

  pair.s@meta.data <- mutate(
    pair.s@meta.data,
    Level1.anno = case_when(
      merged_res_1.0 %in% as.character(cfg$clusters.myeloid) ~ "Myeloid",
      merged_res_1.0 %in% as.character(cfg$clusters.tnk)     ~ "T_NK",
      merged_res_1.0 %in% as.character(cfg$clusters.b)       ~ "B",
      .default = "Blast"
    )
  )

  print(DimPlot(pair.s, group.by = "Sample") +
        (DimPlot(pair.s, group.by = "Level1.anno") + ggtitle("Level1 Anno")))
  print(FeaturePlot(pair.s, features = c("CD19", "CD33", "MS4A1", "CD3E"),
                    order = TRUE) & NoAxes())

  pair.s <- assign_blast_lineage(pair.s)
  print(DimPlot(pair.s, group.by = "Stage") +
        DimPlot(pair.s, group.by = "Blast_Lineage", cols = LINEAGE_COLS))

  pair.objects[[pair_name]] <- pair.s
}
dev.off()


# =============================================================================
# SECTION 3: Transfer CNV metadata to merged integrated object
# =============================================================================

for (pair.s in pair.objects) {
  meta.transfer <- pair.s@meta.data %>% dplyr::select(all_of(METADATA_CNV_COLS))
  merged.all    <- AddMetaData(merged.all, meta.transfer)
}

# Clean up doublet-finder columns from preprocessing
merged.all@meta.data <- merged.all@meta.data[,
  !grepl("^DF|^pANN", colnames(merged.all@meta.data))]
merged.all$Level1.anno.Stage    <- NULL
merged.all$patient.sample.Stage <- NULL

# Rebuild blast lineage and composite grouping variables on updated object
merged.all <- assign_blast_lineage(merged.all)

merged.all@meta.data <- mutate(
  merged.all@meta.data,
  Level1.anno.Sample = case_when(
    Level1.anno == "Blast" ~ paste0("Blast_", Sample),
    .default = as.character(Level1.anno)
  ),
  Blast_Stage = case_when(
    Level1.anno == "Blast" & Stage == "IDX" ~ "Blast_IDX",
    Level1.anno == "Blast" & Stage == "MRD" ~ "Blast_MRD",
    .default = Level1.anno
  ),
  Subtype = case_when(
    patient.sample %in% SUBTYPE_MAP$KMT2Ar ~ "KMT2Ar",
    .default = "Ph-like"
  ),
  Patient.Stage.Healthy = case_when(
    Level1.anno == "Blast" ~ Patient.Stage.All,
    .default = Level1.anno
  ),
  Blast_Lineage_Stage = paste0(Blast_Lineage, "_", Stage)
)

# Set factor levels
merged.all$Blast_Stage <- factor(
  merged.all$Blast_Stage,
  levels = c("Blast_IDX", "Blast_MRD", "B", "T_NK", "Myeloid"))

merged.all$Blast_Lineage_Stage <- factor(
  merged.all$Blast_Lineage_Stage,
  levels = c(paste0(BLAST_LINEAGE_LEVELS[1:4], "_IDX"),
             paste0(BLAST_LINEAGE_LEVELS[1:4], "_MRD"),
             paste0(c("B", "T_NK", "Myeloid"), "_IDX"),
             paste0(c("B", "T_NK", "Myeloid"), "_MRD")))

merged.all$Patient.Stage.Healthy <- factor(
  merged.all$Patient.Stage.Healthy,
  levels = c("B", "Myeloid", "T_NK",
             "PT1_IDX", "PT1_MRD", "PT2_IDX", "PT2_MRD",
             "PT3_IDX", "PT3_MRD", "PT4_IDX", "PT4_MRD",
             "PT5_IDX", "PT5_MRD"))

merged.all <- RunUMAP(merged.all, dims = 1:40)

saveRDS(merged.all, path.output.rds)
cat("Annotated object with CNV metadata saved to:", path.output.rds, "\n")

blasts.all <- subset(merged.all, subset = Level1.anno == "Blast")

blast.freq.colors <- c(
  "#619ccf","#75a8d5","#89b5db","#9cc1e1","#b0cee7","#c4daed",
  "#00ba38","#60d483",
  "#f8766d","#f9877f","#fa9892","#fba9a4","#fbbbb6","#fcccc8",
  "#703ac5","#8050cc","#9066d2","#a07cd9","#b092df","#bfa7e5",
  "#cfbdec","#dfd3f2","#efe9f9")

pdf(file.path(dir.plots, "FrequencyPlots_PostInferCNVAnno.pdf"),
    height = 6, width = 9)
cluster.plot(merged.all@meta.data, Level1.anno.Sample, celltype.short,
             color.scale = blast.freq.colors) + ggtitle("All Cells")
cluster.plot(blasts.all@meta.data, Level1.anno.Sample, celltype.short,
             color.scale = blast.freq.colors) + ggtitle("Blasts")
dev.off()


# =============================================================================
# SECTION 4: ADT heatmaps post-inferCNV
# =============================================================================

pdf(file.path(dir.plots, "ADTHeatmaps_PostInferCNV_ByLineageAndStage.pdf"),
    height = 6, width = 6)
adt_heatmap(merged.all, prots = CORE_ADT_PANEL,
            groupByVar = "Blast_Lineage",      clusterIdents = FALSE,
            title = "Median ADT — Blast Lineage")
adt_heatmap(merged.all, prots = CORE_ADT_PANEL,
            groupByVar = "Blast_Stage",         clusterIdents = FALSE,
            title = "Median ADT — Blast Stage")
adt_heatmap(merged.all, prots = CORE_ADT_PANEL,
            groupByVar = "Blast_Lineage_Stage", clusterIdents = FALSE,
            title = "Median ADT — Blast Lineage × Stage")
adt_heatmap(subset(merged.all, subset = Stage == "IDX"),
            prots = CORE_ADT_PANEL, groupByVar = "Blast_Lineage_Stage",
            clusterIdents = FALSE, title = "Median ADT — IDX only")
adt_heatmap(subset(merged.all, subset = Stage == "MRD"),
            prots = CORE_ADT_PANEL, groupByVar = "Blast_Lineage_Stage",
            clusterIdents = FALSE, title = "Median ADT — MRD only")
dev.off()

merged.samples <- SplitObject(merged.all, split.by = "Sample")
pdf(file.path(dir.plots, "ADTHeatmaps_PostInferCNV_BySample.pdf"),
    height = 8, width = 5)
for (sname in names(merged.samples)) {
  adt_heatmap(merged.samples[[sname]], groupByVar = "Level1.anno.Sample",
              clusterIdents = FALSE,
              title = paste0("Median ADT — ", sname, " Level1 anno"))
  adt_heatmap(merged.samples[[sname]], groupByVar = "Blast_Lineage",
              clusterIdents = FALSE,
              title = paste0("Median ADT — ", sname, " Lineage"))
}
dev.off()


# =============================================================================
# SECTION 5: Supplementary UMAP plots
# =============================================================================

patient.stage.colors <- c(
  "#944939", "#7F7A3A", "#1B4A4A",
  "#f25e60", "#e41a1c", "#6ea8d6", "#377eb8",
  "#7fd977", "#4daf4a", "#c285cd", "#984ea3",
  "#ffb54d", "#ff7f00"
)

cnv.highlight.params <- list(
  list(sample = "6099R", cnv_col = "has_dupli_chr21", label = "chr21 Duplication"),
  list(sample = "6099R", cnv_col = "has_dupli_chrX",  label = "chrX Duplication"),
  list(sample = "6146",  cnv_col = "has_dupli_chr21", label = "chr21 Duplication"),
  list(sample = "6146",  cnv_col = "has_dupli_chrX",  label = "chrX Duplication")
)

pdf(file.path(dir.plots, "SupplementaryUMAPs.pdf"), height = 5, width = 7)
DimPlot(merged.all, group.by = "Stage")
DimPlot(merged.all, group.by = "patient.sample")
DimPlot(merged.all, group.by = "Subtype")
DimPlot(merged.all, group.by = "Patient.Stage.All")
DimPlot(merged.all, group.by = "Patient.Stage.Healthy",
        cols = patient.stage.colors)
for (params in cnv.highlight.params) {
  print(
    DimPlot(merged.all, reduction = "umap.rpca.pc40",
            cells.highlight = colnames(subset(
              merged.all,
              subset = Sample == params$sample &
                       .data[[params$cnv_col]] == TRUE)),
            sizes.highlight = 0.5) +
      ggtitle(paste0(params$sample, " — ", params$label))
  )
}
dev.off()


# =============================================================================
# SECTION 6: ADT differential expression and volcano plots
# =============================================================================

merged.all <- JoinLayers(merged.all, assay = "ADT")

de.comparisons <- list(
  list(ident1 = "Myeloid-like_IDX", ident2 = "B-like_IDX",       title = "Myeloid vs B (IDX)"),
  list(ident1 = "Myeloid-like_MRD", ident2 = "B-like_MRD",       title = "Myeloid vs B (MRD)"),
  list(ident1 = "B-like_MRD",       ident2 = "B-like_IDX",       title = "B-like: MRD vs IDX"),
  list(ident1 = "Myeloid-like_MRD", ident2 = "Myeloid-like_IDX", title = "Myeloid-like: MRD vs IDX")
)

de.results <- lapply(de.comparisons, function(comp) {
  FindMarkers(merged.all, ident.1 = comp$ident1, ident.2 = comp$ident2,
              assay = "ADT", group.by = "Blast_Lineage_Stage",
              logfc.threshold = 0)
})

make_volcano <- function(markers, title) {
  col_vec <- ifelse(
    markers$p_val_adj < DE_P_CUT & markers$avg_log2FC >  DE_FC_CUT, "red",
    ifelse(markers$p_val_adj < DE_P_CUT & markers$avg_log2FC < -DE_FC_CUT,
           "blue", "grey"))
  names(col_vec) <- ifelse(col_vec == "red", "Up",
                    ifelse(col_vec == "blue", "Down", "NS"))
  EnhancedVolcano(markers,
    lab = rownames(markers), x = "avg_log2FC", y = "p_val_adj",
    title = title, xlab = "log2 fold change",
    ylab = "-log10 adjusted p-value",
    pCutoff = DE_P_CUT, FCcutoff = DE_FC_CUT,
    colCustom = col_vec, colAlpha = 0.4,
    legendLabels = c("NS", "NS", "Up", "Down"),
    gridlines.major = FALSE, gridlines.minor = FALSE,
    border = "full", borderWidth = 1.0, borderColour = "black")
}

volcano.plots <- mapply(make_volcano, de.results,
                         lapply(de.comparisons, `[[`, "title"),
                         SIMPLIFY = FALSE)

pdf(file.path(dir.plots, "ADT_MarkerVolcanoPlots.pdf"),
    height = 12, width = 10)
print(ggarrange(plotlist = volcano.plots, ncol = 2, nrow = 2))
dev.off()

merged.blasts <- subset(merged.all, subset = Level1.anno == "Blast")
ggsave(
  RidgePlot(merged.blasts, features = "CD305", assay = "ADT",
             group.by = "Blast_Lineage_Stage"),
  filename = file.path(dir.plots, "CD305_RidgePlot.pdf")
)

cat("CITE-Seq annotation and visualization complete.\n")
