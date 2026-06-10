# =============================================================================
# Projection Annotation Update and Blast Cell Abundance Analysis
# =============================================================================
#
# Description:
#   Updates cell type annotations on projected blast cells using a revised
#   healthy donor reference annotation, then performs a comprehensive
#   abundance analysis examining how cell type frequencies change between
#   diagnosis (IDX) and MRD timepoints.
#
#   The script is organized into four sections:
#     1. Annotation update  — joins new healthy reference labels onto existing
#                              projection results and adds them to the global
#                              Seurat object
#     2. Cell frequency computation — calculates per-sample blast cell type
#                              frequencies (raw and log2 fold-change IDX→MRD)
#     3. Comparison plots   — paired plots stratified by Subtype, MRD Group,
#                              and trajectory cluster (k-means); Wilcoxon
#                              paired tests and pairwise comparisons
#     4. Cell-type-level pairwise comparisons — per-cell-type boxplots with
#                              all pairwise group comparisons
#
# Inputs (set in USER CONFIGURATION):
#   - Updated healthy donor metadata CSV (new cell type labels)
#   - Global merged Seurat RDS (post-labeling; all subtypes combined)
#   - Previous projection results CSV (from 01_project_onto_healthy_reference.R)
#   - Patient number lookup CSV  (cols: sample_patient, PT.number)
#   - Trajectory k-means cluster CSV (cols: Row.names, k_4_Myeloid)
#
# Outputs:
#   - Updated Seurat RDS and projection metadata CSV
#   - PDF figures: paired abundance plots by Subtype / MRD Group / Traj Cluster
#   - EPS figures: aggregate and per-cell-type paired plots
#   - CSV: cell type frequency table and LFC table
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Input files
path.hd.meta         <- "/path/to/HealthyRef_Reannotation_metadata.csv"
                        # cols: X (barcode), NewCellTypeUpdate_r3
path.global.rds      <- "/path/to/RObjects/merged.postLabeling.rds"
path.projection.csv  <- "/path/to/metadata/Projection_MultiomeHD/AllSamples_proj2healthyCells.csv"
path.pt.number.csv   <- "/path/to/PatientData.csv"
                        # cols: sample_patient, PT.number
path.kcluster.csv    <- "/path/to/Tables/Ctype_aggregated_scaled_LFC_k4clustering.csv"
                        # cols: Row.names (PT.number), k_4_Myeloid

# Output files
dir.plots            <- "/path/to/Plots/Updated_Annotation/"
path.freq.csv        <- "/path/to/Tables/CtypeAggFrequency_UpdatedAnno.csv"
path.projection.updated.csv <- "/path/to/Tables/UpdatedProjectionResults_AllCells.csv"

# Cell type ordering — edit if your reference uses different labels
CTYPE_ORDER <- c("HSC/MPP", "LMPP", "CLP", "Pre-pro-B", "Pro-B", "Pre-B",
                 "Immature-B", "Mature-B-1", "Mature-B-2",
                 "MEP", "GMP", "DC-progenitor", "pDC", "cDC", "Monocyte",
                 "T-cell", "NK")

STAGE_ORDER <- c("IDX", "MRD", "RLP")

# Aggregated cell type groupings — maps fine-grained labels to broader classes
# used for abundance analysis. Edit to match your reference annotation.
CTYPE_AGGREGATE_MAP <- list(
  "HSPC"             = c("HSC/MPP", "LMPP"),
  "Early-Lymphoid"   = c("CLP", "Pre-pro-B"),
  "Pro-B"            = c("Pro-B"),
  "Pre-B"            = c("Pre-B"),
  "Immature-Mature-B"= c("Immature-B", "Mature-B-1", "Mature-B-2"),
  "Myeloid-Progenitor"= c("MEP", "GMP", "DC-progenitor"),
  "Mature-Myeloid"   = c("pDC", "cDC", "Monocyte"),
  "T-NK"             = c("T-cell", "NK")
)

# Ordered aggregate cell type levels for plotting (exclude T-NK from
# abundance comparisons as these are normal bystander cells in blast gate)
AGGREGATE_ORDER <- c("Mature-Myeloid", "Myeloid-Progenitor", "HSPC",
                     "Early-Lymphoid", "Pro-B", "Pre-B", "Immature-Mature-B")

# Patients to exclude from abundance analysis (e.g. insufficient paired data)
EXCLUDE_PATIENTS <- c("PT17", "PT5")

# Significance annotation style for stat_compare_means
SYMNUM_ARGS <- list(
  cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
  symbols   = c("****", "***", "**", "*", "ns")
)


# =============================================================================
# SECTION 1: Update projection annotations
# =============================================================================
# The healthy reference was re-annotated after the initial projection run.
# This section joins the new labels (NewCellTypeUpdate_r3) onto the existing
# projection results via the healthy cell barcode, then adds all projection
# metadata to the global Seurat object.

cat("Loading data...\n")
hd.meta         <- read.csv(path.hd.meta)
global.rna      <- readRDS(path.global.rds)
all.project.data <- read.delim(path.projection.csv, row.names = "patient_bc")
pt.number       <- read.csv(path.pt.number.csv)
k.cluster       <- read.csv(path.kcluster.csv) %>%
  select(Row.names, k_4_Myeloid) %>%
  rename(PT.number = Row.names, k_cluster = k_4_Myeloid)

# Build new annotation lookup from updated healthy metadata
new.anno <- hd.meta %>%
  select(X, NewCellTypeUpdate_r3) %>%
  rename(healthy_bc = X, Projected.Ctype.Wmethod = NewCellTypeUpdate_r3)

# Preserve old annotation for QC comparison, then replace with updated labels
all.project.data$old.projection.anno   <- all.project.data$Projected.Ctype.Wmethod
all.project.data$Projected.Ctype.Wmethod <- NULL

# Join patient numbers and trajectory clusters
all.project.data <- all.project.data %>%
  left_join(pt.number, by = c("Patient.Sample" = "sample_patient")) %>%
  left_join(k.cluster, by = "PT.number")

# Join updated cell type labels
new.project.data            <- left_join(all.project.data, new.anno,
                                          by = "healthy_bc")
rownames(new.project.data)  <- new.project.data$patient_bc

# Verify annotation update — inspect the cross-tabulation to confirm
# the new labels are sensible relative to the old ones
cat("Annotation update cross-tabulation (old vs new):\n")
print(table(new.project.data$old.projection.anno,
            new.project.data$Projected.Ctype.Wmethod))

# Add updated projection metadata to global Seurat object
global.rna <- AddMetaData(global.rna, new.project.data)

# Set factor levels for downstream plotting
global.rna$Projected.Ctype.Wmethod <- factor(global.rna$Projected.Ctype.Wmethod,
                                              levels = CTYPE_ORDER)
global.rna$Stage <- factor(global.rna$Stage, levels = STAGE_ORDER)

# Add aggregated cell type column using the mapping defined in USER CONFIGURATION
# This replaces the original for-loop with a vectorized case_when
build_aggregate <- function(ctype_vec, agg_map) {
  result <- rep(NA_character_, length(ctype_vec))
  for (agg_label in names(agg_map)) {
    result[ctype_vec %in% agg_map[[agg_label]]] <- agg_label
  }
  result
}

global.rna@meta.data$Ctype.Aggregate <- build_aggregate(
  as.character(global.rna$Projected.Ctype.Wmethod),
  CTYPE_AGGREGATE_MAP
)

global.rna$Ctype.Aggregate <- factor(global.rna$Ctype.Aggregate,
                                      levels = c(AGGREGATE_ORDER, "T-NK"))

# Save updated projection metadata
write.csv(global.rna@meta.data, file = path.projection.updated.csv)
cat("Updated projection metadata saved.\n")


# =============================================================================
# SECTION 2: Compute blast cell type frequencies
# =============================================================================
# For each patient sample, computes the fraction of blast cells projected to
# each aggregated cell type at IDX and MRD. Filters to:
#   - Blast cells only (Level1.anno == "Blast")
#   - IDX and MRD timepoints (RLP excluded from paired analysis)
#   - T-NK excluded (normal bystander cells that contaminate blast gate)
#   - Patients excluded per EXCLUDE_PATIENTS
#   - Only patients with both IDX and MRD present (paired = 2)

cat("Computing cell type frequencies...\n")

Ctype_count_byStage <- global.rna@meta.data %>%
  dplyr::filter(
    Level1.anno        == "Blast",
    Stage              %in% c("IDX", "MRD"),
    Ctype.Aggregate    %notin% c("T-NK"),
    PT.number          %notin% EXCLUDE_PATIENTS
  ) %>%
  select(Patient.Sample, Stage, PT.number, Ctype.Aggregate,
         Projected.Ctype.Wmethod, Subtype, k_cluster, MRD.Group) %>%
  group_by(Patient.Sample, Stage) %>%
  mutate(total_n = n()) %>%
  group_by(Patient.Sample, Stage, Ctype.Aggregate) %>%
  mutate(cell_n = n(),
         freq   = cell_n / total_n) %>%
  ungroup() %>%
  select(-Projected.Ctype.Wmethod) %>%
  distinct() %>%
  # Keep only patients with both IDX and MRD for a given cell type (paired)
  group_by(PT.number, Ctype.Aggregate) %>%
  filter(n() == 2) %>%
  ungroup() %>%
  arrange(Ctype.Aggregate, PT.number)

# Set factor levels
Ctype_count_byStage$Ctype.Aggregate <- factor(Ctype_count_byStage$Ctype.Aggregate,
                                               levels = AGGREGATE_ORDER)
Ctype_count_byStage$Stage <- factor(Ctype_count_byStage$Stage,
                                     levels = c("IDX", "MRD"))

write.csv(Ctype_count_byStage, file = path.freq.csv)

# Split into IDX and MRD for group comparison plots
idx.ctype <- Ctype_count_byStage %>% filter(Stage == "IDX")
mrd.ctype <- Ctype_count_byStage %>% filter(Stage == "MRD")

# -- Log2 fold-change (MRD / IDX) per patient per cell type ------------------
# A small pseudocount (the minimum non-zero frequency in the data) is added
# before log-transformation to avoid log(0). This is more principled than
# an arbitrary fixed value.
pseudocount <- min(Ctype_count_byStage$freq[Ctype_count_byStage$freq > 0])

Ctype_count_lfc <- Ctype_count_byStage %>%
  select(Ctype.Aggregate, Stage, Subtype, PT.number, freq, MRD.Group, k_cluster) %>%
  distinct() %>%
  group_by(Ctype.Aggregate, PT.number) %>%
  pivot_wider(names_from = Stage, values_from = freq) %>%
  replace(is.na(.), 0) %>%
  mutate(lfc = log2((MRD + pseudocount) / (IDX + pseudocount)))

Ctype_count_lfc$Ctype.Aggregate <- factor(Ctype_count_lfc$Ctype.Aggregate,
                                           levels = AGGREGATE_ORDER)

# Set factor levels for grouping variables
for (df_name in c("idx.ctype", "mrd.ctype", "Ctype_count_lfc")) {
  df <- get(df_name)
  df$MRD.Group <- factor(df$MRD.Group, levels = c("<1", ">1<5", ">5"))
  df$k_cluster <- factor(df$k_cluster, levels = c(1, 2, 3, 4))
  assign(df_name, df)
}


# =============================================================================
# SECTION 3: Paired abundance comparison plots
# =============================================================================
# Generates paired plots (ggpaired) for three grouping variables:
#   - Subtype (KMT2Ar, iAMP21, Ph-like, Ph-positive)
#   - MRD Group (<1%, >1<5%, >5%)
#   - Trajectory cluster (k-means clusters 1-4)
#
# Each PDF contains: IDX frequencies, MRD frequencies, and IDX→MRD LFC,
# each faceted by aggregated cell type.

# Helper: build a standard paired + LFC plot set for one grouping variable
plot_abundance_by_group <- function(idx.df, mrd.df, lfc.df,
                                    group.var, group.label,
                                    lfc.label.y = -6) {
  list(
    ggpaired(idx.df, x = group.var, y = "freq",
             color = group.var, line.color = "gray",
             facet.by = "Ctype.Aggregate", add = "jitter",
             add.params = list(size = 1),
             title = paste("IDX Cell Type Frequencies by", group.label)) +
      stat_compare_means(label.y = 0.8),

    ggpaired(mrd.df, x = group.var, y = "freq",
             color = group.var, line.color = "gray",
             facet.by = "Ctype.Aggregate", add = "jitter",
             add.params = list(size = 1),
             title = paste("MRD Cell Type Frequencies by", group.label)) +
      stat_compare_means(label.y = 0.8),

    ggpaired(lfc.df, x = group.var, y = "lfc",
             color = group.var, line.color = "gray",
             facet.by = "Ctype.Aggregate", add = "jitter",
             add.params = list(size = 1),
             title = paste("LFC (MRD/IDX) by", group.label)) +
      stat_compare_means(label.y = lfc.label.y)
  )
}

# -- Subtype comparisons ------------------------------------------------------
cat("Plotting by Subtype...\n")
pdf(file.path(dir.plots, "CType_abundance_bySubtype.pdf"), height = 8, width = 10)
for (p in plot_abundance_by_group(idx.ctype, mrd.ctype, Ctype_count_lfc,
                                   "Subtype", "Subtype", lfc.label.y = -7)) {
  print(p)
}
dev.off()

# -- MRD Group comparisons ----------------------------------------------------
cat("Plotting by MRD Group...\n")
pdf(file.path(dir.plots, "CType_abundance_byMRDGroup.pdf"), height = 8, width = 10)
for (p in plot_abundance_by_group(idx.ctype, mrd.ctype, Ctype_count_lfc,
                                   "MRD.Group", "MRD Group", lfc.label.y = -6)) {
  print(p)
}
dev.off()

# -- Trajectory cluster comparisons -------------------------------------------
cat("Plotting by Trajectory Cluster...\n")
pdf(file.path(dir.plots, "CType_abundance_byTrajCluster.pdf"), height = 4, width = 12)
for (p in plot_abundance_by_group(idx.ctype, mrd.ctype, Ctype_count_lfc,
                                   "k_cluster", "Trajectory Cluster",
                                   lfc.label.y = -6)) {
  print(p + facet_wrap(~ Ctype.Aggregate, ncol = 8))
}
dev.off()

# -- All-patient paired plot (aggregate, IDX vs MRD) --------------------------
# Primary summary figure: cell type frequency changes across all patients
cat("Plotting all-patient paired abundance...\n")

p.paired.all <- ggpaired(Ctype_count_byStage, x = "Stage", y = "freq",
                          id = "PT.number", color = "Stage",
                          line.color = "gray", line.size = 0.05,
                          facet.by = "Ctype.Aggregate", nrow = 1,
                          point.size = 0.5, ggtheme = theme_classic(),
                          legend = "none") +
  stat_compare_means(paired = TRUE, label = "p.format",
                     label.x.npc = "left", hide.ns = TRUE, label.y = 0.75) +
  theme(strip.text = element_text(size = 8))

ggsave(p.paired.all, device = "eps", width = 7, height = 3, units = "in",
       filename = file.path(dir.plots, "CellFrequencyPaired_AllPatients.eps"))

# -- Aggregate frequency bar plot by stage ------------------------------------
# Shows overall cell type composition at IDX, MRD, and RLP

blast.freq.stage <- global.rna@meta.data %>%
  dplyr::filter(Level1.anno == "Blast") %>%
  select(Stage, Projected.Ctype.Wmethod) %>%
  mutate(Ctype.Collapsed = build_aggregate(as.character(Projected.Ctype.Wmethod),
                                            CTYPE_AGGREGATE_MAP)) %>%
  group_by(Stage) %>%
  mutate(total_n = n()) %>%
  group_by(Stage, Ctype.Collapsed) %>%
  mutate(freq = n() / total_n) %>%
  distinct()

blast.freq.stage$Ctype.Collapsed <- factor(blast.freq.stage$Ctype.Collapsed,
                                            levels = c(AGGREGATE_ORDER, "T-NK"))

p.bar <- ggbarplot(blast.freq.stage, x = "Stage", y = "freq",
                    fill = "Ctype.Collapsed", size = 0.2, legend = "right")

ggsave(p.bar, device = "eps", width = 5, height = 4.5, units = "in",
       filename = file.path(dir.plots, "AggregateCellFrequency_byStage.eps"))


# =============================================================================
# SECTION 4: Cell-type-level pairwise comparisons
# =============================================================================
# For each aggregated cell type, plots all pairwise group comparisons
# (Subtype × Subtype and TrajCluster × TrajCluster) at both IDX and MRD.
# These supplement Section 3 with per-cell-type statistical detail.

ctype.list       <- AGGREGATE_ORDER
subtype.pairs    <- combn(unique(idx.ctype$Subtype), 2, simplify = FALSE)
cluster.pairs    <- combn(c(1, 2, 3, 4), 2, simplify = FALSE)

# Helper: generate per-cell-type pairwise comparison plots
plot_pairwise <- function(df.idx, df.mrd, df.lfc = NULL,
                          group.var, ctype.list, comparisons,
                          include.lfc = FALSE) {
  plots <- list()
  for (ct in ctype.list) {
    plots[[paste0(ct, "_IDX")]] <- df.idx %>%
      filter(Ctype.Aggregate == ct) %>%
      ggboxplot(x = group.var, y = "freq", color = group.var,
                title = paste("IDX Frequencies:", ct), add = "jitter") +
      stat_compare_means(comparisons = comparisons, hide.ns = "p")

    plots[[paste0(ct, "_MRD")]] <- df.mrd %>%
      filter(Ctype.Aggregate == ct) %>%
      ggboxplot(x = group.var, y = "freq", color = group.var,
                title = paste("MRD Frequencies:", ct), add = "jitter") +
      stat_compare_means(comparisons = comparisons, hide.ns = "p")

    if (include.lfc && !is.null(df.lfc)) {
      plots[[paste0(ct, "_LFC")]] <- df.lfc %>%
        filter(Ctype.Aggregate == ct) %>%
        ggboxplot(x = group.var, y = "lfc", color = group.var,
                  title = paste("LFC:", ct), add = "jitter") +
        stat_compare_means(comparisons = comparisons, hide.ns = "p")
    }
  }
  plots
}

# -- Pairwise by Subtype -------------------------------------------------------
cat("Plotting pairwise comparisons by Subtype...\n")
pdf(file.path(dir.plots, "CType_pairwiseComparisons_bySubtype.pdf"),
    height = 6, width = 6)
for (p in plot_pairwise(idx.ctype, mrd.ctype, group.var = "Subtype",
                         ctype.list = ctype.list,
                         comparisons = subtype.pairs)) {
  print(p)
}
dev.off()

# -- Pairwise by Trajectory Cluster -------------------------------------------
cat("Plotting pairwise comparisons by Trajectory Cluster...\n")
pdf(file.path(dir.plots, "CType_pairwiseComparisons_byTrajCluster.pdf"),
    height = 6, width = 6)
for (p in plot_pairwise(idx.ctype, mrd.ctype, Ctype_count_lfc,
                         group.var = "k_cluster",
                         ctype.list = ctype.list,
                         comparisons = cluster.pairs,
                         include.lfc = TRUE)) {
  print(p)
}
dev.off()

cat("\nAbundance analysis complete. All figures saved to:", dir.plots, "\n")
