# =============================================================================
# IDX Arrest-State Signature Survival Analysis
# =============================================================================
#
# Description:
#   Derives IDX-specific gene signatures from DEGs that distinguish dASP
#   groups at diagnosis, then validates their prognostic relevance in an
#   independent B-ALL cohort (St. Jude / COG).
#
#   Pipeline steps:
#     1. Identify IDX-specific dASP signatures by FindAllMarkers on IDX blasts,
#        keeping genes unique to each dASP (not shared across dASPs)
#     2. Score the external COG cohort (bulk RNA-seq, rlog values) using GSVA
#     3. Run Cox regression (EFS) per signature adjusting for CNS status,
#        WBC count, MRD positivity, and risk level
#     4. Plot Kaplan-Meier curves and forest plots of hazard ratios
#     5. Validate TC_2 signature at the per-gene level (gene-level Cox EFS)
#     6. GO BP + Hallmark pathway enrichment on IDX signatures (compareCluster)
#     7. Correlation of IDX signature scores with arrest-state LFC changes
#     8. Correlation of TC_2 gene expression with HSPC-like LFC at MRD
#
# Inputs (set in USER CONFIGURATION):
#   - Pre-computed IDX dASP DEG RData file (tc.idx.degs)
#   - COG cohort Seurat object (bulk rlog RNA-seq values)
#   - HTAN blast Seurat object (for intra-cohort validation plots)
#   - Updated projection metadata CSV (for cell type frequency correlations)
#   - Cell type frequency table CSV (for correlation plots)
#   - AUCell gene set scores RData (for gene-level expression plots)
#   - Healthy donor RNA object (for GO enrichment universe)
#
# Outputs:
#   - Survival PDFs (KM curves, forest plots)
#   - GO/Hallmark enrichment dot plots
#   - Association plots (signature vs MRD group, arrest-state LFC, cell type)
#   - TC_2 IDX poor-outcome gene table CSV
#
# Dependencies:
#   Seurat, survival, ggfortify, cowplot, ggpubr, clusterProfiler,
#   msigdbr, org.Hs.eg.db, GSVA, AUCell, dplyr, tidyverse, data.table,
#   RColorBrewer, conflicted
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(tibble)
library(grid)
library(conflicted)
library(tidyr)
library(enrichplot)
library(survival)
library(ggfortify)
library(cowplot)
library(readxl)
library(ggpubr)
library(data.table)
library(scales)
library(stringr)
library(AUCell)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA)
library(SummarizedExperiment)
library(purrr)

conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("union",  "base")
conflicted::conflict_prefer("setdiff","base")

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ────────────────────────────────────────────────────────

# Pre-computed IDX dASP DEGs (FindAllMarkers on IDX blasts by dASP)
# If this file does not exist, set RERUN_IDX_DEGS = TRUE to re-run Section 1
path.idx.degs    <- "/path/to/Survival_Analysis/Objects/TC_IDX_DEGs.RData"
RERUN_IDX_DEGS   <- FALSE   # set TRUE to re-run FindAllMarkers

# HTAN blast Seurat object (IDX/MRD; used for Section 1 re-run and Sections 7–8)
path.global.rna  <- "/path/to/R_Objects/global.rna.blasts.onlyidxmrd.rds"

# Updated projection metadata CSV (adds Projected.Ctype.Wmethod, k_cluster, etc.)
path.projection.meta <- "/path/to/Tables/UpdatedProjectionResults_AllCells.csv"

# COG / St. Jude cohort Seurat object (bulk rlog RNA-seq)
path.cog.rds     <- "/path/to/Survival_Analysis/Objects/BALL_StJudesCohort_SeuratObject_RLogValues.rds"

# Cell type frequency table (from arrest-state abundance analysis)
path.celltype.freq <- "/path/to/Tables/CtypeAggFrequency_UpdatedAnno.csv"

# AUCell gene set scores for HTAN blasts (Hallmark gene sets, by dASP × stage)
path.aucell.rdata <- "/path/to/Tables/AUCell/geneSets_H_ExpFiltered_dsKmeansStage.RData"

# Healthy donor RNA object (used to define expression universe for GO)
path.hd.rna      <- "/path/to/SeuratObj/HDMultiome/global_downsample50K.rds"

# Output directories
dir.survival.plots <- "/path/to/Survival_Analysis/Plots/StJudes_COG/"
dir.tables         <- "/path/to/Survival_Analysis/Tables/"

# ── IDX DEG filtering thresholds ─────────────────────────────────────────────
IDX_LFC_CUTOFF  <- 0.25   # min avg_log2FC for IDX signature genes
IDX_PADJ_CUTOFF <- 0.05   # max p_val_adj
IDX_PCT_CUTOFF  <- 0.20   # min pct.1

# ── Survival analysis parameters ─────────────────────────────────────────────
UPPER_QUARTILE  <- 0.5    # top fraction for "high" signature group
LOWER_QUARTILE  <- 0.5    # bottom fraction for "low" signature group

# Risk group assignment: favorable vs. intermediate/unfavorable subtypes
FAVORABLE_SUBTYPES <- c("ETV6-RUNX1", "Hyperdiploid", "DUX4",
                          "NUTM1", "PAX5 P80R")
UNFAVORABLE_SUBTYPES <- c("ZNF384", "PAX5alt", "Ph-like_other",
                            "TCF3-PBX1", "Near haploid", "ETV6-RUNX1-like",
                            "Ph", "Ph-like_CRLF2", "KMT2A",
                            "Low hypodiploid", "MEF2D", "BCL2/MYC",
                            "HLF", "iAMP21")

# ── Enrichment parameters ─────────────────────────────────────────────────────
GO_PVAL         <- 0.2
GO_QVAL         <- 0.2
GO_SIMPLIFY     <- 0.4
ENRICHR_DBS     <- c("MSigDB_Hallmark_2020",
                      "GO_Biological_Process_2023",
                      "GO_Molecular_Function_2023",
                      "Reactome_2022")

# Pseudocount for LFC computation in correlation plots
LFC_PSEUDOCOUNT  <- 0.0003632401


# ── Shared helper functions ───────────────────────────────────────────────────

#' Score a Seurat object using GSVA
gsva_score <- function(seurat.obj, gene.sets) {
  counts <- GetAssayData(seurat.obj, assay = "RNA", layer = "counts")
  params <- gsvaParam(exprData = counts, geneSets = gene.sets)
  scores <- t(gsva(params, verbose = FALSE))
  AddMetaData(seurat.obj, scores)
}

#' Kaplan–Meier + Cox EFS plot for a single signature score.
#' Returns the Cox model coefficient table (invisibly).
km_efs_plot <- function(seurat.obj, score.name,
                          upper.quartile = UPPER_QUARTILE,
                          lower.quartile = LOWER_QUARTILE,
                          risk = FALSE, mrd = FALSE, scale = TRUE) {
  surv.data <- seurat.obj@meta.data
  surv.data$rank <- rank(surv.data[[score.name]])
  n <- nrow(surv.data)

  surv.data <- surv.data %>%
    mutate(high.mpp = case_when(
      rank > n * upper.quartile ~ "Sig High",
      rank < n * lower.quartile ~ "Sig Low"
    )) %>%
    filter(!is.na(high.mpp))

  km.fit   <- survfit(Surv(time.EFS, EFS.indicator) ~ high.mpp,
                       data = surv.data)
  cox.uni  <- coxph(Surv(time.EFS, EFS.indicator) ~ high.mpp,
                     data = surv.data)
  p.val    <- summary(cox.uni)$logtest[3]

  if (scale) surv.data$score_scaled <- scale(surv.data[[score.name]])[, 1]
  else       surv.data$score_scaled <- surv.data[[score.name]]

  # Multivariable Cox model with clinical covariates
  cox.multi <- if (risk) {
    coxph(Surv(time.EFS, EFS.indicator) ~
            score_scaled + CNS + WBC + Day.29.MRD.Positivity,
          data = surv.data)
  } else if (mrd) {
    coxph(Surv(time.EFS, EFS.indicator) ~
            score_scaled + CNS + WBC + Risk_Level,
          data = surv.data)
  } else {
    coxph(Surv(time.EFS, EFS.indicator) ~
            score_scaled + CNS + WBC + Day.29.MRD.Positivity + Risk_Level,
          data = surv.data)
  }

  cox.summ        <- summary(cox.multi)
  p.val.signature <- cox.summ$coefficients[1, 5]
  hr              <- cox.summ$coefficients[1, 2]

  p.km <- autoplot(km.fit) +
    labs(x = "EFS (Days)", y = "Survival Probability") +
    ylim(0, 1) + theme_classic() +
    scale_fill_manual(values  = c("#BE433D", "#459DA5")) +
    scale_color_manual(values = c("#BE433D", "#459DA5")) +
    geom_step(aes(color = group), alpha = 0.8) +
    geom_point(shape = "+") +
    annotate("text", label = paste0("p = ", round(p.val, 5)),
             x = 2500, y = 0.4) +
    annotate("text", label = paste0("p.sig = ", round(p.val.signature, 5)),
             x = 2500, y = 0.50) +
    ggtitle(score.name,
            subtitle = paste0("n = ", nrow(surv.data)))

  surv.box <- surv.data %>%
    filter(EFS.indicator %in% c(0, 1)) %>%
    mutate(event = if_else(EFS.indicator == 0, "No Event", "Event"))

  p.box <- ggplot(surv.box, aes(x = event, y = .data[[score.name]])) +
    geom_jitter() +
    geom_violin(draw_quantiles = 0.5, alpha = 0.5) +
    stat_compare_means() +
    theme_bw() +
    ggtitle(score.name, subtitle = paste0("n = ", nrow(surv.box))) +
    labs(y = score.name)

  print(cowplot::plot_grid(p.km, p.box))

  invisible(cox.summ$coefficients)
}

#' Scatter plot of IDX signature score vs. MRD cell type frequency (LFC)
plot_signature_lfc <- function(sig_name, meta_df, freq_df,
                                 pseudocount = LFC_PSEUDOCOUNT) {
  sig_col <- paste0(sig_name, "_avg_IDX")

  sig_sample <- meta_df %>%
    filter(Stage == "IDX") %>%
    group_by(PT.number) %>%
    summarise(!!sig_col := mean(.data[[sig_name]], na.rm = TRUE), .groups = "drop")

  freq_df$k_cluster <- as.factor(freq_df$k_cluster)

  plot_df <- freq_df %>%
    pivot_wider(id_cols     = c(PT.number, Ctype.Aggregate, k_cluster),
                names_from  = Stage, values_from = cell.freq,
                values_fill = 0) %>%
    mutate(LFC = log2((MRD + pseudocount) / (IDX + pseudocount))) %>%
    left_join(sig_sample, by = "PT.number")

  ggplot(plot_df, aes(y = LFC, x = .data[[sig_col]])) +
    geom_point(size = 1, aes(color = k_cluster)) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    stat_cor(aes(y = LFC), label.x.npc = 0.2, label.y.npc = 0.05,
             vjust = -1, size = 3, method = "spearman") +
    facet_wrap(~ Ctype.Aggregate, scales = "free", ncol = 7) +
    labs(y = "LFC (MRD/IDX cell freq)",
         x = paste(sig_name, "score (IDX)"),
         title = paste(sig_name, "IDX score vs arrest-state LFC")) +
    theme_bw()
}

#' Gene-level Pearson correlation with a cell-type-specific LFC vector
gene_corr_assoc <- function(gene_lfc_mat, celltype_lfc, sig_genes = NULL,
                              alpha_thresh = 0.05) {
  patients <- base::intersect(colnames(gene_lfc_mat), names(celltype_lfc))
  gene_lfc_mat <- as.matrix(gene_lfc_mat[, patients, drop = FALSE])
  y <- as.numeric(celltype_lfc[patients])
  if (is.null(sig_genes)) sig_genes <- rownames(gene_lfc_mat)
  sig_genes <- sig_genes[sig_genes %in% rownames(gene_lfc_mat)]

  res <- lapply(sig_genes, function(g) {
    x    <- as.numeric(gene_lfc_mat[g, ])
    keep <- is.finite(x) & is.finite(y)
    if (sum(keep) < 3) return(NULL)
    ct   <- cor.test(x[keep], y[keep], method = "pearson")
    data.frame(gene = g, R = unname(ct$estimate),
               p_val = ct$p.value, n = sum(keep))
  })

  df         <- bind_rows(res)
  df$p_adj   <- p.adjust(df$p_val, method = "fdr")
  df[order(df$p_adj), ]
}


# =============================================================================
# SECTION 1: Derive IDX-specific dASP signatures
# =============================================================================
# Identifies genes upregulated at IDX in each dASP relative to all other dASPs.
# Only genes unique to one dASP are retained to maximise specificity.

if (RERUN_IDX_DEGS) {
  global.rna       <- readRDS(path.global.rna)
  metadata_to_add  <- read.csv(path.projection.meta, row.names = 1)
  global.rna       <- AddMetaData(global.rna, metadata = metadata_to_add)

  Idents(global.rna) <- "Ctype.Aggregate"
  global.rna <- subset(global.rna, idents = "T-NK", invert = TRUE)

  idx.blasts  <- subset(global.rna, subset = Stage == "IDX")
  Idents(idx.blasts) <- "k_cluster"

  tc.idx.degs <- FindAllMarkers(idx.blasts,
                                  logfc.threshold     = 0,
                                  max.cells.per.ident = 20000,
                                  latent.vars         = "PT.number",
                                  test.use            = "MAST")
  save(tc.idx.degs, file = path.idx.degs)
  cat("IDX DEGs saved to:", path.idx.degs, "\n")
} else {
  load(path.idx.degs)
}

# Filter and build unique per-dASP signature lists
tc.idx.degs.filt <- tc.idx.degs %>%
  filter(avg_log2FC > IDX_LFC_CUTOFF,
         p_val_adj  < IDX_PADJ_CUTOFF,
         pct.1      > IDX_PCT_CUTOFF,
         !grepl("^MT-", gene))

tc.idx.degs.list <- split(tc.idx.degs.filt, f = tc.idx.degs.filt$cluster)
tc.idx.degs.list <- lapply(tc.idx.degs.list, pull, gene)

# Keep only genes unique to one dASP (not shared across dASPs)
unique.tc.list <- lapply(names(tc.idx.degs.list), function(nm) {
  others <- unlist(tc.idx.degs.list[setdiff(names(tc.idx.degs.list), nm)])
  setdiff(tc.idx.degs.list[[nm]], others)
})
names(unique.tc.list) <- paste0("TC_", names(tc.idx.degs.list))

cat("Unique IDX signature sizes:\n")
print(sapply(unique.tc.list, length))

# Write signature gene table for supplement
max_len    <- max(sapply(unique.tc.list, length))
padded     <- lapply(unique.tc.list, function(x) { length(x) <- max_len; x })
sig.df     <- as.data.frame(padded) %>%
  relocate(TC_1, TC_2, TC_3, TC_4)
write.table(sig.df,
            file.path(dir.tables, "TC_IDX_unique_signatures.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)


# =============================================================================
# SECTION 2: Set up and score the COG / St. Jude external cohort
# =============================================================================

cat("Loading COG cohort...\n")
cog.s <- readRDS(path.cog.rds)

# Clean metadata: strip redundant "0:No / 1:Yes" prefixes from factor columns
cols.to.clean <- c("Treatment.Arm", "WBC.count.at.baseline",
                    "Time.to.Relapse.early.intermediate.late.",
                    "CNS.status.at.baseline", "Bone.Marrow.status.at.baseline",
                    "First.event.type", "Day.8.Percent.MRD",
                    "Day.29.MRD.category")

cog.s@meta.data <- cog.s@meta.data %>%
  mutate(across(all_of(cols.to.clean),
                ~ str_split_fixed(., ":", 2)[, 2]),
         EFS.indicator = if_else(`If.had.a.relapse` == "0:No", 0, 1),
         time.EFS      = as.numeric(EFS.years.from.Dx) * 365,
         time.OS       = as.numeric(OS.years.from.Dx)  * 365,
         status.OS     = as.numeric(OS.indictor.from.Dx))

# Standardise key clinical variables
cog.s@meta.data <- cog.s@meta.data %>%
  mutate(
    Day.29.MRD.category = case_when(
      Day.29.Percent.MRD < 0.01                       ~ "<0.01%",
      Day.29.Percent.MRD >= 0.01 & Day.29.Percent.MRD < 0.1  ~ "0.01-0.1%",
      Day.29.Percent.MRD >= 0.1  & Day.29.Percent.MRD < 1    ~ "0.1-1%",
      Day.29.Percent.MRD >= 1    & Day.29.Percent.MRD < 5    ~ "1-5%",
      Day.29.Percent.MRD >= 5                         ~ ">5%"
    ),
    Day.29.MRD.Positivity = if_else(Day.29.Percent.MRD < 0.01,
                                     "Negative", "Positive"),
    CNS.status.at.baseline = na_if(CNS.status.at.baseline, "Unknown/missing"),
    WBC.count.at.baseline  = na_if(WBC.count.at.baseline,  "Unknown/missing"),
    WBC.count.at.baseline  = if_else(WBC.count.at.baseline == "<50",
                                      "<50 k", WBC.count.at.baseline),
    Risk_Level  = case_when(
      Primary.sample.subtype %in% FAVORABLE_SUBTYPES    ~ "Favorable",
      Primary.sample.subtype %in% UNFAVORABLE_SUBTYPES  ~ "Intermediate/Unfavorable"
    ),
    WBC = case_when(
      WBC.count.at.baseline == "<50 k"               ~ "<50K",
      WBC.count.at.baseline %in% c("50-100 k", ">100 k") ~ ">50K"
    ),
    CNS = case_when(
      CNS.status.at.baseline == "CNS 1"              ~ "CNS 1",
      CNS.status.at.baseline %in% c("CNS 2", "CNS 3") ~ "CNS 2+3"
    )
  )

# Set reference factor levels for Cox regression
cog.s$WBC                  <- relevel(factor(cog.s$WBC,  c("<50K", ">50K")),  ref = "<50K")
cog.s$CNS                  <- relevel(factor(cog.s$CNS,  c("CNS 1", "CNS 2+3")), ref = "CNS 1")
cog.s$Day.29.MRD.Positivity <- relevel(factor(cog.s$Day.29.MRD.Positivity,
                                               c("Negative", "Positive")), ref = "Negative")
cog.s$Risk_Level            <- relevel(factor(cog.s$Risk_Level,
                                               c("Favorable", "Intermediate/Unfavorable")),
                                        ref = "Favorable")

cat("COG cohort:", ncol(cog.s), "patients\n")
cat("Events (EFS):", sum(cog.s$EFS.indicator), "of", ncol(cog.s), "\n")

# Build COG patient subsets for stratified survival analysis
cog.mrdPos       <- subset(cog.s, subset = Day.29.MRD.Positivity == "Positive")
cog.hr           <- subset(cog.s, subset = Risk_Level == "Intermediate/Unfavorable")
cog.byTreatment  <- SplitObject(cog.s, "TreatmentProtocol")
cog.byTreatment[c("AALL08B1", "AALL0031")] <- NULL

cog.objs <- c(list(AllPts = cog.s,
                    MRDPositive = cog.mrdPos,
                    HighRiskSubtypes = cog.hr),
               cog.byTreatment)

# Score all COG subsets using GSVA
cat("Scoring with GSVA...\n")
cog.objs <- lapply(cog.objs, function(obj) gsva_score(obj, unique.tc.list))


# =============================================================================
# SECTION 3: EFS survival analysis per signature × patient group
# =============================================================================

score.names <- rev(names(unique.tc.list))
groups      <- c("All Patients", "MRD+ Patients", "High-Risk Subtypes")
efs.overall <- data.frame()

for (i in 1:3) {
  file.name <- file.path(dir.survival.plots,
                           paste0("EFS_SurvivalPlots_",
                                  names(cog.objs)[i], "_GSVA_IDXSigs.pdf"))

  pdf(file.name, height = 3, width = 7)

  efs.mtx <- matrix(NA, nrow = length(score.names), ncol = 5,
                     dimnames = list(names(unique.tc.list),
                                     c("coef", "exp(coef)", "se(coef)", "z", "Pr(>|z|)")))

  for (j in seq_along(score.names)) {
    score.name  <- score.names[j]
    is.mrd.grp  <- names(cog.objs)[i] %in% c("MRDPositive", "P9906")
    is.risk.grp <- names(cog.objs)[i] == "HighRiskSubtypes"

    efs.mtx[score.name, ] <- km_efs_plot(
      cog.objs[[i]], score.name,
      mrd  = is.mrd.grp,
      risk = is.risk.grp
    )[1, ]
  }

  dev.off()

  # Build forest plot data
  efs.df <- efs.mtx %>%
    as.data.frame() %>%
    dplyr::select("exp(coef)", "se(coef)", "Pr(>|z|)") %>%
    dplyr::rename(HR = "exp(coef)", HRse = "se(coef)", pvalue = "Pr(>|z|)") %>%
    mutate(HRup   = HR + 1.96 * HRse,
           HRlow  = HR - 1.96 * HRse,
           Association = case_when(
             HRlow > 1 & pvalue < 0.05 ~ "Adverse",
             HRup  < 1 & pvalue < 0.05 ~ "Favorable",
             TRUE                       ~ "N.S."
           ),
           PatientGroup      = groups[i],
           Signature_Name    = factor(rownames(efs.mtx),
                                      levels = c("TC_1","TC_2","TC_3","TC_4")))

  efs.overall <- rbind(efs.overall, efs.df)
}

p.forest <- efs.overall %>%
  ggplot(aes(x = Signature_Name, y = HR, ymax = HRup, ymin = HRlow,
             color = Association)) +
  geom_hline(yintercept = 1, lty = 2, color = "darkgrey") +
  geom_pointrange(size = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("Adverse" = "red", "N.S." = "black")) +
  facet_wrap(~ PatientGroup, ncol = 1) +
  theme(strip.background = element_rect(fill = "lightgray"))

pdf(file.path(dir.survival.plots,
               "UniqueIDX_TCSignatures_EFS_HRPlot_allGroups.pdf"),
    height = 5, width = 4)
p.forest
dev.off()


# =============================================================================
# SECTION 4: TC_2 signature per-gene Cox analysis
# =============================================================================
# Tests each gene in the TC_2 IDX signature for individual EFS association
# in the full COG cohort (FDR-adjusted, top 25%/bottom 25% split).

cat("Running per-gene Cox analysis for TC_2 signature...\n")
tc.2.sig     <- unique.tc.list[["TC_2"]]

# Temporarily downgrade RNA assay for signatureScore_zscore compatibility
cog.s[["RNA"]] <- as(object = cog.s[["RNA"]], Class = "Assay")

for (gene in tc.2.sig) {
  if (gene %notin% rownames(cog.s)) next
  cog.s <- signatureScore_zscore(cog.s, gene, score.name = gene)
}

genes.scored <- tc.2.sig[tc.2.sig %in% rownames(cog.s)]

pdf(file.path(dir.survival.plots, "COG_tc2sig_EFS_perGene.pdf"),
    height = 3, width = 7)

efs.gene.mtx <- matrix(NA, nrow = length(genes.scored), ncol = 5,
                         dimnames = list(genes.scored,
                                         c("coef","exp(coef)","se(coef)","z","Pr(>|z|)")))

for (j in seq_along(genes.scored)) {
  efs.gene.mtx[j, ] <- km_efs_plot(cog.s, genes.scored[j],
                                     upper.quartile = 0.75,
                                     lower.quartile = 0.25,
                                     scale = FALSE)[1, ]
}

dev.off()

efs.gene.df <- efs.gene.mtx %>%
  as.data.frame() %>%
  dplyr::select("exp(coef)", "se(coef)", "Pr(>|z|)") %>%
  dplyr::rename(HR = "exp(coef)", HRse = "se(coef)", pvalue = "Pr(>|z|)") %>%
  mutate(HRup  = HR + 1.96 * HRse,
         HRlow = HR - 1.96 * HRse,
         fdr   = p.adjust(pvalue, method = "fdr"),
         Association = case_when(
           HRlow > 1 & fdr < 0.05 ~ "Adverse",
           HRup  < 1 & fdr < 0.05 ~ "Favorable",
           TRUE                    ~ "N.S."
         )) %>%
  rownames_to_column("gene") %>%
  left_join(tc.idx.degs.filt %>% filter(cluster == 2),
            by = "gene")

write.csv(efs.gene.df,
          file.path(dir.tables, "TC2_IDX_poorOutcomeGenes.csv"),
          row.names = FALSE)

# Forest plot of significant TC_2 genes
p.tc2.genes <- efs.gene.df %>%
  filter(Association %in% c("Adverse", "Favorable")) %>%
  ggplot(aes(x = reorder(gene, -HR), y = HR,
             ymax = HRup, ymin = HRlow, color = Association)) +
  geom_hline(yintercept = 1, lty = 2, color = "darkgrey") +
  geom_pointrange(size = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("Adverse" = "red", "Favorable" = "black")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("TC_2 IDX genes — EFS association")

pdf(file.path(dir.survival.plots, "TC2_IDX_GeneLevel_EFS_ForestPlot.pdf"),
    height = 3, width = 17)
p.tc2.genes
dev.off()


# =============================================================================
# SECTION 5: GO BP + Hallmark enrichment on IDX signatures (compareCluster)
# =============================================================================

cat("Running pathway enrichment on IDX signatures...\n")

hd.rna  <- readRDS(path.hd.rna)
universe <- rownames(hd.rna@assays$RNA)

hallmark_results <- compareCluster(
  unique.tc.list,
  fun        = "enricher",
  TERM2GENE  = msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::select(gs_name, gene_symbol),
  pAdjustMethod = "BH",
  pvalueCutoff  = GO_PVAL,
  qvalueCutoff  = GO_QVAL,
  minGSSize     = 3,
  universe      = universe
)

go_bp_results <- compareCluster(
  unique.tc.list,
  fun           = "enrichGO",
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = GO_PVAL,
  qvalueCutoff  = GO_QVAL,
  readable      = TRUE,
  minGSSize     = 3,
  universe      = universe
)

go_bp_simp <- clusterProfiler::simplify(go_bp_results,
                                          cutoff = GO_SIMPLIFY,
                                          by = "p.adjust",
                                          select_fun = min,
                                          measure = "Wang")

# Combined enrichment data frame
combined_df <- rbind(
  as.data.frame(go_bp_simp)  %>% mutate(source = "GO_BP"),
  as.data.frame(hallmark_results) %>% mutate(source = "Hallmark")
) %>%
  mutate(logP = pmin(-log10(as.numeric(p.adjust)), 10))

write.csv(combined_df,
          file.path(dir.tables,
                     "GOBP_Hallmark_combined_IDXSignatures.csv"),
          row.names = FALSE)

# Filter to top 5 per cluster + shared terms
top5 <- combined_df %>%
  group_by(Cluster) %>%
  slice_max(logP, n = 5, with_ties = FALSE) %>%
  ungroup()

shared_terms <- combined_df %>%
  group_by(Description) %>%
  summarise(n_clusters = n_distinct(Cluster)) %>%
  filter(n_clusters > 1) %>%
  pull(Description)

plot_df <- combined_df %>%
  filter(Description %in% union(shared_terms, top5$Description)) %>%
  arrange(Cluster, qvalue) %>%
  mutate(rank = row_number(),
         Cluster = factor(Cluster, levels = rev(names(unique.tc.list))))

p.enrich <- ggplot(plot_df,
                    aes(x = Cluster, y = reorder(Description, -rank),
                        size = Count, color = logP)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "blue", high = "red",
                        name = "-log10(p.adj)", limits = c(0, 10)) +
  labs(title = "GO BP & Hallmark — IDX dASP Signatures",
       x = "Signature", y = "Enriched Terms") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(size = 10, angle = 45, hjust = 1))

ggsave(file.path(dir.survival.plots,
                  "GOBP_Hallmark_IDXSignatures_enrichment.pdf"),
       p.enrich, width = 7, height = 6)


# =============================================================================
# SECTION 6: IDX signature score vs. arrest-state LFC (correlation plots)
# =============================================================================

cat("Building signature–LFC correlation plots...\n")

# Load HTAN blast RNA with AUCell scores
load(path.aucell.rdata)

global.rna      <- readRDS(path.global.rna)
metadata_to_add <- read.csv(path.projection.meta, row.names = 1)
global.rna      <- AddMetaData(global.rna, metadata = metadata_to_add)

kmeans.subset   <- subset(global.rna,
                            subset = patient_bc %in%
                              geneSets_H_allBlast@colData@rownames)

auc_mat <- t(getAUC(geneSets_H_allBlast))
kmeans.subset <- AddMetaData(kmeans.subset, as.data.frame(auc_mat))

# Compute AUCell scores for IDX signatures on HTAN blasts
expr_mat <- GetAssayData(kmeans.subset, assay = "RNA", slot = "data")
rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE)
auc_sig  <- AUCell_calcAUC(geneSets = unique.tc.list, rankings = rankings)
kmeans.subset <- AddMetaData(kmeans.subset,
                               as.data.frame(t(assay(auc_sig, "AUC"))))

celltype.freq <- read.csv(path.celltype.freq, row.names = 1)

sig.lfc.plots <- imap(unique.tc.list,
                       ~ plot_signature_lfc(.y, kmeans.subset@meta.data,
                                             celltype.freq))

pdf(file.path(dir.survival.plots, "TC_IDX_SignatureLFC_Correlations.pdf"),
    height = 8, width = 14)
ggarrange(plotlist = sig.lfc.plots[c("TC_1", "TC_2", "TC_3")],
           nrow = 3, common.legend = TRUE, legend = "right")
dev.off()


# =============================================================================
# SECTION 7: TC_2 gene expression vs. HSPC-like LFC at MRD (per-gene correlation)
# =============================================================================

cat("Computing TC_2 gene-level correlations with HSPC LFC...\n")

# Average expression per patient at IDX
avg.per.patient <- AverageExpression(
  object    = subset(kmeans.subset, subset = Stage == "IDX"),
  assays    = "RNA",
  group.by  = "PT.number",
  return.seurat = FALSE
)$RNA

# HSPC-like LFC vector per patient
lfc.mtx   <- sig.lfc.plots[["TC_2"]] %>%
  filter(Ctype.Aggregate == "HSPC") %>%
  pull(value, PT.number)

# Per-gene Pearson correlation
sig.genes <- unique.tc.list[["TC_2"]]
patients  <- names(lfc.mtx)

gene_corr_df <- gene_corr_assoc(
  gene_lfc_mat  = avg.per.patient[, patients, drop = FALSE],
  celltype_lfc  = lfc.mtx,
  sig_genes     = sig.genes[sig.genes %in% rownames(avg.per.patient)]
)

sig.genes.for.plot <- gene_corr_df %>%
  filter(!is.na(R), !is.na(p_adj)) %>%
  slice_max(R, n = 30) %>%
  mutate(sig = factor(p_adj < 0.1,
                       levels = c(FALSE, TRUE),
                       labels = c("p_adj >= 0.1", "p_adj < 0.1")))

pdf(file.path(dir.survival.plots,
               "TC2_IDX_HSPC_LFC_GeneAssociation.pdf"),
    height = 6, width = 4)

ggplot(sig.genes.for.plot, aes(x = R, y = reorder(gene, R))) +
  geom_vline(xintercept = 0, lty = 2, color = "gray50") +
  geom_point(aes(size = -log10(p_adj), color = sig), shape = 16) +
  scale_color_manual(values = c("p_adj >= 0.1" = "steelblue",
                                 "p_adj < 0.1"  = "firebrick")) +
  scale_size(name = "-log10(p_adj)") +
  labs(x     = "Pearson R (gene expression vs HSPC LFC)",
       y     = "TC_2 IDX gene",
       title = "TC_2 gene-level correlation\nwith HSPC-like arrest-state LFC") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 8))

dev.off()

cat("All survival analyses complete.\n")
