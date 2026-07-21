# =============================================================================
# Dual Hybrid Batch 2 — Variant Cell Mapping, Fisher Tests, and QC
# =============================================================================
#
# Description:
#   Loads vartrix per-cell genotype matrices (consensus and coverage modes)
#   from each sample, maps variant calls onto Seurat single-cell RNA objects,
#   and generates comprehensive diagnostic, statistical, and QC outputs.
#
#   Pipeline sections:
#     1. Shared data loading — mutation tables, Seurat object, VCF labels
#     2. Per-sample Seurat preparation — subset, normalize, UMAP
#     3. Per-sample variant mapping — vartrix → Seurat metadata → UMAP plots,
#        contingency tables, Fisher's exact tests, per-lineage VAF,
#        coverage-weighted pseudobulk VAF, and QC metrics
#     4. Cross-sample summary — combined heatmaps, forest plots, VAF barplots,
#        expression-vs-genotyping correlations, and saturation curves
#
# Inputs (set in USER CONFIGURATION):
#   - Vartrix output directories (from 02_vartrix_job.sh)
#   - Global RNA Seurat RDS (post-labeling)
#   - Updated projection metadata CSV
#   - Filtered mutation TSV (from Mission Bio scDNA analysis)
#   - MRD variant panel VCF (for SNP row labels)
#
# Outputs:
#   - Per-sample UMAP + contingency table PDFs
#   - Per-sample and combined Fisher test CSV tables
#   - Per-sample QC PDFs (saturation, coverage, expression-vs-genotyping)
#   - Cross-sample summary PDF (heatmaps, forest plots, VAF barplots)
#   - Combined CSVs: Fisher results, genotyping QC, lineage VAF, weighted VAF
#
# Dependencies:
#   tidyverse, Matrix, vcfR, VariantAnnotation, snpStats, Seurat, spgs,
#   grid, gridExtra, data.table, ggrepel, ggpubr, scales
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(tidyverse)
library(Matrix)
library(vcfR)
library(VariantAnnotation)
library(snpStats)
library(Seurat)
library(spgs)
library(grid)
library(gridExtra)
library(data.table)
library(ggrepel)
library(ggpubr)

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Root directory for vartrix results (must match results_root from 02_vartrix_job.sh)
results_root <- "/path/to/Scripts/dual_hybrid_batch2"

# Output directory for all plots
plot_dir <- "/path/to/Plots"

# Sample IDs to process (must match subdirectory names in results_root)
sample_ids <- c("PAVEPS_IDX", "PAVEPS_MRD", "PAREZH_IDX", "PAXYMC_MRD", "3840_IDX")

# Mutation evidence tables (from Mission Bio scDNA analysis)
path.mut.dna  <- "/path/to/Tables/mut_evidence_onlyscDNA.csv"
path.mut.dual <- "/path/to/Tables/mut_evidence_MissionBioMultiome.csv"
path.mut.all  <- "/path/to/Tables/filtered_mutations_gt0_meta.tsv"

# Global RNA Seurat object (post-labeling; must have sample.name in metadata)
path.global.rds     <- "/path/to/R_Objects/merged.postLabeling.rds"
path.projection.csv <- "/path/to/Tables/UpdatedProjectionResults_AllCells.csv"

# MRD variant panel VCF (GRCh38; used for SNP row labels)
path.vcf <- "/path/to/ref/MRD_allFilteredMuts.GRCh38.vcf"

# Fisher test minimum alt cells (mutations with fewer alt cells are skipped)
MIN_ALT_CELLS <- 3

# Coverage-weighted VAF: minimum reads per cell to include
MIN_READS_PER_CELL <- 2

# Variants to highlight in filtered barplots (set after first run)
# Leave as NULL for all variants with alt evidence, or specify a vector:
VARIANTS_TO_PLOT <- c(
  "KRAS:c.38G>A",
  "NRAS:c.38G>A",
  "PTPN11:c.215C>T",
  "KRAS:c.35_36insAGG"
)

# Vartrix genotype encoding (do not change unless vartrix encoding changes)
GENO_LABELS <- c("0" = "No Call", "1" = "ref/ref",
                  "2" = "alt/alt", "3" = "alt/ref")

# Blast lineage assignment
MYELOID_CTYPES <- c("cDC", "DC-progenitor", "GMP", "MEP", "Monocyte", "pDC")
HSC_CTYPES     <- c("HSC/MPP", "LMPP")
NK_T_CTYPES    <- c("NK", "T-cell")
B_CTYPES       <- c("CLP", "Immature-B", "Mature-B-1", "Mature-B-2",
                     "Pre-B", "Pre-pro-B", "Pro-B")

# Lineage colour palette (used across all plots)
LINEAGE_COLS <- c("HSC-like"     = "#00BFC4",
                   "B-lineage"    = "#F8766D",
                   "Myeloid-like" = "#C77CFF",
                   "Normal"       = "#7CAE00",
                   "Healthy"      = "#7CAE00")


# ── Helper functions ──────────────────────────────────────────────────────────

# Blast lineage assignment using the cell type maps defined in USER CONFIGURATION
assign_blast_lineage <- function(meta_df) {
  meta_df %>%
    mutate(Blast_Lineage = case_when(
      Projected.Ctype.Wmethod %in% MYELOID_CTYPES & Level1.anno == "Blast" ~ "Myeloid-like",
      Projected.Ctype.Wmethod %in% HSC_CTYPES     & Level1.anno == "Blast" ~ "HSC-like",
      Projected.Ctype.Wmethod %in% NK_T_CTYPES    & Level1.anno == "Blast" ~ "NK_T-cell",
      Projected.Ctype.Wmethod %in% B_CTYPES       & Level1.anno == "Blast" ~ "B-lineage",
      .default = "Healthy"
    ))
}

# Fisher's exact test: enrichment of ALT alleles in blast compartments vs. healthy
# Returns a data.frame of test results with BH-adjusted p-values.
run_fisher_per_sample <- function(seurat_obj, snv_matrix, mut.filtered,
                                   mut.names.table, sample_id,
                                   min_alt_cells = MIN_ALT_CELLS) {

  # Define all pairwise comparisons to test
  comparisons <- list(
    list(group = "All Blasts",         ref = "Healthy",
         group_fn = function(bl) bl %in% c("B-lineage", "Myeloid-like", "HSC-like"),
         ref_fn   = function(bl) bl == "Healthy"),
    list(group = "B-lineage Blast",    ref = "Healthy",
         group_fn = function(bl) bl == "B-lineage",
         ref_fn   = function(bl) bl == "Healthy"),
    list(group = "Myeloid-like Blast", ref = "Healthy",
         group_fn = function(bl) bl == "Myeloid-like",
         ref_fn   = function(bl) bl == "Healthy"),
    list(group = "HSC-like Blast",     ref = "Healthy",
         group_fn = function(bl) bl == "HSC-like",
         ref_fn   = function(bl) bl == "Healthy"),
    list(group = "Myeloid-like Blast", ref = "B-lineage Blast",
         group_fn = function(bl) bl == "Myeloid-like",
         ref_fn   = function(bl) bl == "B-lineage"),
    list(group = "HSC-like Blast",     ref = "B-lineage Blast",
         group_fn = function(bl) bl == "HSC-like",
         ref_fn   = function(bl) bl == "B-lineage")
  )

  results <- list()

  for (j in seq_along(mut.filtered)) {
    if (!mut.filtered[j] %in% rownames(snv_matrix)) next

    tryCatch({
      mut.df <- snv_matrix[mut.filtered[j], ] %>% t()
      rownames(mut.df) <- paste0(sample_id, "_", rownames(mut.df), "-1")
      shared <- intersect(rownames(mut.df), colnames(seurat_obj))
      if (length(shared) == 0) next

      genotypes <- mut.df[shared, 1]
      meta      <- seurat_obj@meta.data[shared, ]
      has_cov   <- genotypes != 0
      meta      <- meta[has_cov, ]
      genotypes <- genotypes[has_cov]
      has_alt   <- genotypes %in% c(2, 3)

      if (sum(has_alt) < min_alt_cells) next
      mut_label <- mut.names.table$mutWithName[j]

      for (comp in comparisons) {
        is_group <- comp$group_fn(meta$Blast_Lineage)
        is_ref   <- comp$ref_fn(meta$Blast_Lineage)
        keep     <- is_group | is_ref
        if (sum(is_group) == 0 || sum(is_ref) == 0) next

        tab <- table(Group = is_group[keep], Alt = has_alt[keep])
        if (!all(dim(tab) == c(2, 2))) next

        ft <- fisher.test(tab)
        comp_label <- paste0(comp$group, " vs. ", comp$ref)
        key <- paste0(mut_label, "__", gsub(" ", "_", comp_label))

        results[[key]] <- data.frame(
          sample     = sample_id,
          mutation   = mut_label,
          comparison = comp_label,
          group      = comp$group,
          reference  = comp$ref,
          odds_ratio = ft$estimate,
          p_value    = ft$p.value,
          ci_lower   = ft$conf.int[1],
          ci_upper   = ft$conf.int[2],
          n_group_alt = sum(is_group[keep] & has_alt[keep]),
          n_group_ref = sum(is_group[keep] & !has_alt[keep]),
          n_ref_alt   = sum(!is_group[keep] & has_alt[keep]),
          n_ref_ref   = sum(!is_group[keep] & !has_alt[keep]),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      cat("  Fisher: skipping", mut.filtered[j], ":", conditionMessage(e), "\n")
    })
  }

  results_df <- bind_rows(results)
  if (nrow(results_df) > 0) {
    results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
    results_df$sig_label <- case_when(
      results_df$p_adj < 0.001 ~ "***",
      results_df$p_adj < 0.01  ~ "**",
      results_df$p_adj < 0.05  ~ "*",
      TRUE                     ~ "ns"
    )
  }
  results_df
}

# Forest plot of Fisher OR results for one sample
plot_fisher_forest <- function(fisher_df, sample_id) {
  if (nrow(fisher_df) == 0) return(NULL)

  plot_df <- fisher_df %>%
    mutate(or_plot       = pmin(odds_ratio, 500),
           ci_upper_plot = pmin(ci_upper, 1000),
           ci_lower_plot = pmax(ci_lower, 0.001),
           label         = paste0(mutation, "\n", comparison),
           panel         = if_else(grepl("vs\\. Healthy", comparison),
                                    "vs. Healthy Cells", "vs. B-lineage Blasts"))

  ggplot(plot_df, aes(x = or_plot, y = reorder(label, or_plot),
                       color = comparison)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower_plot, xmax = ci_upper_plot), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    geom_text(aes(label = sig_label), hjust = -0.5, size = 4, show.legend = FALSE) +
    scale_x_log10() +
    facet_wrap(~ panel, scales = "free_y", ncol = 1) +
    theme_bw() +
    theme(strip.text  = element_text(size = 11, face = "bold"),
          axis.text.y = element_text(size = 8),
          legend.position = "none") +
    labs(title    = paste0("Fisher's Exact Test — Somatic Mutation Enrichment\n", sample_id),
         subtitle = paste0("Mutations with >= ", MIN_ALT_CELLS,
                            " alt-carrying cells | * p.adj<0.05, ** p.adj<0.01, *** p.adj<0.001"),
         x = "Odds Ratio (log scale)", y = NULL)
}

# Per-lineage VAF (fraction of genotyped cells carrying ALT)
compute_lineage_vaf <- function(seurat_obj, snv_matrix, mut.filtered,
                                 mut.names.table, sample_id) {
  vaf_rows <- list()
  for (j in seq_along(mut.filtered)) {
    if (!mut.filtered[j] %in% rownames(snv_matrix)) next
    tryCatch({
      mut.df <- snv_matrix[mut.filtered[j], ] %>% t()
      rownames(mut.df) <- paste0(sample_id, "_", rownames(mut.df), "-1")
      shared <- intersect(rownames(mut.df), colnames(seurat_obj))
      if (length(shared) == 0) next
      genotypes <- mut.df[shared, 1]
      meta <- seurat_obj@meta.data[shared, ]
      has_cov <- genotypes != 0
      meta <- meta[has_cov, ]; genotypes <- genotypes[has_cov]
      has_alt <- genotypes %in% c(2, 3)
      mut_label <- mut.names.table$mutWithName[j]

      for (lin in c("B-lineage", "HSC-like", "Myeloid-like", "Healthy")) {
        cells <- meta$Blast_Lineage == lin
        n_lin <- sum(cells)
        if (n_lin == 0) next
        vaf_rows[[paste0(sample_id, "_", mut_label, "_", lin)]] <- data.frame(
          sample = sample_id, mutation = mut_label, variant_id = mut.filtered[j],
          lineage = lin, n_cells = n_lin, n_alt = sum(cells & has_alt),
          frac_mutated = sum(cells & has_alt) / n_lin, stringsAsFactors = FALSE
        )
      }
    }, error = function(e) cat("  VAF: skipping", mut.filtered[j], "\n"))
  }
  bind_rows(vaf_rows)
}

# Coverage-weighted pseudobulk VAF per lineage (from alt/ref count matrices)
compute_weighted_vaf <- function(seurat_obj, alt_matrix, ref_matrix,
                                  mut.filtered, mut.names.table, sample_id,
                                  min_reads = MIN_READS_PER_CELL) {
  vaf_rows <- list()
  for (j in seq_along(mut.filtered)) {
    if (!mut.filtered[j] %in% rownames(alt_matrix)) next
    tryCatch({
      alt_counts   <- as.numeric(alt_matrix[mut.filtered[j], ])
      ref_counts   <- as.numeric(ref_matrix[mut.filtered[j], ])
      total_counts <- alt_counts + ref_counts
      names(alt_counts) <- names(ref_counts) <- names(total_counts) <-
        paste0(sample_id, "_", colnames(alt_matrix), "-1")

      shared <- intersect(names(alt_counts), colnames(seurat_obj))
      if (length(shared) == 0) next
      alt_counts <- alt_counts[shared]; ref_counts <- ref_counts[shared]
      total_counts <- total_counts[shared]
      meta <- seurat_obj@meta.data[shared, ]
      has_cov <- total_counts >= min_reads
      meta <- meta[has_cov, ]; alt_counts <- alt_counts[has_cov]
      ref_counts <- ref_counts[has_cov]; total_counts <- total_counts[has_cov]
      mut_label <- mut.names.table$mutWithName[j]

      for (lin in c("B-lineage", "HSC-like", "Myeloid-like", "Healthy")) {
        cells <- meta$Blast_Lineage == lin; n_lin <- sum(cells)
        if (n_lin == 0) next
        lin_alt <- sum(alt_counts[cells]); lin_total <- sum(total_counts[cells])
        cell_vafs <- alt_counts[cells] / total_counts[cells]
        cell_vafs <- cell_vafs[is.finite(cell_vafs)]
        vaf_rows[[paste0(sample_id, "_", mut_label, "_", lin)]] <- data.frame(
          sample = sample_id, mutation = mut_label, variant_id = mut.filtered[j],
          lineage = lin, n_cells = n_lin,
          n_cells_with_alt = sum(alt_counts[cells] > 0),
          total_alt_reads = lin_alt, total_ref_reads = sum(ref_counts[cells]),
          total_reads = lin_total,
          mean_coverage = round(lin_total / n_lin, 1),
          pseudobulk_vaf = round(if (lin_total > 0) lin_alt / lin_total else 0, 4),
          mean_cell_vaf = round(mean(cell_vafs, na.rm = TRUE), 4),
          median_cell_vaf = round(median(cell_vafs, na.rm = TRUE), 4),
          frac_any_alt = round(sum(alt_counts[cells] > 0) / n_lin, 4),
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) cat("  Weighted VAF: skipping", mut.filtered[j], "\n"))
  }
  bind_rows(vaf_rows)
}

# Gene and UMI saturation curves (subsampling)
plot_saturation <- function(summary_file, sample_id,
                             subsample_fracs = seq(0.1, 1.0, by = 0.1),
                             n_reps = 3) {
  if (!file.exists(summary_file)) return(NULL)
  reads <- fread(summary_file, header = TRUE) %>%
    filter(CB != "" & CB != "." & gene != "" & gene != ".")
  if (nrow(reads) == 0) return(NULL)
  n_total <- nrow(reads)

  sat_rows <- list()
  for (frac in subsample_fracs) {
    for (rep in 1:n_reps) {
      sub <- reads[sample(n_total, round(n_total * frac), replace = FALSE), ]
      per_cell <- sub %>% group_by(CB) %>%
        summarise(n_reads = n(), n_genes = n_distinct(gene),
                  n_umis = n_distinct(UMI), .groups = "drop")
      sat_rows[[paste0(frac, "_", rep)]] <- data.frame(
        frac = frac, rep = rep,
        median_reads = median(per_cell$n_reads),
        median_genes = median(per_cell$n_genes),
        median_umis  = median(per_cell$n_umis),
        stringsAsFactors = FALSE
      )
    }
  }
  sat_avg <- bind_rows(sat_rows) %>% group_by(frac) %>%
    summarise(across(starts_with("median"), mean), .groups = "drop")

  p_gene <- ggplot(sat_avg, aes(x = median_reads, y = median_genes)) +
    geom_line(color = "#E69F00", linewidth = 1.2) + geom_point(color = "#E69F00") +
    theme_bw() + labs(title = paste0("Gene Saturation — ", sample_id),
                       x = "Median reads/cell", y = "Genes/cell")
  p_umi <- ggplot(sat_avg, aes(x = median_reads, y = median_umis)) +
    geom_line(color = "#8B008B", linewidth = 1.2) + geom_point(color = "#8B008B") +
    theme_bw() + labs(title = paste0("UMI Saturation — ", sample_id),
                       x = "Median reads/cell", y = "UMIs/cell")
  list(gene = p_gene, umi = p_umi, data = sat_avg)
}


# =============================================================================
# SECTION 1: Load shared data
# =============================================================================

cat("Loading shared data...\n")

mut.dna  <- read.csv(path.mut.dna, row.names = 1)
mut.dual <- read.csv(path.mut.dual, row.names = 1)
mut.all  <- read.delim(path.mut.all)

global.rna      <- readRDS(path.global.rds)
metadata_to_add <- read_csv(path.projection.csv) %>% assign_blast_lineage()
global.rna      <- AddMetaData(global.rna, metadata = metadata_to_add)

# VCF SNP labels (used as row names for vartrix matrices)
snps <- fread(path.vcf, sep = "\t", header = TRUE, skip = "#CHROM") %>%
  dplyr::select(`#CHROM`, POS, REF, ALT) %>%
  mutate(mut = paste0(`#CHROM`, ":", POS, ":", REF, ":", ALT)) %>%
  pull(mut)


# =============================================================================
# SECTION 2: Per-sample Seurat preparation
# =============================================================================

cat("Building per-sample Seurat objects...\n")

seurat.list <- setNames(lapply(sample_ids, function(sid) {
  s <- subset(global.rna, subset = sample.name == sid) %>%
    NormalizeData() %>% FindVariableFeatures() %>%
    ScaleData() %>% RunPCA(npcs = 50) %>%
    FindNeighbors() %>% RunUMAP(dims = 1:40) %>%
    subset(subset = Blast_Lineage != "NK_T-cell")
  s
}), sample_ids)


# =============================================================================
# SECTION 3: Per-sample variant mapping, Fisher tests, VAF, and QC
# =============================================================================

all_fisher_results <- list()
all_genotyping_qc  <- list()
all_exp_vs_gt      <- list()
all_lineage_vaf    <- list()
all_weighted_vaf   <- list()

for (sid in sample_ids) {

  cat("\n========================================\n")
  cat("  Processing:", sid, "\n")
  cat("========================================\n")

  vartrix_dir <- file.path(results_root, sid, "vartrix_results")

  # -- Patient-specific mutations -------------------------------------------
  mut.filtered <- mut.all %>%
    filter(multiome_sample == sid) %>%
    pull(scDNA_variant_GRCh38) %>% unique()

  if (length(mut.filtered) == 0) {
    cat("  No mutations found for", sid, "— skipping.\n"); next
  }

  mut.names.table <- data.frame(mutations = mut.filtered) %>%
    left_join(mut.all %>%
                dplyr::select(scDNA_variant_GRCh38, scDNA_chr_GRCh38,
                              scDNA_pos_GRCh38, VEP_SYMBOL, VEP_HGVSc,
                              VEP_Amino_acids) %>%
                distinct(),
              by = c("mutations" = "scDNA_variant_GRCh38")) %>%
    separate(VEP_HGVSc, into = c(NA, "codonChange"), sep = ":") %>%
    mutate(mutWithName = paste0(VEP_SYMBOL, ":", codonChange))

  current.rna <- seurat.list[[sid]]

  # -- Load vartrix consensus matrix ----------------------------------------
  consensus_file <- file.path(vartrix_dir, "vartrix_10x_consensus.mtx")
  if (!file.exists(consensus_file)) {
    cat("  Consensus matrix not found — skipping.\n"); next
  }

  snv_matrix <- as.data.frame(as.matrix(readMM(consensus_file)))
  barcodes   <- spgs::reverseComplement(
    read.table(file.path(vartrix_dir, "all.cell.barcodes.tsv"))$V1,
    case = "upper"
  )
  colnames(snv_matrix) <- barcodes

  # Handle variant count mismatch (vartrix may filter variants)
  if (nrow(snv_matrix) != length(snps)) {
    cat("  WARNING: row count mismatch; using vartrix variants file.\n")
    vartrix_snps <- read.table(
      file.path(vartrix_dir, "vartrix_10x_variants_consensus.tsv"))$V1
    rownames(snv_matrix) <- snps[vartrix_snps]
  } else {
    rownames(snv_matrix) <- snps
  }

  # -- Load vartrix coverage matrices (alt + ref) --------------------------
  alt_file <- file.path(vartrix_dir, "vartrix_alt.mtx")
  ref_file <- file.path(vartrix_dir, "vartrix_ref.mtx")
  has_coverage_mtx <- file.exists(alt_file) && file.exists(ref_file)

  if (has_coverage_mtx) {
    alt_mtx <- as.data.frame(as.matrix(readMM(alt_file)))
    ref_mtx <- as.data.frame(as.matrix(readMM(ref_file)))
    colnames(alt_mtx) <- colnames(ref_mtx) <- barcodes
    rownames(alt_mtx) <- rownames(ref_mtx) <- rownames(snv_matrix)

    weighted_vaf <- compute_weighted_vaf(
      current.rna, alt_mtx, ref_mtx,
      mut.filtered, mut.names.table, sid
    )
    if (nrow(weighted_vaf) > 0) all_weighted_vaf[[sid]] <- weighted_vaf
  }

  # -- Plot mutations on UMAP (one page per mutation) -----------------------
  pdf(file.path(plot_dir, paste0(sid, "_DualHybridBatch2_MappingToSeurat.pdf")),
      height = 14, width = 12)

  for (j in seq_along(mut.filtered)) {
    if (mut.filtered[j] %notin% rownames(snv_matrix)) next
    tryCatch({
      mut.df <- snv_matrix[mut.filtered[j], ] %>% t()
      mut.df[, mut.filtered[j]] <- GENO_LABELS[as.character(mut.df[, mut.filtered[j]])]
      rownames(mut.df) <- paste0(sid, "_", rownames(mut.df), "-1")
      current.rna <- AddMetaData(current.rna, mut.df)

      alt.barcodes <- names(mut.df[mut.df[, 1] %in% c("alt/alt", "alt/ref"), ])

      # Panel 1: genotype UMAP
      p1 <- DimPlot(current.rna, group.by = mut.filtered[j],
                     order = c("ref/ref", "alt/ref", "alt/alt"),
                     na.value = "gray90", shuffle = TRUE) +
        ggtitle(paste0(mut.names.table$mutWithName[j], "\n", sid))

      # Panel 2: per-cell VAF (or binary highlight if no coverage matrices)
      if (has_coverage_mtx && mut.filtered[j] %in% rownames(alt_mtx) &&
          sum(as.numeric(alt_mtx[mut.filtered[j], ])) > 0) {
        cell_alt   <- as.numeric(alt_mtx[mut.filtered[j], ])
        cell_ref   <- as.numeric(ref_mtx[mut.filtered[j], ])
        cell_total <- cell_alt + cell_ref
        cell_vaf   <- ifelse(cell_total > 0, cell_alt / cell_total, NA)
        vaf_df     <- data.frame(row.names = paste0(sid, "_", barcodes, "-1"),
                                  cell_vaf = cell_vaf)
        current.rna <- AddMetaData(current.rna, vaf_df)
        p2 <- FeaturePlot(current.rna, features = "cell_vaf",
                           order = TRUE, pt.size = 0.5) +
          scale_color_gradientn(
            colors   = c("grey90", "#FFEDA0", "#FD8D3C", "#E31A1C", "#800026"),
            limits   = c(0, 1), na.value = "grey90", name = "VAF") +
          ggtitle(paste0("Per-cell VAF: ", mut.names.table$mutWithName[j]))
      } else {
        p2 <- DimPlot(current.rna, cells.highlight = alt.barcodes) +
          ggtitle("Alt-carrying cells")
      }

      p3 <- DimPlot(current.rna, group.by = "Level1.anno")
      p4 <- DimPlot(current.rna, group.by = "Blast_Lineage")

      # Contingency tables (raw counts + fraction mutated)
      tab1 <- tableGrob(table(
        current.rna@meta.data[, c("Blast_Lineage", mut.filtered[j])]))
      tab2 <- tableGrob(table(
        current.rna@meta.data[, c("Level1.anno", mut.filtered[j])]))

      build_frac_table <- function(meta, group_col, geno_col) {
        df <- meta[, c(group_col, geno_col)]
        colnames(df) <- c("Group", "genotype")
        df <- df %>% filter(genotype != "No Call")
        if (nrow(df) == 0) return(tableGrob(data.frame(Note = "No genotyped cells"),
                                              rows = NULL))
        df %>% group_by(Group) %>%
          summarise(n_genotyped  = n(),
                    n_mutated    = sum(genotype %in% c("alt/alt", "alt/ref")),
                    frac_mutated = round(n_mutated / n_genotyped * 100, 1),
                    .groups      = "drop") %>%
          mutate(label = paste0(n_mutated, "/", n_genotyped,
                                " (", frac_mutated, "%)")) %>%
          tableGrob(rows = NULL)
      }

      tab_frac1 <- build_frac_table(current.rna@meta.data,
                                      "Blast_Lineage", mut.filtered[j])
      tab_frac2 <- build_frac_table(current.rna@meta.data,
                                      "Level1.anno", mut.filtered[j])

      print(ggarrange(p1, p2, p3, p4,
                       tab1, tab_frac1, tab2, tab_frac2,
                       ncol = 2, nrow = 4, heights = c(2, 2, 1.2, 1.2)))
    }, error = function(e) {
      cat("  SKIPPING", mut.filtered[j], ":", conditionMessage(e), "\n")
    })
  }

  # -- Fisher's exact tests -------------------------------------------------
  fisher_results <- run_fisher_per_sample(
    current.rna, snv_matrix, mut.filtered, mut.names.table, sid)

  if (nrow(fisher_results) > 0) {
    write.csv(fisher_results,
              file.path(results_root, sid, paste0(sid, "_fisher_results.csv")),
              row.names = FALSE)
    all_fisher_results[[sid]] <- fisher_results
    p_forest <- plot_fisher_forest(fisher_results, sid)
    if (!is.null(p_forest)) print(p_forest)
  }

  dev.off()

  # -- Per-lineage VAF ------------------------------------------------------
  lineage_vaf <- compute_lineage_vaf(
    current.rna, snv_matrix, mut.filtered, mut.names.table, sid)
  if (nrow(lineage_vaf) > 0) all_lineage_vaf[[sid]] <- lineage_vaf

  # -- Genotyping QC --------------------------------------------------------
  snv_filt <- snv_matrix[mut.filtered[mut.filtered %in% rownames(snv_matrix)], ,
                          drop = FALSE]
  rownames(snv_filt) <- mut.names.table$mutWithName[
    match(rownames(snv_filt), mut.names.table$mutations)]

  nocall.pct <- rowSums(snv_filt == 0) / ncol(snv_filt)
  refref.pct <- rowSums(snv_filt == 1) / ncol(snv_filt)
  alt.pct    <- rowSums(snv_filt == 2 | snv_filt == 3, na.rm = TRUE) / ncol(snv_filt)

  gt_qc <- data.frame(
    sample = sid, mutation = names(nocall.pct),
    pct_nocall = as.numeric(nocall.pct) * 100,
    pct_ref = as.numeric(refref.pct) * 100,
    pct_alt = as.numeric(alt.pct) * 100,
    pct_genotyped = (1 - as.numeric(nocall.pct)) * 100,
    n_cells_total = ncol(snv_filt), stringsAsFactors = FALSE
  )
  gt_qc$gene <- sapply(strsplit(gt_qc$mutation, ":"), `[`, 1)
  all_genotyping_qc[[sid]] <- gt_qc

  # -- Expression vs. genotyping correlation --------------------------------
  genes <- sapply(strsplit(names(refref.pct), ":"), `[`, 1)
  pctExp <- DotPlot(current.rna, features = unique(genes))$data
  pct.mtx <- as.data.frame(1 - nocall.pct) %>% rownames_to_column("gene") %>%
    separate(gene, sep = ":", into = c("gene", "mut")) %>%
    arrange(gene, desc(`1 - nocall.pct`)) %>% filter(!duplicated(gene))
  colnames(pct.mtx) <- c("gene", "mut", "pct_genotyped")
  pctExp <- left_join(pctExp, pct.mtx, by = c("features.plot" = "gene"))
  pctExp$sample <- sid
  all_exp_vs_gt[[sid]] <- pctExp

  # -- QC PDF (saturation, coverage, expression-vs-genotyping) ---------------
  sample_dir   <- file.path(results_root, sid)
  summary_tsv  <- file.path(sample_dir, paste0(sid, ".annotated.intersect.summary.tsv"))
  target_file  <- file.path(sample_dir, paste0(sid, ".target_stats.regions.bed.gz"))
  bait_file    <- file.path(sample_dir, paste0(sid, ".bait_stats.regions.bed.gz"))

  sat_plots <- plot_saturation(summary_tsv, sid)

  qc_pdf <- file.path(plot_dir, paste0(sid, "_DualHybridBatch2_QC.pdf"))
  pdf(qc_pdf, height = 5, width = 6)

  if (!is.null(sat_plots)) {
    print(ggarrange(sat_plots$gene, sat_plots$umi, ncol = 2))
  }

  # Mosdepth coverage plots (if available)
  if (file.exists(target_file) && file.exists(bait_file)) {
    targets <- fread(target_file, col.names = c("chrom","start","end","name","mean_depth"))
    baits   <- fread(bait_file,   col.names = c("chrom","start","end","name","mean_depth"))
    targets$set <- "Target"; baits$set <- "Bait"
    cov_df <- bind_rows(targets, baits) %>% mutate(region_len = end - start)

    log1p10 <- scales::trans_new("log1p10",
                                  transform = function(x) log10(x + 1),
                                  inverse   = function(x) 10^x - 1)

    print(ggplot(cov_df, aes(x = mean_depth, fill = set)) +
            geom_histogram(position = "identity", alpha = 0.4, bins = 60) +
            scale_x_continuous("Mean depth", trans = log1p10,
                               breaks = c(0, 1, 10, 100, 1000, 10000)) +
            ylab("Intervals") + theme_bw() +
            labs(title = paste0("Target vs Bait Coverage: ", sid)))

    # Depth vs expression correlation
    goi <- mut.all %>% filter(multiome_sample == sid) %>% pull(VEP_SYMBOL) %>% unique()
    goi_counts <- rowSums(current.rna@assays[["RNA"]]@counts)
    goi_counts <- goi_counts[names(goi_counts) %in% goi]
    goi_df <- data.frame(gene = names(goi_counts), gene_counts = as.numeric(goi_counts))
    target_cov <- cov_df %>% filter(set == "Target", name %in% goi) %>%
      left_join(goi_df, by = c("name" = "gene")) %>% distinct()
    target_cov$mean_depth <- target_cov$mean_depth + 1

    print(ggplot(target_cov, aes(x = mean_depth, y = gene_counts)) +
            geom_point(aes(color = name)) +
            geom_text_repel(aes(label = name, color = name), max.overlaps = 15) +
            scale_x_log10() + scale_y_log10() + theme_bw() +
            geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
            stat_cor() + theme(legend.position = "none") +
            labs(title = paste0("Depth vs Expression: ", sid),
                 x = "Depth at mutation sites", y = "scRNA counts"))
  }

  # Expression vs genotyping scatter
  print(ggplot(pctExp, aes(x = pct_genotyped, y = pct.exp)) +
          geom_point(aes(color = avg.exp.scaled), size = 2.5) +
          geom_text_repel(aes(label = features.plot)) +
          theme_bw() +
          geom_smooth(method = "lm", se = FALSE, linetype = "dashed") +
          stat_cor(label.x = 0.2, label.y = 25) +
          labs(title = paste0("% Genotyped vs % Expressing: ", sid),
               x = "% Cells Genotyped", y = "% Cells Expressing Transcript"))

  dev.off()
  cat("  Done:", sid, "\n")
}


# =============================================================================
# SECTION 4: Cross-sample summary
# =============================================================================

combined_fisher       <- bind_rows(all_fisher_results)
combined_gt_qc        <- bind_rows(all_genotyping_qc)
combined_exp_gt       <- bind_rows(all_exp_vs_gt)
combined_vaf          <- bind_rows(all_lineage_vaf)
combined_weighted_vaf <- bind_rows(all_weighted_vaf)

write.csv(combined_fisher,
          file.path(plot_dir, "AllSamples_Fisher_Results_Combined.csv"),
          row.names = FALSE)
write.csv(combined_gt_qc,
          file.path(plot_dir, "AllSamples_Genotyping_QC_Combined.csv"),
          row.names = FALSE)
write.csv(combined_vaf,
          file.path(plot_dir, "AllSamples_Lineage_VAF_Combined.csv"),
          row.names = FALSE)
write.csv(combined_weighted_vaf,
          file.path(plot_dir, "AllSamples_WeightedVAF_Combined.csv"),
          row.names = FALSE)

# Variant selection for filtered barplots
variants_with_alt <- combined_vaf %>% filter(n_alt >= 1) %>%
  distinct(sample, mutation, variant_id)
write.csv(variants_with_alt,
          file.path(plot_dir, "AllSamples_VariantsWithAltEvidence.csv"),
          row.names = FALSE)

# Apply manual variant filter if specified in USER CONFIGURATION
combined_vaf_plot <- combined_vaf %>% filter(n_alt >= 1) %>%
  semi_join(variants_with_alt, by = c("sample", "mutation", "variant_id"))
combined_vaf_plot_filt <- if (!is.null(VARIANTS_TO_PLOT)) {
  combined_vaf %>% filter(mutation %in% VARIANTS_TO_PLOT)
} else {
  combined_vaf_plot
}


# ── Cross-sample summary PDF ──────────────────────────────────────────────────

summary_pdf <- file.path(plot_dir, "AllSamples_DualHybridBatch2_Summary.pdf")
pdf(summary_pdf, height = 8, width = 12)

# Page 1: Genotyping rate heatmap
if (nrow(combined_gt_qc) > 0) {
  print(ggplot(combined_gt_qc, aes(x = sample, y = mutation, fill = pct_genotyped)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = paste0(round(pct_genotyped, 1), "%")), size = 2.5) +
    scale_fill_gradient(low = "white", high = "darkgreen", name = "% Genotyped") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6)) +
    labs(title = "Genotyping Rate per Mutation Across Samples", x = NULL, y = NULL))
}

# Page 2: Stacked bar — No Call / Ref / Alt breakdown
if (nrow(combined_gt_qc) > 0) {
  gt_long <- combined_gt_qc %>%
    dplyr::select(sample, mutation, pct_nocall, pct_ref, pct_alt) %>%
    pivot_longer(cols = c(pct_nocall, pct_ref, pct_alt),
                 names_to = "category", values_to = "pct") %>%
    mutate(category = recode(category, pct_nocall = "No Call",
                             pct_ref = "Ref/Ref", pct_alt = "Alt (het + hom)"))

  print(ggplot(gt_long, aes(x = mutation, y = pct, fill = category)) +
    geom_col() + facet_wrap(~ sample, ncol = 5) + coord_flip() +
    scale_fill_manual(values = c("No Call" = "gray80", "Ref/Ref" = "steelblue",
                                  "Alt (het + hom)" = "firebrick")) +
    theme_bw() + theme(axis.text.y = element_text(size = 7),
                        legend.position = "bottom") +
    labs(title = "Genotyping Breakdown", x = NULL, y = "% Cells", fill = NULL))
}

# Page 3: Expression vs genotyping (all samples)
if (nrow(combined_exp_gt) > 0) {
  print(ggplot(combined_exp_gt, aes(x = pct_genotyped, y = pct.exp)) +
    geom_point(aes(color = sample), size = 2, alpha = 0.7) +
    geom_text_repel(aes(label = features.plot, color = sample),
                    size = 2.5, max.overlaps = 10, show.legend = FALSE) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray30") +
    stat_cor(size = 3.5) + theme_bw() + theme(legend.position = "bottom") +
    labs(title = "Expression vs Genotyping Rate — All Samples",
         x = "% Genotyped", y = "% Expressing Transcript"))
}

# Page 4: Fisher OR heatmap vs. Healthy
if (nrow(combined_fisher) > 0) {
  vs_healthy <- combined_fisher %>%
    filter(grepl("vs\\. Healthy", comparison)) %>%
    mutate(log2_or = log2(pmin(odds_ratio, 500)),
           facet_label = gsub(" vs\\. Healthy", "", comparison))
  if (nrow(vs_healthy) > 0) {
    print(ggplot(vs_healthy, aes(x = sample, y = mutation, fill = log2_or)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sig_label), size = 4) +
      scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                            midpoint = 0, name = "log2(OR)") +
      facet_wrap(~ facet_label, scales = "free_y") + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "bold")) +
      labs(title = "Mutation Enrichment vs Healthy", x = NULL, y = NULL))
  }
}

# Page 5: Lineage VAF barplot (filtered)
if (nrow(combined_vaf_plot_filt) > 0) {
  plot_df <- combined_vaf_plot_filt %>%
    mutate(lineage = factor(lineage,
      levels = c("HSC-like", "B-lineage", "Myeloid-like", "Healthy"),
      labels = c("HSC-like", "B-lineage", "Myeloid-like", "Normal")))

  print(ggplot(plot_df, aes(x = mutation, y = frac_mutated, fill = lineage)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    facet_grid(cols = vars(sample), scales = "free_x", space = "free_x") +
    scale_fill_manual(values = LINEAGE_COLS, name = NULL) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                        strip.text = element_text(face = "bold"),
                        legend.position = "top") +
    labs(title = "Fraction of Cells with Mutated Transcripts",
         x = NULL, y = "Fraction mutated"))
}

# Page 6: Pseudobulk VAF vs consensus (filtered)
if (nrow(combined_weighted_vaf) > 0 && nrow(combined_vaf) > 0) {
  vaf_compare <- left_join(
    combined_vaf %>% dplyr::select(sample, mutation, lineage, frac_mutated),
    combined_weighted_vaf %>% dplyr::select(sample, mutation, lineage,
                                             pseudobulk_vaf, mean_coverage),
    by = c("sample", "mutation", "lineage")
  ) %>% filter(!is.na(pseudobulk_vaf))

  if (!is.null(VARIANTS_TO_PLOT))
    vaf_compare <- vaf_compare %>% filter(mutation %in% VARIANTS_TO_PLOT)

  print(ggplot(vaf_compare, aes(x = frac_mutated, y = pseudobulk_vaf,
                                 color = lineage, size = mean_coverage)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, linetype = "dashed", color = "gray40") +
    facet_wrap(~ sample) + scale_color_manual(values = LINEAGE_COLS) +
    theme_bw() + theme(strip.text = element_text(face = "bold")) +
    labs(title = "Consensus vs Coverage-Weighted VAF",
         x = "Fraction mutated (consensus)", y = "Pseudobulk VAF"))

  # Weighted VAF barplot
  weighted_plot <- combined_weighted_vaf %>%
    filter(n_cells_with_alt >= 1) %>%
    mutate(lineage = factor(lineage,
      levels = c("B-lineage", "HSC-like", "Myeloid-like", "Healthy"),
      labels = c("B-lineage", "HSC-like", "Myeloid-like", "Normal")))
  if (!is.null(VARIANTS_TO_PLOT))
    weighted_plot <- weighted_plot %>% filter(mutation %in% VARIANTS_TO_PLOT)

  if (nrow(weighted_plot) > 0) {
    print(ggplot(weighted_plot,
                 aes(x = mutation, y = pseudobulk_vaf, fill = lineage)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      facet_grid(cols = vars(sample), scales = "free_x", space = "free_x") +
      scale_fill_manual(values = LINEAGE_COLS, name = NULL) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          strip.text = element_text(face = "bold"),
                          legend.position = "top") +
      labs(title = paste0("Pseudobulk VAF (min ", MIN_READS_PER_CELL, " reads/cell)"),
           x = NULL, y = "Pseudobulk VAF"))
  }
}

# Page 7: Bubble plot — VAF + fraction mutated combined
if (nrow(combined_weighted_vaf) > 0 && nrow(combined_vaf) > 0) {
  bubble_df <- left_join(
    combined_vaf %>% dplyr::select(sample, mutation, lineage, frac_mutated, n_cells, n_alt),
    combined_weighted_vaf %>% dplyr::select(sample, mutation, lineage,
                                             pseudobulk_vaf, mean_coverage),
    by = c("sample", "mutation", "lineage")
  ) %>%
    filter(!is.na(pseudobulk_vaf), n_alt >= 1) %>%
    mutate(lineage = factor(lineage,
      levels = c("B-lineage", "HSC-like", "Myeloid-like", "Healthy"),
      labels = c("B-lineage", "HSC-like", "Myeloid-like", "Normal")))

  if (!is.null(VARIANTS_TO_PLOT))
    bubble_df <- bubble_df %>% filter(mutation %in% VARIANTS_TO_PLOT)

  if (nrow(bubble_df) > 0) {
    print(ggplot(bubble_df, aes(x = mutation, y = pseudobulk_vaf,
                                 size = frac_mutated, color = lineage)) +
      geom_point(alpha = 0.8) +
      facet_grid(cols = vars(sample), scales = "free_x", space = "free_x") +
      scale_size_continuous(name = "Fraction\nmutated", range = c(0, 10),
                             labels = scales::percent) +
      scale_color_manual(values = LINEAGE_COLS, name = NULL) +
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                          strip.text = element_text(face = "bold")) +
      labs(title = "Somatic Mutation Evidence: Pseudobulk VAF + Fraction Mutated",
           x = NULL, y = "Pseudobulk VAF"))
  }
}

dev.off()
cat("\nSummary PDF saved:", summary_pdf, "\n")
cat("All samples complete.\n")
