# =============================================================================
# PacBio Hybrid Capture Pilot — Variant Cell Mapping and QC
# =============================================================================
#
# Description:
#   Loads vartrix per-cell genotype matrices from each demultiplexed barcode,
#   maps variant calls onto Seurat single-cell RNA objects, and generates
#   diagnostic and QC plots for the hybridization capture pilot experiment.
#
#   Pipeline sections:
#     1. Mutation reference preparation — builds cleaned VCF CSV and a
#        human-readable mutation name table from the filtered mutation list
#     2. Seurat object preparation — subsets the global RNA object by patient
#        sample and runs PCA/UMAP for per-sample visualization
#     3. Variant-to-cell mapping — reads vartrix consensus matrices, maps
#        genotype calls onto Seurat UMAPs, and saves per-barcode PDFs
#     4. Vartrix call rate QC — histograms of no-call, ref/ref, and alt %
#        per mutation per barcode
#     5. Gene expression vs. ref/ref correlation — scatter plots comparing
#        per-gene expression (DotPlot) against vartrix ref/ref frequency
#     6. IsoSeq performance QC — paired plots of FLNC % and primer orientation
#        % across probe concentrations
#     7. Target region coverage — cumulative coverage curves across mutation
#        loci from bedtools coverage histograms
#     8. Gene mapping summary — bar charts of unique read-to-gene mapping
#        counts per barcode
#     9. qPCR target expression — DotPlot of off-target transcript expression
#        in blast cells by stage
#
# Inputs (set in USER CONFIGURATION):
#   - Global RNA Seurat RDS (post-labeling, with Level1.anno and Blast_Lineage)
#   - Projection metadata CSV (UpdatedProjectionResults)
#   - Filtered mutation TSV (from Mission Bio scDNA analysis)
#   - Barcode-to-sample CSV (pacbio_pilot_barcodes.csv)
#   - VCF of MRD variant panel (GRCh38, sorted)
#   - IsoSeq performance stats CSV
#   - Coverage histogram files (*.hist.all.txt, from bedtools coverage -hist)
#   - Per-barcode vartrix output directories (containing vartrix_consensus.mtx
#     and all.cell.barcodes.tsv)
#
# Outputs:
#   - Cleaned mutation VCF CSV
#   - Per-barcode mutation UMAP PDFs
#   - Vartrix call rate QC PDF
#   - IsoSeq QC and expression correlation PDF
#   - Coverage curve PNG
#   - Gene mapping bar chart PDF
#   - qPCR target expression PDF
#
# Dependencies:
#   tidyverse, Matrix, vcfR, VariantAnnotation, snpStats, Seurat, spgs,
#   grid, gridExtra, data.table, RColorBrewer, ggpubr
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
library(RColorBrewer)
library(ggpubr)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Working directory — must contain one subdirectory per barcode with vartrix
# output files (vartrix_consensus.mtx, all.cell.barcodes.tsv)
cwd <- "/path/to/Scripts/hybrid_pilot"

# Barcode-to-biosample mapping CSV
# Expected columns: Barcodes, Bio.Sample
path.barcode.csv <- file.path(cwd, "pacbio_pilot_barcodes.csv")

# Filtered mutation file (tab-separated) from Mission Bio scDNA pipeline
# Expected columns: scDNA_variant_GRCh38, scDNA_chr_GRCh38, scDNA_pos_GRCh37,
#   scDNA_pos_GRCh38, scDNA_ref, scDNA_alt, VEP_SYMBOL, VEP_variant,
#   VEP_HGVSc, multiome_sample
path.mut.all <- "/path/to/MissionBio/Tables/filtered_mutations_gt0_meta.tsv"

# Global RNA Seurat object (post-labeling; requires Level1.anno, Blast_Lineage,
# sample.name metadata columns)
path.global.rds <- "/path/to/R_Objects/merged.postLabeling.rds"

# Updated projection metadata CSV (adds Projected.Ctype.Wmethod etc.)
path.projection.csv <- "/path/to/Tables/UpdatedProjectionResults_AllCells.csv"

# Sorted VCF of MRD variant panel (GRCh38); used to build SNP row labels
path.vcf.sorted <- "/path/to/ref/MRD_variant_list.GRCh38.sorted.vcf"

# IsoSeq pilot performance stats CSV
# Expected columns: Samples, Reads.input, Correct.primer.orientation,
#   FLNC.with.poly.A, Stage
path.isoseq.stats <- "/path/to/Tables/PacBio_HybridPilot_PerformanceStats.csv"

# Patient samples corresponding to each barcode folder (in the same order as
# barcode.folders). Maps barcode folder index → patient sample ID.
# Each patient sample may appear more than once (e.g. at different probe
# concentrations). Format: named list, name = barcode label, value = sample ID.
sample.id.conversion <- list(
  "PAXFFK_IDX_150" = "PAXFFK_IDX",
  "PAXFFK_MRD_150" = "PAXFFK_MRD",
  "PAXFFK_IDX_300" = "PAXFFK_IDX",
  "PAXFFK_MRD_300" = "PAXFFK_MRD"
)

# Patient samples to filter mutations for (multiome_sample column values)
PATIENT_SAMPLES_OF_INTEREST <- c("PAXFFK_IDX", "PAXFFK_MRD")

# Samples to exclude from IsoSeq QC paired plots (e.g. failed runs)
ISOSEQ_QC_EXCLUDE <- c("IL7R-PT14-IDX")

# Coverage histogram file indices to use (subset of all hist.all files found).
# Index into the sorted list returned by list.files().
COV_FILE_INDICES <- c(2, 4, 6, 8)

# Sample labels for coverage plot legend (in the same order as COV_FILE_INDICES)
COV_LABELS <- c("PAXFFK_IDX_150", "PAXFFK_MRD_150",
                 "PAXFFK_IDX_300", "PAXFFK_MRD_300")

# Probe concentrations corresponding to each sample in sample.id.conversion
# (used for IsoSeq QC paired plot x-axis)
CONCENTRATIONS <- c(150, 150, 300, 300)

# qPCR target genes to plot expression for in blast cells
QPCR_TARGETS <- c("NRAS", "CREBBP", "NT5C2", "MALAT1", "RPL13A", "MT-ATP6")

# Vartrix genotype encoding (do not change unless vartrix encoding changes)
VARTRIX_CODES <- c("0" = "No Call", "1" = "ref/ref",
                   "2" = "alt/alt",  "3" = "alt/ref")

# Output paths
dir.plots      <- "/path/to/Plots/"
path.vcf.out   <- "/path/to/ref/all.filtered.mut.cleaned.csv"


# =============================================================================
# SECTION 1: Mutation reference preparation
# =============================================================================
# Loads the full filtered mutation table and builds:
#   (a) A cleaned VCF CSV replacing NA gene symbols and '*' alt alleles with
#       VCF-compatible '.' placeholders — used as input to vartrix/cellsnp-lite
#   (b) A human-readable mutation name table (SYMBOL:codonChange) for plot labels

cat("Loading mutation data...\n")
mut.all <- read.delim(path.mut.all)

# Build cleaned VCF CSV (single output; NA and '*' replaced with '.')
mut.csv <- mut.all %>%
  dplyr::select(scDNA_chr_GRCh38, scDNA_pos_GRCh37, scDNA_pos_GRCh38,
                scDNA_ref, scDNA_alt, VEP_SYMBOL, VEP_variant) %>%
  distinct() %>%
  mutate(
    VEP_SYMBOL = if_else(is.na(VEP_SYMBOL), ".", VEP_SYMBOL),
    scDNA_alt  = if_else(scDNA_alt == "*",   ".", scDNA_alt)
  )

write.csv(mut.csv, path.vcf.out, row.names = FALSE, quote = FALSE)
cat("Cleaned mutation VCF CSV written to:", path.vcf.out, "\n")

# Subset mutations for patient(s) of interest
mut.filtered <- mut.all %>%
  dplyr::filter(multiome_sample %in% PATIENT_SAMPLES_OF_INTEREST) %>%
  pull(scDNA_variant_GRCh38) %>%
  unique()

# Build human-readable mutation name table for plot titles
mut.names.table <- data.frame(mutations = mut.filtered) %>%
  left_join(
    mut.all %>%
      dplyr::select(scDNA_variant_GRCh38, VEP_SYMBOL, VEP_HGVSc) %>%
      distinct(),
    by = c("mutations" = "scDNA_variant_GRCh38")
  ) %>%
  separate(VEP_HGVSc, into = c(NA, "codonChange"), sep = ":") %>%
  mutate(mutWithName = paste0(VEP_SYMBOL, ":", codonChange))


# =============================================================================
# SECTION 2: Seurat object preparation
# =============================================================================
# Loads the global RNA object, adds projection metadata, then subsets and
# re-reduces each unique patient sample for per-sample UMAP visualization.
# One Seurat object is created per unique sample (not per barcode), since
# multiple barcodes can map to the same patient sample (different concentrations).

cat("Loading Seurat object...\n")
global.rna      <- readRDS(path.global.rds)
metadata_to_add <- read_csv(path.projection.csv)
global.rna      <- AddMetaData(global.rna, metadata = metadata_to_add)

cat("Preparing per-sample Seurat objects...\n")
unique.samples <- unique(unlist(sample.id.conversion))

seurat.objs <- setNames(
  lapply(unique.samples, function(sname) {
    s <- subset(global.rna, subset = sample.name == sname)
    s <- s %>%
      NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA(npcs = 50) %>%
      FindNeighbors() %>%
      RunUMAP(dims = 1:40)
    s
  }),
  unique.samples
)


# =============================================================================
# SECTION 3: Variant-to-cell mapping
# =============================================================================
# For each barcode folder:
#   1. Reads the vartrix consensus matrix and cell barcode list
#   2. Reverse-complements barcodes (PacBio long-read barcodes are stored as
#      reverse complement relative to 10x Genomics Seurat cell names)
#   3. Assigns VCF-derived SNP IDs as row names
#   4. For each mutation of interest, adds genotype calls to the Seurat object
#      metadata and plots a 4-panel UMAP (genotype / alt highlight /
#      Level1.anno / Blast_Lineage) plus contingency tables

cat("Loading VCF SNP labels...\n")
snps <- fread(path.vcf.sorted, sep = "\t", header = TRUE, skip = "#CHROM") %>%
  dplyr::select(`#CHROM`, POS, REF, ALT) %>%
  mutate(mut = paste0(`#CHROM`, ":", POS, ":", REF, ":", ALT)) %>%
  pull(mut)

barcode.folders <- list.files(cwd, pattern = "^bc")
sample.names    <- read.csv(path.barcode.csv)

for (i in seq_along(barcode.folders)) {
  pacbio.barcode <- barcode.folders[i]
  patient.sample <- sample.id.conversion[[i]]
  cat("\nMapping variants for:", pacbio.barcode, "→", patient.sample, "\n")

  setwd(file.path(cwd, pacbio.barcode))

  # Load vartrix genotype matrix
  # vartrix_10x_consensus.mtx: output from vartrix run on 10x-format barcodes
  snv_matrix <- readMM("vartrix_10x_consensus.mtx") %>%
    as.matrix() %>%
    as.data.frame()

  # Reverse-complement barcodes to match Seurat cell names
  barcodes <- read.table("all.cell.barcodes.tsv", header = FALSE)
  barcodes <- spgs::reverseComplement(barcodes$V1, case = "upper")

  colnames(snv_matrix) <- barcodes
  rownames(snv_matrix) <- snps

  out.pdf <- file.path(dir.plots,
                       paste0(names(sample.id.conversion)[i],
                              "_HybridPilot_PacBio_MappingToSeurat.pdf"))
  pdf(out.pdf, height = 8, width = 10)

  for (j in seq_along(mut.filtered)) {
    mut.id   <- mut.filtered[j]
    mut.label <- mut.names.table$mutWithName[j]

    # Extract and recode genotype calls for this mutation
    mut.df <- t(snv_matrix[mut.id, , drop = FALSE]) %>% as.data.frame()
    mut.df[[mut.id]] <- recode(as.character(mut.df[[mut.id]]),
                               !!!VARTRIX_CODES)

    # Prefix barcodes with sample name to match Seurat cell names
    rownames(mut.df) <- paste0(seurat.objs[[patient.sample]]$sample.name[1],
                                "_", rownames(mut.df), "-1")

    seurat.objs[[patient.sample]] <- AddMetaData(
      seurat.objs[[patient.sample]], mut.df
    )

    alt.barcodes <- rownames(mut.df)[mut.df[[1]] %in% c("alt/alt", "alt/ref")]

    print(
      (DimPlot(seurat.objs[[patient.sample]],
               group.by  = mut.id,
               order     = c("ref/ref", "alt/ref", "alt/alt"),
               na.value  = "gray90") +
         ggtitle(paste0(mut.label, "\n", pacbio.barcode))) +
      DimPlot(seurat.objs[[patient.sample]],
              cells.highlight = alt.barcodes) +
      DimPlot(seurat.objs[[patient.sample]], group.by = "Level1.anno") +
      DimPlot(seurat.objs[[patient.sample]], group.by = "Blast_Lineage")
    )

    tab1 <- tableGrob(table(
      seurat.objs[[patient.sample]]@meta.data[, c("Blast_Lineage", mut.id)]
    ))
    tab2 <- tableGrob(table(
      seurat.objs[[patient.sample]]@meta.data[, c("Level1.anno", mut.id)]
    ))
    grid.arrange(tab1, tab2, ncol = 2)
  }

  dev.off()
  cat("Saved:", out.pdf, "\n")
}


# =============================================================================
# SECTION 4: Vartrix call rate QC
# =============================================================================
# For each barcode, computes the fraction of cells with no call, ref/ref, and
# alt allele detection per mutation. Histograms are written to a single PDF.
# Uses vartrix_consensus.mtx (standard vartrix output, not 10x-format).

pdf(file.path(dir.plots, "HybridPilot_VartrixCallFrequency.pdf"),
    height = 2.5, width = 7)

for (i in seq_along(barcode.folders)) {
  pacbio.barcode <- barcode.folders[i]
  patient.sample <- sample.id.conversion[[i]]
  sample.label   <- sample.names$Bio.Sample[
    match(pacbio.barcode, sample.names$Barcodes)
  ]

  setwd(file.path(cwd, pacbio.barcode))

  # vartrix_consensus.mtx: standard vartrix output (not 10x-format)
  snv_matrix <- readMM("vartrix_consensus.mtx") %>%
    as.matrix() %>%
    as.data.frame()

  barcodes <- read.table("all.cell.barcodes.tsv", header = FALSE)
  barcodes <- spgs::reverseComplement(barcodes$V1, case = "upper")

  colnames(snv_matrix) <- barcodes
  rownames(snv_matrix) <- snps

  snv_matrix           <- snv_matrix[mut.filtered, , drop = FALSE]
  rownames(snv_matrix) <- mut.names.table$mutWithName

  nocall.pct <- rowSums(snv_matrix == 0) / ncol(snv_matrix)
  refref.pct <- rowSums(snv_matrix == 1) / ncol(snv_matrix)
  alt.pct    <- rowSums(snv_matrix == 2 | snv_matrix == 3,
                         na.rm = TRUE) / ncol(snv_matrix)

  par(mfrow = c(1, 3))
  hist(nocall.pct, main = paste0(sample.label, " No Call %"), xlim = c(0.6, 1))
  hist(refref.pct, main = paste0(sample.label, " ref/ref %"), xlim = c(0, 0.4))
  hist(alt.pct,    main = paste0(sample.label, " alt %"),     xlim = c(0, 0.001))
}

dev.off()


# =============================================================================
# SECTION 5: Gene expression vs. ref/ref correlation
# =============================================================================
# Compares per-gene RNA expression (pct.exp from DotPlot) against the fraction
# of cells with ref/ref calls in the vartrix matrix for each mutation's gene.
# refref.pct is carried over from the last barcode in Section 4; if a
# per-barcode breakdown is needed, move this block inside the Section 4 loop.

genes <- sapply(strsplit(names(refref.pct), ":"), `[`, 1)

# Build ref/ref percent tables keyed by gene (keeping highest refref.pct per
# gene when the same gene has multiple mutations)
build_ref_pct_table <- function(refref.pct) {
  as.data.frame(refref.pct) %>%
    rownames_to_column("gene") %>%
    separate(gene, sep = ":", into = c("gene", "mut")) %>%
    arrange(gene, desc(refref.pct)) %>%
    filter(!duplicated(gene))
}

idx.ref.pct.mtx <- build_ref_pct_table(refref.pct)
mrd.ref.pct.mtx <- build_ref_pct_table(refref.pct)   # same source; update if per-sample

idx.pctExp <- DotPlot(seurat.objs[[unique.samples[1]]], features = unique(genes))$data
mrd.pctExp <- DotPlot(seurat.objs[[unique.samples[2]]], features = unique(genes))$data

idx.pctExp <- left_join(idx.pctExp, idx.ref.pct.mtx, by = c("features.plot" = "gene"))
mrd.pctExp <- left_join(mrd.pctExp, mrd.ref.pct.mtx, by = c("features.plot" = "gene"))


# =============================================================================
# SECTION 6: IsoSeq performance QC
# =============================================================================
# Plots FLNC % and correct primer orientation % comparing probe concentrations
# across samples. Samples in ISOSEQ_QC_EXCLUDE are dropped before plotting.

stats <- read.csv(path.isoseq.stats) %>%
  mutate(
    PercentCorrectPrimers = Correct.primer.orientation / Reads.input,
    PercentFLNC           = FLNC.with.poly.A / Correct.primer.orientation,
    Concentration         = CONCENTRATIONS
  )

pdf(file.path(dir.plots, "PilotPacBio_QC_ComparisonByConcentration.pdf"),
    height = 5, width = 6)

stats %>%
  dplyr::filter(Samples %notin% ISOSEQ_QC_EXCLUDE) %>%
  ggpaired(x = "Concentration", y = "PercentFLNC",
           id = "Stage", color = "Concentration", line.color = "gray",
           title = "FLNC % Comparison between Concentrations") +
  stat_compare_means()

stats %>%
  dplyr::filter(Samples %notin% ISOSEQ_QC_EXCLUDE) %>%
  ggpaired(x = "Concentration", y = "PercentCorrectPrimers",
           id = "Stage", color = "Concentration", line.color = "gray",
           title = "Correct Primers % Comparison between Concentrations") +
  stat_compare_means()

ggscatter(idx.pctExp, x = "pct.exp", y = "refref.pct",
          label = "features.plot", repel = TRUE, color = "avg.exp.scaled",
          font.label = c(9, "plain"),
          title = "Gene Expression vs. ref/ref % — IDX")

ggscatter(mrd.pctExp, x = "pct.exp", y = "refref.pct",
          label = "features.plot", repel = TRUE, color = "avg.exp.scaled",
          font.label = c(9, "plain"),
          title = "Gene Expression vs. ref/ref % — MRD")

dev.off()


# =============================================================================
# SECTION 7: Target region coverage curves
# =============================================================================
# Reads per-sample bedtools coverage histograms and plots cumulative coverage
# curves (fraction of loci profiled ≥ depth) on a single panel.
# COV_FILE_INDICES selects the relevant files from all hist.all files found.

cov.files <- list.files(cwd, pattern = "hist.all", recursive = TRUE,
                         full.names = TRUE)[COV_FILE_INDICES]

cov       <- lapply(cov.files, read.table)
cov_cumul <- lapply(cov, function(x) 1 - cumsum(x[, 5]))

cols <- brewer.pal(length(cov), "Dark2")

png(file.path(dir.plots, "CoverageAcrossMutations_HybridPilot.png"),
    h = 1000, w = 1000, pointsize = 20)

plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400],
     type = "n",
     xlab = "Reads Across Loci (Whole Sample)",
     ylab = "Fraction of total loci profiled \u2265 depth",
     ylim = c(0, 1.0),
     main = "Target Region Coverage")

abline(v = 250,  col = "gray60")
abline(h = c(0.10, 0.20, 0.50, 0.90), col = "gray60")
axis(1, at = 250,  labels = 250)
axis(2, at = c(0.10, 0.20, 0.50, 0.90),
     labels = c(0.10, 0.20, 0.50, 0.90))

for (i in seq_along(cov)) {
  points(cov[[i]][2:401, 2], cov_cumul[[i]][1:400],
         type = "l", lwd = 3, col = cols[i])
}

legend("topright", legend = COV_LABELS, col = cols, lty = 1, lwd = 4)
dev.off()


# =============================================================================
# SECTION 8: Gene mapping summary
# =============================================================================
# For each barcode, reads the gene overlap summary file produced by the
# bedtools intersect step in iso-seq.sh and plots the top 20 mapped genes
# alongside the top 20 genes that overlap the mutation target panel.
# Mitochondrial and ribosomal genes are collapsed to "mito" and "ribo".

summary.files <- list.files(cwd, pattern = "annotated.intersect.geneOverlap.txt",
                              recursive = TRUE, full.names = TRUE)

pdf(file.path(dir.plots, "HybridPilot_MappingOfTargetsBySample.pdf"),
    height = 5, width = 10)

for (i in seq_along(summary.files)) {
  mapping <- read.table(summary.files[i], header = FALSE)
  colnames(mapping) <- c("n", "gene")

  mapping <- mapping %>%
    mutate(gene = case_when(
      grepl("^MT",  gene) ~ "mito",
      grepl("^RPL", gene) ~ "ribo",
      grepl("^RPS", gene) ~ "ribo",
      TRUE                ~ gene
    )) %>%
    group_by(gene) %>%
    dplyr::summarize(n = sum(n), .groups = "drop") %>%
    arrange(desc(n))

  make_bar <- function(data, title) {
    ggplot(data %>% dplyr::slice(1:20),
           aes(x = reorder(gene, n), y = n, label = n)) +
      geom_col() +
      geom_text(hjust = -0.25) +
      coord_flip() +
      theme_bw() +
      labs(x = "Gene", y = "Count", title = title)
  }

  p1 <- make_bar(mapping,
                  paste0(names(sample.id.conversion)[i], " — All Genes"))
  p2 <- make_bar(mapping %>% dplyr::filter(gene %in% mut.all$VEP_SYMBOL),
                  paste0(names(sample.id.conversion)[i], " — Target Genes Only"))

  print(ggarrange(p1, p2, ncol = 2))
}

dev.off()


# =============================================================================
# SECTION 9: qPCR target expression in blast cells
# =============================================================================
# Plots expression of non-specific transcripts (potential qPCR off-targets)
# in blast cells by disease stage to assess contamination or background signal.

blasts.rna <- subset(global.rna, subset = Level1.anno == "Blast")

pdf(file.path(dir.plots, "Expression_qPCRTargets_byStage.pdf"),
    height = 4, width = 6)

DotPlot(blasts.rna, features = QPCR_TARGETS,
        group.by = "Stage", scale = FALSE) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  ggtitle("Expression of qPCR Targets by Stage")

dev.off()

cat("Pipeline complete.\n")
