# =============================================================================
# Transcriptional Regulatory Network (TRN) Construction
# Myeloid-like vs B-lineage Blast Cells by Stage
# =============================================================================
#
# Description:
#   Constructs enhancer-to-promoter (EP) transcriptional regulatory networks
#   for Myeloid-like and B-lineage blast cells at each stage (IDX, MRD)
#   using linked gene-peak regression on metacell-smoothed RNA and ATAC data.
#
#   The pipeline runs in three major steps:
#
#   STEP 1 — Metacell construction
#     Builds metacells separately for RNA and ATAC modalities using hdWGCNA,
#     grouping by blast lineage and subtype. Metacells are normalised and
#     filtered for genes/peaks expressed above a minimum fraction threshold.
#
#   STEP 2 — Gene-peak regression and motif scanning
#     Regresses gene expression against peaks within 500 kb, adjusting for
#     subtype. Identifies peak-gene links via FDR-corrected linear regression.
#     Scans peaks for TF motif matches using matchMotifs (JASPAR/CISBP).
#     Runs chromVAR and FindMarkers (LR test) to identify differential TFs.
#
#   STEP 3 — TRN assembly
#     Constructs per-cell-state regulons by intersecting:
#       - Predicted EP peak-gene links (from Step 2 regression)
#       - TF motif-peak matches
#       - Differentially expressed TFs and target genes (DEGs)
#     Writes regulon TSVs and Metascape-formatted gene lists per cell state.
#
# Usage:
#   The outer loop iterates over stages (IDX, MRD). Each call to trn()
#   runs the full pipeline for one stage and writes all outputs to a
#   stage-specific subdirectory under Scripts/TRN/.
#
#   To run in parallel across stages, replace the for loop at the bottom
#   with: bplapply(ctypes.list, trn, BPPARAM = MulticoreParam())
#
# Inputs (set in USER CONFIGURATION at the top of the trn() function):
#   - Combined multiome Seurat object (RNA + ATAC assays)
#   - Gene annotation text file (hg38 TSS positions)
#   - hg38 genome annotation RDS (for TSS and GRange operations)
#   - Signac hg38 annotation RDS (for peak annotation)
#
# Outputs (written to Scripts/TRN/<stage>/):
#   - <stage>_rna_metacell_<date>.rds       : smoothed RNA metacell object
#   - <stage>_atac_metacell_<date>.rds      : smoothed ATAC metacell object
#   - filtered_gene_peak_links_metacell_<date>.tsv  : EP regression links
#   - gene_peak_enhancerpeaksonly_links_metacell_<date>.tsv
#   - diff_TFactivity_<date>.rds            : differential TF motif activity
#   - degs_<date>.RDS                       : DEGs per blast lineage
#   - motif_peak_match_1e-05_<date>.rds     : TF motif × peak binary matrix
#   - <stage>.EP_Prediction_TRN_useDEGs_TF_regulons_*.tsv  : per-state regulons
#   - <stage>.EP_Prediction_TRN_useDEGs_nodes_*.tsv         : node tables
#   - <stage>.*_MetascapeInput.csv          : Metascape-formatted gene lists
#   - DEGs.csv                              : DEG table for Metascape
#
# Dependencies:
#   Seurat, Signac, hdWGCNA, WGCNA, chromVAR, chromVARmotifs, motifmatchr,
#   JASPAR2020, TFBSTools, BSgenome.Hsapiens.UCSC.hg38, data.table,
#   GenomicRanges, clusterProfiler, org.Hs.eg.db, BiocParallel, dplyr,
#   patchwork, RColorBrewer, Biobase, cowplot, tidyverse, conflicted
#
#   External function file:
#     EP_TRN_function.R  — provides get_gene2peak_map(), annPeak2Gene(),
#                           geneOverlapPeak() from the lab TRN pipeline
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(pheatmap)
  library(dplyr)
  library(patchwork)
  library(RColorBrewer)
  library(data.table)
  library(Biobase)
  library(GenomicRanges)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(conflicted)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(BiocParallel)
  library(tidyverse)
  library(cowplot)
  library(WGCNA)
  library(hdWGCNA)
  library(chromVAR)
  library(chromVARmotifs)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
})

conflicted::conflicts_prefer(base::intersect)
conflicted::conflicts_prefer(base::union)
conflicted::conflicts_prefer(base::setdiff)
conflicted::conflicts_prefer(dplyr::select)
conflicted::conflicts_prefer(dplyr::filter)

`%notin%` <- Negate(`%in%`)
theme_set(theme_cowplot())
set.seed(2023)
options(error = traceback)


# =============================================================================
# TRN pipeline function
# =============================================================================
# All parameters that vary by run are set in the USER CONFIGURATION block
# inside this function. The outer loop calls trn("IDX") and trn("MRD").

trn <- function(ctype_TC) {

  cat("\n=== Starting TRN pipeline for stage:", ctype_TC, "===\n")
  cat(format(Sys.time()), "\n")

  # ── USER CONFIGURATION (inside function) ────────────────────────────────────

  current.wd   <- "/path/to/Myeloid-Like_Blasts"
  current.date <- "013125"     # date tag used in all output filenames

  # Combined multiome Seurat object (RNA + ATAC assays, blast cells only,
  # B-lineage and Myeloid-like lineages, IDX or MRD stage)
  path.multiome <- "/path/to/Objects/CombinedTC12_RNA_ATAC_CISBP_chromVAR_onlyBlasts_recallPeaksbyBlastLineage.rds"

  # Gene annotation (hg38 TSS positions; tab-separated: gene_name, chr, strand, start, end)
  path.gene.ann <- "/path/to/Scripts/TRN/gene_ann_hg38.txt"

  # hg38 annotation RDS for Signac (used in peak-to-gene annotation)
  path.hg38.ann <- "/path/to/data/annotation_hg38.rds"

  grouping.ident   <- "Blast_Lineage"  # metadata column defining cell groups
  vars.to.regress  <- "Subtype"        # covariate to regress in regression

  # Metacell parameters
  k_metacell  <- 15    # nearest-neighbours for metacell construction
  peak.frac   <- 0.05  # min fraction of cells a peak must be accessible in
  gene.frac   <- 0.05  # min fraction of cells a gene must be expressed in
  ep.dist     <- 5e5   # enhancer-promoter search window (500 kb)
  motif.pv    <- 1e-5  # p-value cutoff for motif matching

  # TRN construction parameters
  includeNegativeDEGs  <- FALSE
  filterWeakPeak       <- TRUE
  useTopDscoreTF       <- FALSE   # use top chromVAR d-score TFs (vs. diff TFs)
  useDAP               <- FALSE   # filter peaks by DAPs (vs. open-fraction cutoff)
  useAllDEGs           <- TRUE    # use all DEGs or strict sig-only subset
  useEnhancerPeakOnly  <- TRUE    # exclude promoter-promoter interactions
  peakFracCutoff       <- 0.1     # min peak-open fraction per cell state
  tfexprfracCutoff     <- 0.10    # min TF gene expression fraction per state

  # ── Output directory ────────────────────────────────────────────────────────

  stage.dir <- file.path(current.wd, "Scripts/TRN", ctype_TC)
  if (!dir.exists(stage.dir)) {
    dir.create(stage.dir, recursive = TRUE)
    cat("Created output directory:", stage.dir, "\n")
  }
  setwd(stage.dir)

  # Source external TRN helper functions (get_gene2peak_map, annPeak2Gene, etc.)
  source("/path/to/EP_TRN_function.R")


  # ===========================================================================
  # STEP 1: Load multiome object and subset to current stage
  # ===========================================================================

  cat("Loading multiome object...\n")
  multi.obj <- readRDS(path.multiome)

  # Recompute ATAC QC metrics removed during object construction
  multi.obj$nCount_ATAC   <- colSums(multi.obj@assays$ATAC@counts)
  multi.obj$nFeature_ATAC <- colSums(multi.obj@assays$ATAC@counts > 0)

  # Subset to the current stage and to blast lineages of interest
  Idents(multi.obj) <- "Stage"
  multi.ctype.obj <- subset(multi.obj, idents = ctype_TC) %>%
    subset(subset = Blast_Lineage %in% c("B-lineage", "Myeloid-like"))

  cat("Cells in stage", ctype_TC, ":", ncol(multi.ctype.obj), "\n")
  print(table(multi.ctype.obj$Blast_Lineage))


  # ===========================================================================
  # STEP 1b: RNA dimensionality reduction and hdWGCNA setup
  # ===========================================================================

  multi.ctype.obj <- FindVariableFeatures(multi.ctype.obj,
                                           selection.method = "vst",
                                           assay = "RNA")
  multi.ctype.obj <- RunPCA(multi.ctype.obj, assay = "RNA")

  multi.ctype.obj <- RunTFIDF(multi.ctype.obj, assay = "ATAC")
  multi.ctype.obj <- FindTopFeatures(multi.ctype.obj, min.cutoff = 100,
                                      assay = "ATAC")
  multi.ctype.obj <- RunSVD(multi.ctype.obj,
                              features = VariableFeatures(multi.ctype.obj),
                              assay    = "ATAC")
  multi.ctype.obj <- FindMultiModalNeighbors(
    multi.ctype.obj,
    reduction.list = list("pca", "lsi"),
    dims.list      = list(1:40, 2:40)
  )
  multi.ctype.obj <- RunUMAP(multi.ctype.obj, nn.name = "weighted.nn",
                               reduction.name = "UMAP_joint")

  # Set up hdWGCNA on variable genes expressed in >= gene.frac of cells
  DefaultAssay(multi.ctype.obj) <- "RNA"
  multi.ctype.obj <- FindVariableFeatures(multi.ctype.obj, nfeatures = 2000)
  vegs   <- VariableFeatures(multi.ctype.obj)
  gfreq  <- rowMeans(multi.ctype.obj@assays$RNA@counts > 0)
  vegs   <- setdiff(vegs, names(which(gfreq < gene.frac)))

  multi.ctype.obj <- SetupForWGCNA(
    multi.ctype.obj,
    gene_select = "fraction",
    fraction    = gene.frac,
    wgcna_name  = "misc",
    features    = vegs
  )


  # ===========================================================================
  # STEP 1c: Metacell construction (RNA and ATAC separately)
  # ===========================================================================

  cat("Building ATAC metacells...\n")
  DefaultAssay(multi.ctype.obj) <- "ATAC"
  multi.ctype.obj <- MetacellsByGroups(
    seurat_obj  = multi.ctype.obj,
    min_cells   = k_metacell,
    group.by    = c("Blast_Lineage", "Subtype"),
    reduction   = "UMAP_joint",
    k           = k_metacell,
    max_shared  = 3,
    ident.group = "Blast_Lineage",
    wgcna_name  = "misc",
    verbose     = TRUE
  )
  metacell.atac <- GetMetacellObject(multi.ctype.obj)
  cat("ATAC metacells:", ncol(metacell.atac), "\n")

  cat("Building RNA metacells...\n")
  DefaultAssay(multi.ctype.obj) <- "RNA"
  multi.ctype.obj <- MetacellsByGroups(
    seurat_obj  = multi.ctype.obj,
    min_cells   = k_metacell,
    group.by    = c("Blast_Lineage", "Subtype"),
    reduction   = "UMAP_joint",
    k           = k_metacell,
    max_shared  = 3,
    ident.group = "Blast_Lineage"
  )
  metacell.rna <- GetMetacellObject(multi.ctype.obj)
  cat("RNA metacells:", ncol(metacell.rna), "\n")

  # Normalise metacell matrices
  metacell.atac <- JoinLayers(metacell.atac) %>% RunTFIDF(assay = "ATAC")
  metacell.rna  <- JoinLayers(metacell.rna)  %>% NormalizeData()

  smooth.rna  <- as.matrix(GetAssayData(metacell.rna,  assay = "RNA",  layer = "data"))
  smooth.atac <- as.matrix(GetAssayData(metacell.atac, assay = "ATAC", layer = "data"))
  stopifnot(all(colnames(smooth.rna) == colnames(smooth.atac)))

  saveRDS(metacell.rna,
          paste0(ctype_TC, "_rna_metacell_",  current.date, ".rds"))
  saveRDS(metacell.atac,
          paste0(ctype_TC, "_atac_metacell_", current.date, ".rds"))


  # ===========================================================================
  # STEP 1d: Peak and gene filtering for regression
  # ===========================================================================

  cat("Filtering peaks and genes...\n")
  peaks.mean.ctype <- sapply(unique(multi.ctype.obj$Blast_Lineage), function(x)
    rowMeans(multi.ctype.obj@assays$ATAC@counts[, multi.ctype.obj$Blast_Lineage == x] > 0)
  )
  filtered.peaks <- names(which(apply(peaks.mean.ctype, 1, max) > peak.frac))

  Gene.mean.ctype <- sapply(unique(multi.ctype.obj$Blast_Lineage), function(x)
    rowMeans(multi.ctype.obj@assays$RNA@data[, multi.ctype.obj$Blast_Lineage == x] > 0)
  )
  final.genes <- names(which(apply(Gene.mean.ctype, 1, max) > gene.frac))
  final.genes <- intersect(final.genes, rownames(multi.ctype.obj@assays$RNA@counts))

  metacell.atac <- metacell.atac[rownames(metacell.atac) %in% filtered.peaks, ]
  metacell.rna  <- metacell.rna[final.genes, ]
  saveRDS(filtered.peaks, paste0("filtered_peaks_", current.date, ".rds"))


  # ===========================================================================
  # STEP 2a: Gene-peak affinity map (500 kb window)
  # ===========================================================================

  cat("Building gene-peak affinity map...\n")
  gene_ann <- fread(path.gene.ann)
  gene_ann[, Tss := ifelse(strand == "+", start, end)]

  gene2peak.map <- get_gene2peak_map(
    gene.list   = data.table(gene_name = final.genes),
    peak_names  = filtered.peaks,
    gene_ann    = gene_ann,
    distal_dist = ep.dist
  )

  peaks.used  <- colnames(gene2peak.map)
  sele.genes  <- final.genes[final.genes %in% rownames(gene2peak.map)]
  smooth.rna  <- as.matrix(smooth.rna)[sele.genes, ]
  smooth.atac <- as.matrix(smooth.atac)[peaks.used, ]

  save(metacell.rna, metacell.atac, gene2peak.map,
       file = paste0("rna_atac_metacell_matrices_", current.date, ".RData"))


  # ===========================================================================
  # STEP 2b: Gene-peak regression (EP links)
  # ===========================================================================

  cat("Running gene-peak regression...\n")
  smooth.atac.df <- data.frame(t(smooth.atac))
  smooth.rna.df  <- data.frame(t(smooth.rna))

  # Rename peak columns: "." → "-" (R colname conversion reversal)
  peaknames <- colnames(smooth.atac.df)
  peaknames <- gsub(".", "-", peaknames, fixed = TRUE)
  colnames(smooth.rna.df)  <- sele.genes
  colnames(smooth.atac.df) <- peaknames

  metacells.df <- cbind(smooth.rna.df, smooth.atac.df)

  # Curate per-gene data frames for regression
  curate.data <- lapply(sele.genes, function(gene0) {
    peaks   <- names(which(gene2peak.map[gene0, ] == 1))
    data    <- subset(metacells.df, select = c(gene0, peaks))
    names(data)[1] <- "y"
    data$Subtype   <- as.factor(
      sapply(strsplit(rownames(data), "[._]+"), `[`, 2)
    )
    data
  })
  names(curate.data) <- sele.genes

  stime <- Sys.time()
  regr_list <- lapply(sort(sele.genes), function(x) {
    coef(summary(lm(y ~ Subtype + . - 1, data = curate.data[[x]])))
  })
  names(regr_list) <- sort(sele.genes)
  cat("Regression time:", format(Sys.time() - stime), "\n")

  # Summarise and filter regression results
  regr.sum <- rbindlist(lapply(names(regr_list), function(t) {
    x <- data.table(regr_list[[t]][, c(1, 4), drop = FALSE],
                    keep.rownames = TRUE)
    x$gene_name <- t
    x
  }))
  names(regr.sum)[c(1, 3)] <- c("peak_name", "P_value")
  regr.sum[, p_val_bonf := pmin(1, P_value * .N)]
  regr.sum$fdr <- p.adjust(regr.sum$P_value, method = "fdr")
  regr.sum$peak_name <- gsub("`", "", regr.sum$peak_name, fixed = TRUE)

  predicted_links <- regr.sum[P_value < 0.01 & abs(Estimate) > 0.05 &
                                grepl("^chr", peak_name)] %>%
    setorder(-Estimate)

  fwrite(predicted_links,
         file = paste0("filtered_gene_peak_links_metacell_", current.date, ".tsv"),
         sep  = "\t")
  cat("EP links:", nrow(predicted_links), "| unique genes:",
      length(unique(predicted_links$gene_name)), "\n")


  # ===========================================================================
  # STEP 2c: Promoter-enhancer annotation (optional; removes PPIs)
  # ===========================================================================

  annotation <- readRDS(path.hg38.ann)
  tss_ann    <- Signac:::CollapseToLongestTranscript(annotation)
  tss_ann    <- as.data.frame(resize(tss_ann, width = 1, fix = "start"))
  names(tss_ann)[1] <- "chr"
  tss_ann$Tss <- tss_ann$start
  gene_ann    <- tss_ann <- data.table(tss_ann)

  gene_ann_sele <- gene_ann[gene_name %in% predicted_links$gene_name]
  setkey(gene_ann_sele, gene_name)
  predicted_links[, chr          := gene_ann_sele[J(predicted_links$gene_name)]$chr]
  predicted_links[, Tss          := gene_ann_sele[J(predicted_links$gene_name)]$Tss]
  predicted_links[, promoter_pos := paste(chr, Tss - 1000, Tss + 1000, sep = "-")]

  peak.ann   <- annPeak2Gene(unique(predicted_links$peak_name), tss_ann, 1000)
  setkey(peak.ann, peak_name)
  peaks.prom <- peak.ann[nchar(gene_name) > 0]$peak_name
  peaks.nprom <- peak.ann[nchar(gene_name) == 0]$peak_name

  predicted_links.ep <- predicted_links[peak_name %in% peaks.nprom]

  # Assign nearest promoter peak to each gene
  gene_list <- unique(predicted_links.ep[, .(gene_name, chr, Tss)])
  gene_list[, start := Tss - 1000][, end := Tss + 1000]
  gene_list2peak <- geneOverlapPeak(gene_list,
                                     peak_names = predicted_links$peak_name,
                                     mid_dist   = 1000)
  setkey(gene_list2peak, gene_name)
  predicted_links.ep[, promoter_peak := gene_list2peak[J(predicted_links.ep$gene_name)]$peak_name]
  predicted_links.ep <- predicted_links.ep[promoter_peak != "Not_Found"]

  predicted_links.ep[, start   := as.integer(tstrsplit(peak_name, "-")[[2]]), by = peak_name]
  predicted_links.ep[, end     := as.integer(tstrsplit(peak_name, "-")[[3]]), by = peak_name]
  predicted_links.ep[, ep_dist := abs(Tss - (start + end) / 2)]

  fwrite(predicted_links.ep,
         file = paste0("gene_peak_enhancerpeaksonly_links_metacell_", current.date, ".tsv"),
         sep  = "\t")


  # ===========================================================================
  # STEP 2d: Differential peak accessibility (DAPs)
  # ===========================================================================

  cat("Finding DAPs...\n")
  DefaultAssay(multi.ctype.obj) <- "ATAC"
  Idents(multi.ctype.obj)       <- grouping.ident

  daps <- FindMarkers(
    object              = multi.ctype.obj,
    test.use            = "LR",
    ident.1             = "Myeloid-like",
    ident.2             = "B-lineage",
    max.cells.per.ident = 10000,
    logfc.threshold     = 0.25,
    min.pct             = 0.05,
    latent.vars         = c("nCount_ATAC", vars.to.regress)
  ) %>%
    filter(p_val_adj < 0.05) %>%
    rownames_to_column("peak") %>%
    mutate(cluster = if_else(avg_log2FC > 0, "Myeloid-like", "B-lineage"),
           avg_log2FC = abs(avg_log2FC))

  saveRDS(daps, paste0("DAPeaks_filtered_", current.date, ".rds"))


  # ===========================================================================
  # STEP 2e: chromVAR + differential TF motif activity
  # ===========================================================================

  cat("Running chromVAR...\n")
  DefaultAssay(multi.ctype.obj) <- "ATAC"
  multi.ctype.obj <- RunChromVAR(multi.ctype.obj,
                                  genome = BSgenome.Hsapiens.UCSC.hg38)

  DefaultAssay(multi.ctype.obj) <- "chromvar"
  multi.ctype.obj <- ScaleData(multi.ctype.obj)

  # Average z-score and d-score per cell state
  compute_avg_scores <- function(obj, assay_slot) {
    mtx <- if (assay_slot == "scale.data")
      obj@assays$chromvar@scale.data
    else
      obj@assays$chromvar@data
    avg <- sapply(unique(obj$Blast_Lineage), function(x)
      rowMeans(mtx[, obj$Blast_Lineage == x]))
    avg.df <- as.data.frame(avg)
    avg.df$gene <- sapply(strsplit(rownames(avg.df), "-"), `[`, 1)
    tf.conv  <- bitr(avg.df$gene, fromType = "ENSEMBL",
                      toType = "SYMBOL", OrgDb = org.Hs.eg.db, drop = FALSE)
    avg.df   <- left_join(avg.df, tf.conv, by = c("gene" = "ENSEMBL")) %>%
      filter(!is.na(SYMBOL)) %>%
      dplyr::select(-gene) %>%
      group_by(SYMBOL) %>%
      summarize(across(everything(), mean)) %>%
      tibble::column_to_rownames("SYMBOL")
    avg.df
  }

  avg.zscore.state <- compute_avg_scores(multi.ctype.obj, "scale.data")
  avg.dscore.state <- compute_avg_scores(multi.ctype.obj, "data")

  saveRDS(avg.zscore.state, paste0("avg_chromvar_zscore_",    current.date, ".rds"))
  saveRDS(avg.dscore.state, paste0("avg_chromvar_datascore_", current.date, ".rds"))

  # Differential TF activity (MAST; patient-regressed)
  Idents(multi.ctype.obj) <- grouping.ident
  diff.tf <- FindAllMarkers(
    object              = multi.ctype.obj,
    test.use            = "LR",
    max.cells.per.ident = 10000,
    latent.vars         = vars.to.regress,
    logfc.threshold     = 0.25,
    only.pos            = TRUE
  ) %>%
    mutate(gene = sapply(strsplit(rownames(.), "-"), `[`, 3)) %>%
    filter(p_val_adj < 0.05)

  saveRDS(diff.tf, paste0("diff_TFactivity_", current.date, ".rds"))

  # DEGs per blast lineage
  DefaultAssay(multi.ctype.obj) <- "RNA"
  Idents(multi.ctype.obj)       <- grouping.ident

  degs <- FindAllMarkers(
    object              = multi.ctype.obj,
    test.use            = "MAST",
    max.cells.per.ident = 10000,
    latent.vars         = vars.to.regress,
    logfc.threshold     = 0.25,
    only.pos            = TRUE
  )
  saveRDS(degs, paste0("degs_", current.date, ".RDS"))

  # TFs that are also differentially expressed genes
  tfs.deg <- degs[degs$gene %in% gsub(".1", "", diff.tf$gene) &
                    abs(degs$avg_log2FC) > 0.1 &
                    degs$p_val_adj < 0.05, ]
  saveRDS(tfs.deg, paste0("multi_tfs_deg_state_", current.date, ".rds"))


  # ===========================================================================
  # STEP 2f: Motif scan (TF × peak binary matrix)
  # ===========================================================================

  cat("Running motif scan...\n")
  metacell.atac <- readRDS(paste0(ctype_TC, "_atac_metacell_", current.date, ".rds"))
  peaks <- rownames(metacell.atac)

  pks <- data.table(tidyr::separate(
    data    = data.table(peak_name = peaks),
    col     = "peak_name",
    into    = c("chr", "start", "end"),
    remove  = FALSE
  ))
  pks$start <- as.integer(pks$start)
  pks$end   <- as.integer(pks$end)

  peaks.gr <- GenomicRanges::GRanges(
    seqnames = pks$chr,
    ranges   = IRanges::IRanges(start = pks$start, end = pks$end)
  )

  motif_ix <- matchMotifs(human_pwms_v2, peaks.gr,
                           genome    = BSgenome.Hsapiens.UCSC.hg38,
                           p.cutoff  = motif.pv)

  motif_pk_match <- motif_ix@assays@data$motifMatches
  rownames(motif_pk_match) <- paste0(
    motif_ix@rowRanges@seqnames, "-", motif_ix@rowRanges@ranges)
  colnames(motif_pk_match) <- motif_ix@colData$name

  saveRDS(motif_pk_match,
          paste0("motif_peak_match_", motif.pv, "_", current.date, ".rds"))

  # Per-cell expression and peak accessibility summaries
  mtx <- multi.ctype.obj@assays$RNA@counts
  expr.frac.state <- sapply(unique(multi.ctype.obj$Blast_Lineage), function(x)
    rowMeans(mtx[, multi.ctype.obj$Blast_Lineage == x] > 0))
  expr.median.state <- sapply(unique(multi.ctype.obj$Blast_Lineage), function(x)
    sparseMatrixStats::rowMedians(
      multi.ctype.obj@assays$RNA@data[, multi.ctype.obj$Blast_Lineage == x]))
  rownames(expr.median.state) <- rownames(mtx)
  expr.avg.state <- sapply(unique(multi.ctype.obj$Blast_Lineage), function(x)
    rowMeans(multi.ctype.obj@assays$RNA@data[, multi.ctype.obj$Blast_Lineage == x]))

  mtx.atac <- multi.ctype.obj@assays$ATAC@counts
  open.frac.pks <- sapply(unique(multi.ctype.obj$Blast_Lineage), function(x)
    rowMeans(mtx.atac[, multi.ctype.obj$Blast_Lineage == x] > 0))
  pm <- matrixStats::rowMaxs(open.frac.pks)
  weak.pks <- rownames(mtx.atac)[pm < peak.frac]

  saveRDS(expr.frac.state,   paste0("expr_frac_state_",   current.date, ".rds"))
  saveRDS(expr.median.state, paste0("median_expr_state_", current.date, ".rds"))
  saveRDS(expr.avg.state,    paste0("avg_expr_state_",    current.date, ".rds"))
  saveRDS(open.frac.pks,     paste0("peak_openFrac_byState_", current.date, ".rds"))
  saveRDS(weak.pks,          paste0("weak_peaks_",        current.date, ".rds"))


  # ===========================================================================
  # STEP 3: TRN assembly per cell state
  # ===========================================================================

  cat("Assembling TRNs...\n")

  # Load all intermediate files for network construction
  diff.tfs        <- readRDS(paste0("diff_TFactivity_",      current.date, ".rds")) %>%
    filter(abs(avg_log2FC) > 0.5)
  daps            <- readRDS(paste0("DAPeaks_filtered_",      current.date, ".rds"))
  degs            <- readRDS(paste0("degs_",                  current.date, ".RDS"))
  motif_pk_match  <- as.matrix(readRDS(paste0("motif_peak_match_", motif.pv, "_", current.date, ".rds")))
  open.frac.pks   <- readRDS(paste0("peak_openFrac_byState_", current.date, ".rds"))
  expr.frac.state <- readRDS(paste0("expr_frac_state_",       current.date, ".rds"))
  expr.median.state <- readRDS(paste0("median_expr_state_",   current.date, ".rds"))
  expr.avg.state  <- readRDS(paste0("avg_expr_state_",        current.date, ".rds"))
  tfs.deg         <- readRDS(paste0("multi_tfs_deg_state_",   current.date, ".rds"))
  weak.pks        <- readRDS(paste0("weak_peaks_",            current.date, ".rds"))
  avg.zscore.state <- readRDS(paste0("avg_chromvar_zscore_",  current.date, ".rds"))
  avg.dscore.state <- readRDS(paste0("avg_chromvar_datascore_", current.date, ".rds"))

  predicted_links <- if (useEnhancerPeakOnly) {
    fread(paste0("gene_peak_enhancerpeaksonly_links_metacell_", current.date, ".tsv"))
  } else {
    fread(paste0("filtered_gene_peak_links_metacell_",          current.date, ".tsv"))
  }

  filtered_links <- NULL
  tfs.list       <- list()

  for (cell_state0 in colnames(expr.frac.state)) {

    # -- Select target genes -----------------------------------------------
    deg.filter <- if (useAllDEGs) {
      if (includeNegativeDEGs)
        degs[degs$cluster == cell_state0 & abs(degs$avg_log2FC) > 0.25 &
               degs$p_val_adj < 0.05, ]
      else
        degs[degs$cluster == cell_state0 & degs$avg_log2FC  > 0.25 &
               degs$p_val_adj < 0.05 & degs$pct.1 > 0.05, ]
    } else {
      degs[degs$cluster == cell_state0 & degs$p_val_adj < 0.01 &
             degs$avg_log2FC > 0.5 & degs$pct.1 > 0.05, ]
    }
    interested_genes <- unique(deg.filter$gene)

    # -- Select TFs --------------------------------------------------------
    avg.zscore.state0 <- setNames(avg.zscore.state[[cell_state0]],
                                   rownames(avg.zscore.state))

    if (useTopDscoreTF) {
      avg.dscore.state0 <- setNames(avg.dscore.state[[cell_state0]],
                                     rownames(avg.dscore.state))
      interested_tfs <- names(sort(avg.dscore.state0[avg.dscore.state0 > 0],
                                    decreasing = TRUE))
    } else {
      interested_tfs <- unique(c(
        diff.tfs[diff.tfs$avg_log2FC > 0.5 & diff.tfs$cluster == cell_state0, ]$gene,
        tfs.deg[tfs.deg$cluster == cell_state0, ]$gene
      ))
    }

    tfs.list[[length(tfs.list) + 1]] <- interested_tfs

    # -- Select peaks -------------------------------------------------------
    open.frac.pks0 <- open.frac.pks[, cell_state0]
    filtered.peaks <- if (useDAP) {
      rownames(daps[daps$cluster == cell_state0, ])
    } else {
      names(which(open.frac.pks0 > peakFracCutoff))
    }

    # -- Filter TFs by expression -------------------------------------------
    tf.expr.frac <- expr.frac.state[interested_tfs, cell_state0]
    tf.expr.frac <- tf.expr.frac[!is.na(tf.expr.frac)]
    interested_tfs <- names(which(tf.expr.frac > tfexprfracCutoff))

    if (length(interested_tfs) <= 5) {
      interested_tfs <- intersect(
        names(sort(tf.expr.frac, decreasing = TRUE))[1:10],
        names(which(tf.expr.frac > 0.05))
      )
      message("  < 5 TFs for ", cell_state0, " — using top 10 by expression")
    }

    # -- Intersect peaks, TFs, and predicted links --------------------------
    predicted_links0 <- predicted_links[predicted_links$gene_name %in% interested_genes, ]
    interested_tfs   <- interested_tfs[interested_tfs %in% colnames(motif_pk_match)]
    sele.tf.mat      <- motif_pk_match[, interested_tfs, drop = FALSE]
    sele.peaks       <- intersect(names(which(rowSums(sele.tf.mat > 0) > 0)),
                                   filtered.peaks)
    predicted_links0 <- predicted_links0[predicted_links0$peak_name %in% sele.peaks, ]

    if (filterWeakPeak) {
      predicted_links0 <- predicted_links0[peak_name %notin% weak.pks]
    }

    predicted_links0 <- predicted_links0[Estimate > 0.2 & fdr < 0.05]
    predicted_links0$cell_state <- cell_state0
    filtered_links <- rbind(filtered_links, predicted_links0)

    # -- Build per-TF regulons ----------------------------------------------
    expr.avg.state0 <- expr.avg.state[, cell_state0]
    tf_regulons     <- list()

    for (TF0 in interested_tfs) {
      peaks0 <- names(which(sele.tf.mat[, TF0] > 0))
      ep0    <- predicted_links0[peak_name %in% peaks0]
      if (nrow(ep0) == 0) next

      ep0 <- ep0 %>%
        dplyr::select(gene_name, Estimate, fdr, peak_name) %>%
        mutate(peak_strength = open.frac.pks0[peak_name],
               score         = -log10(fdr),
               TF            = TF0) %>%
        group_by(gene_name) %>%
        summarize(TF            = TF0,
                  edge_strength  = sum(Estimate * peak_strength),
                  Estimate       = sum(Estimate),
                  score          = mean(score),
                  peaks          = paste(peak_name,     collapse = ","),
                  peak_strengths = paste(peak_strength, collapse = ","),
                  npeaks         = n(),
                  .groups        = "drop") %>%
        mutate(TF_strength   = avg.zscore.state0[TF0],
               gene_strength = expr.avg.state0[gene_name])

      tf_regulons[[TF0]] <- ep0[!duplicated(ep0), ]
    }

    tf_regulons.comb <- rbindlist(tf_regulons) %>% setorder(edge_strength)
    message(cell_state0, ": ", nrow(tf_regulons.comb), " regulons")

    # Build node table
    nodes1 <- tf_regulons.comb %>%
      dplyr::select(gene = TF, node_strength = TF_strength) %>%
      mutate(node_strength = node_strength / max(abs(node_strength)),
             node_type = "TF")
    nodes2 <- tf_regulons.comb %>%
      dplyr::select(gene = gene_name, node_strength = gene_strength) %>%
      filter(!gene %in% nodes1$gene) %>%
      mutate(node_strength = node_strength / max(abs(node_strength)),
             node_type = "gene")
    nodes <- unique(rbind(nodes1, nodes2))

    # -- Construct output filenames -----------------------------------------
    cell_state1 <- gsub("/", "-", gsub(" ", "-", cell_state0, fixed = TRUE),
                         fixed = TRUE)

    build_fname <- function(suffix, ext) {
      nm <- paste0(ctype_TC, ".EP_Prediction_TRN_useDEGs_", suffix, "_",
                   cell_state1, "_", ctype_TC,
                   if (useDAP) "_useDAPs" else paste0("_peakFracCutoff", peakFracCutoff),
                   "_ctypeMetacells", ext)
      if (!useEnhancerPeakOnly) nm <- sub(ext, paste0("_IncludePPI", ext), nm)
      if (useAllDEGs)           nm <- sub(ext, paste0("_useAllDEGs", ext), nm)
      if (includeNegativeDEGs)  nm <- sub(ext, paste0("_NegativeDEGs", ext), nm)
      nm
    }

    tf_regulon_fname <- build_fname("TF_regulons", ".tsv")
    nodes_fname      <- build_fname("nodes",       ".tsv")
    tf_rank_fname    <- build_fname("ConnectRank_TF_regulons", ".pdf")

    write.table(tf_regulons.comb, file = tf_regulon_fname,
                sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(nodes, file = nodes_fname,
                sep = "\t", row.names = FALSE, quote = FALSE)

    # -- Top-15 TF bar chart -----------------------------------------------
    pdata <- tf_regulons.comb %>%
      count(TF) %>%
      arrange(desc(n)) %>%
      slice_head(n = 15) %>%
      mutate(TF = factor(TF, levels = rev(TF)))

    p_rank <- ggplot(pdata, aes(x = TF, y = n, fill = n)) +
      geom_bar(stat = "identity") + coord_flip() +
      theme_classic() + scale_fill_viridis_c() + xlab("") +
      ggtitle(paste0(cell_state0, " — ", nrow(tf_regulons.comb), " targets"))

    ggsave(p_rank, filename = tf_rank_fname, device = "pdf",
           width = 4, height = 4)

    # -- Metascape input CSV (TF → gene list, up to 30 TFs) ----------------
    current.trn <- tf_regulons.comb %>%
      dplyr::select(TF, gene_name) %>%
      split(.$TF) %>%
      lapply(pull, gene_name) %>%
      .[order(sapply(., length), decreasing = TRUE)]

    max.len <- max(lengths(current.trn))
    current.trn <- lapply(current.trn, function(x) { length(x) <- max.len; x })
    current.trn <- bind_cols(current.trn)
    colnames(current.trn) <- paste0("Transcription Factor:", colnames(current.trn))
    current.trn[is.na(current.trn)] <- ""
    if (ncol(current.trn) > 30) current.trn <- current.trn[, 1:30]

    write.csv(current.trn, row.names = FALSE, na = "",
              file = gsub(".tsv", "_MetascapeInput.csv", tf_regulon_fname))
    write.csv(degs, row.names = FALSE, na = "", file = "DEGs.csv")

  }  # end cell state loop

  # Save filtered EP links
  ep_file <- build_fname("filtered_regrRes4_gene_peak_links_metacell", ".tsv")
  fwrite(filtered_links, file = ep_file, sep = "\t")

  cat("=== Completed TRN for stage:", ctype_TC, "===\n")
  cat(format(Sys.time()), "\n")
}


# =============================================================================
# Run pipeline for each stage
# =============================================================================
# To run in parallel: bplapply(ctypes.list, trn, BPPARAM = MulticoreParam())

ctypes.list <- c("IDX", "MRD")

for (stage in ctypes.list) {
  trn(stage)
}
