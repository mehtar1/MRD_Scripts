# =============================================================================
# CITE-Seq inferCNV — Per-Sample Copy Number Inference
# =============================================================================
#
# Description:
#   Runs inferCNV on a single CITE-Seq patient sample to identify somatic
#   copy-number alterations in leukemic blast cells, using normal B and T/NK
#   cells from the same sample as reference. After inferCNV completes, the
#   CNV scores are added back into the Seurat object via add_to_seurat() and
#   the CNV heatmap is plotted.
#
#   Designed for HPC array submission: the sample index is passed as a
#   command-line argument, so one job per sample can run in parallel.
#
# Usage (interactive):
#   Set SAMPLE_INDEX in USER CONFIGURATION below and source the script.
#
# Usage (HPC array job, e.g. SLURM):
#   Rscript citeseq_infercnv.R <sample_index>
#   where <sample_index> is a 1-based integer into the sample list.
#
#   Example SLURM array header:
#     #SBATCH --array=1-10
#     Rscript citeseq_infercnv.R $SLURM_ARRAY_TASK_ID
#
# Inputs (set in USER CONFIGURATION):
#   - Merged filtered Seurat RDS (all samples, post-QC)
#   - Hg38 gene position file (tab-separated: gene, chr, start, end)
#
# Outputs (written to dir.infercnv/<sample>/):
#   - inferCNV result RDS
#   - inferCNV heatmap PDF
#   - Seurat RDS with CNV scores added (via add_to_seurat)
#
# Dependencies:
#   infercnv, Seurat, data.table, tidyverse, BiocParallel
#
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────

library(infercnv)
library(Seurat)
library(data.table)
library(tidyverse)
library(BiocParallel)

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────

# Merged Seurat object (all samples combined, with Level1.anno and
# patient.sample in metadata)
path.merged.rds  <- "/path/to/Objects/merged_filtered.rds"

# Hg38 gene position file (tab-separated, no header):
#   col1 = gene name, col2 = chromosome, col3 = start, col4 = end
path.gene.order  <- "/path/to/Tools/InferCNV/Hg38_Gene_Position.txt"

# Output directories
dir.infercnv     <- "/path/to/infercnv_results/"   # per-sample subdirs created here
dir.objects      <- "/path/to/Objects/"

# Reference cell type labels used as normal baseline in inferCNV
# (must match Level1.anno values in your object)
REF_GROUPS <- c("B", "T_NK")

# Chromosomes to exclude from CNV analysis
CHR_EXCLUDE <- c("chrY", "chrM")

# inferCNV run parameters
INFERCNV_CUTOFF     <- 0.1    # 0.1 for 10x Genomics; 1.0 for Smart-seq2
INFERCNV_THREADS    <- 16
INFERCNV_ANALYSIS_MODE <- "subclusters"

# UMI normalization parameters for per-sample preprocessing
N_PCA_DIMS <- 40

# Sample index: set here for interactive use, or pass as command-line argument
# for HPC array jobs (see Usage above). 1-based index into sample list.
SAMPLE_INDEX <- 1   # overridden by command-line arg if provided

options(scipen = 100)
options("Seurat.object.assay.version" = "v3")


# ── Resolve sample index from command line (HPC mode) ────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  SAMPLE_INDEX <- as.integer(args[1])
  cat("Running in HPC mode. Sample index:", SAMPLE_INDEX, "\n")
}


# =============================================================================
# SECTION 1: Load data and select sample
# =============================================================================

cat("Loading merged object...\n")
merged.obj <- readRDS(path.merged.rds)
merged.obj <- JoinLayers(merged.obj, assay = "RNA",
                          layers = c("counts", "data", "scale.data"))
merged.obj <- JoinLayers(merged.obj, assay = "ADT")

# Split by patient sample and run basic dimensionality reduction per sample
# (used for UMAP visualization of each sample independently)
obj.list     <- SplitObject(merged.obj, split.by = "patient.sample")
sample.names <- names(obj.list)

if (SAMPLE_INDEX < 1 || SAMPLE_INDEX > length(sample.names)) {
  stop("SAMPLE_INDEX ", SAMPLE_INDEX,
       " is out of range (1 to ", length(sample.names), ")")
}

obj.list <- bplapply(obj.list, function(x) {
  FindVariableFeatures(x) %>%
    RunPCA(npcs = N_PCA_DIMS) %>%
    RunUMAP(dims = 1:N_PCA_DIMS)
})

current.sample <- sample.names[SAMPLE_INDEX]
seurat.rna     <- obj.list[[SAMPLE_INDEX]]
cat("Processing sample:", current.sample, "\n")


# =============================================================================
# SECTION 2: Prepare inferCNV inputs
# =============================================================================
# Builds the gene order matrix required by inferCNV, sorting genes by
# chromosome order (chr1–22, chrX) and removing duplicates.

cat("Preparing gene order file...\n")
gene.list <- read.table(path.gene.order, sep = "\t",
                         stringsAsFactors = FALSE)
gene.list <- gene.list[!duplicated(gene.list$V1), ]
rownames(gene.list) <- gene.list$V1

chr.order <- c(paste0("chr", 1:22), "chrX")
genes.ann <- do.call(rbind, lapply(chr.order,
                                    function(chr) gene.list[gene.list$V2 == chr, ]))
genes.ann$V1 <- NULL
names(genes.ann) <- NULL

# Extract count matrix and cell annotations (Level1.anno only)
mtx      <- GetAssayData(seurat.rna, assay = "RNA", layer = "counts")
cell.ann <- seurat.rna@meta.data["Level1.anno"]
names(cell.ann) <- NULL


# =============================================================================
# SECTION 3: Run inferCNV
# =============================================================================

out.dir <- file.path(dir.infercnv, current.sample)
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

cat("Creating inferCNV object...\n")
infercnv.obj <- CreateInfercnvObject(
  raw_counts_matrix = as.matrix(mtx),
  annotations_file  = cell.ann,
  delim             = "\t",
  gene_order_file   = genes.ann,
  ref_group_names   = REF_GROUPS,
  chr_exclude       = CHR_EXCLUDE
)

cat("Running inferCNV (this may take several hours)...\n")
infercnv.obj <- infercnv::run(
  infercnv.obj,
  cutoff               = INFERCNV_CUTOFF,
  out_dir              = out.dir,
  cluster_by_groups    = TRUE,
  denoise              = TRUE,
  output_format        = "pdf",
  analysis_mode        = INFERCNV_ANALYSIS_MODE,
  num_threads          = INFERCNV_THREADS,
  no_plot              = TRUE,     # heatmap plotted separately below
  no_prelim_plot       = TRUE,
  HMM                  = TRUE
)

# Save inferCNV result object
path.infercnv.rds <- file.path(dir.infercnv,
                                paste0(current.sample, "_finalInferCNV.rds"))
saveRDS(infercnv.obj, file = path.infercnv.rds)
cat("inferCNV result saved to:", path.infercnv.rds, "\n")


# =============================================================================
# SECTION 4: Add CNV scores to Seurat object and save
# =============================================================================
# add_to_seurat() reads the inferCNV output directory and appends per-cell
# CNV proportion columns (has_dupli_chrN, proportion_dupli_chrN, etc.)
# to the Seurat object metadata.

seurat.rna <- NormalizeData(seurat.rna)
seurat.rna <- add_to_seurat(seurat_obj        = seurat.rna,
                             infercnv_output_path = out.dir)

path.seurat.out <- file.path(dir.objects,
                              paste0(current.sample, "_InferCNV.rds"))
saveRDS(seurat.rna, file = path.seurat.out)
cat("Seurat object with CNV scores saved to:", path.seurat.out, "\n")


# =============================================================================
# SECTION 5: Plot inferCNV heatmap
# =============================================================================
# Generates the standard inferCNV heatmap PDF. Run this after inferCNV
# completes; it can also be re-run independently by loading the saved RDS.

infercnv.obj <- readRDS(path.infercnv.rds)

plot_cnv(
  infercnv.obj,
  out_dir               = out.dir,
  title                 = current.sample,
  obs_title             = "Observations (Cells)",
  ref_title             = "References (Cells)",
  cluster_by_groups     = TRUE,
  cluster_references    = TRUE,
  contig_cex            = 1,
  x.center              = mean(infercnv.obj@expr.data),
  x.range               = "auto",
  hclust_method         = "ward.D",
  color_safe_pal        = FALSE,
  output_filename       = "infercnv",
  output_format         = "pdf",
  dynamic_resize        = 1,
  ref_contig            = NULL,
  write_expr_matrix     = FALSE
)

cat("inferCNV heatmap saved to:", out.dir, "\n")
cat("Sample", current.sample, "complete.\n")
