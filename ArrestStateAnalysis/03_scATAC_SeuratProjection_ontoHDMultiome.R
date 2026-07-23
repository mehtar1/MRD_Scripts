# =============================================================================
# Projection of patient scATAC-seq onto the healthy-donor multiome reference
# =============================================================================
#
# Description:
#   For each patient sample, this script builds a Signac ChromatinAssay from the
#   scATAC-pro peak-by-cell matrix (re-quantified on the reference peak set),
#   runs TF-IDF + LSI, and projects the cells onto the healthy-donor (HD)
#   multiome reference by anchor-based label transfer (FindTransferAnchors /
#   MapQuery) in the reference LSI space. Transferred cell-type labels and
#   B-trajectory pseudotime are recorded, the query is embedded in the reference
#   UMAP, and a per-sample projection UMAP is exported.
#
# Inputs:
#   - HD multiome reference Seurat object (.rds) with an "ATAC" assay,
#     an "lsi" reduction, and a "UMAP_lsi25" reduction model.
#   - Per-sample scATAC-pro re-quantified matrices (matrix.rds) on the
#     reference peak set.
#   - Annotation table (CSV) mapping barcodes to samples and prior labels.
#   - hg38 gene annotation (GRanges .rds).
#
# Outputs:
#   - Per-sample projected Seurat objects (.rds).
#   - Per-sample projection UMAP figures (.eps).
#   - A combined cell-to-reference mapping table (.tsv).
#
# Dependencies:
#   Seurat, Signac, ggplot2, data.table, RColorBrewer
# =============================================================================


# ── Libraries ────────────────────────────────────────────────────────────────
library(Seurat)
library(Signac)
library(ggplot2)
library(data.table)
library(RColorBrewer)

`%notin%` <- Negate(`%in%`)


# ── USER CONFIGURATION ───────────────────────────────────────────────────────
# Edit these paths to match your environment. All downstream paths are derived
# from the project root below.

dir.project <- "/path/to/project"                       # project root
dir.scATACpro <- file.path(dir.project, "scATACpro")    # per-sample scATAC-pro output
dir.figures   <- file.path(dir.project, "Figures/Projection_HDMultiome/SeuratMethod_ATAConUMAPlsi25")
dir.objects   <- file.path(dir.project, "SeuratObj/Projection_SeuratMethod_ATACtoUMAPlsi25")
dir.metadata  <- file.path(dir.project, "metadata")

# HD multiome reference object (must contain the "ATAC" assay, "lsi" reduction,
# and "UMAP_lsi25" reduction model)
path.reference <- file.path(dir.project,
                            "SeuratObj/HDMultiome/Downsample_Btraj",
                            "Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels.rds")

# hg38 gene annotation (GRanges)
path.annotation <- file.path(dir.project, "reference/GRanges_hg38.rds")

# Barcode-to-sample annotation table (prior labels for the query cells)
path.query.annotation <- file.path(dir.metadata,
                                    "HTAN_MRD_AllSamples_RNAproj2healthyCells_to_MultiomeRef.csv")

# Combined output mapping table
path.output.map <- file.path(dir.metadata,
                             "HTAN_MRD_AllSamples_ATACproj2healthyCells_to_MultiomeRef_UMAPlsi25.tsv")

# Reference reductions / dimensions
LSI_DIMS       <- 2:25
REF_REDUCTION  <- "lsi"
REF_UMAP_MODEL <- "UMAP_lsi25"

for (d in c(dir.figures, dir.objects, dir.metadata)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ── Load the healthy-donor multiome reference ────────────────────────────────
seurat.multi <- readRDS(path.reference)
DefaultAssay(seurat.multi) <- "ATAC"

# Reference embedding used as the grey background of every projection plot
HDUMAP <- data.frame(seurat.multi@reductions[[REF_UMAP_MODEL]]@cell.embeddings)
HDUMAP$Ctype <- "Healthy"
names(HDUMAP) <- c("refUMAP_1", "refUMAP_2", "Ctype")

annotation <- readRDS(path.annotation)
query.annotation <- fread(path.query.annotation)


# ── Sample list to project ───────────────────────────────────────────────────
# One entry per patient sample (<PatientID>_<Timepoint>). Edit as needed.
sampleID <- c(
  "PAPBCW_IDX", "PATJRM_IDX", "PATKFT_IDX", "PATSBP_IDX", "PAVALZ_MRD",
  "PAWIYT_IDX", "PAYLZC_IDX", "PAYNDU_IDX", "PATJRM_MRD", "PATKFT_MRD",
  "PAVALZ_IDX", "PAVRGL_IDX", "PAWIYT_MRD", "PAYLZC_MRD", "PAYNDU_MRD",
  "PARCAH_IDX", "PAREZH_IDX", "PARNLW_IDX", "PATSBP_MRD", "PAVEPS_MRD",
  "PAVTKL_IDX", "PAVXLK_IDX", "PAWARW_IDX", "PAWVAB_IDX", "PAXPPX_IDX",
  "PARCAH_MRD", "PAREZH_MRD", "PARNLW_MRD", "PAVEPS_IDX", "PAVRGL_MRD",
  "PAVTKL_MRD", "PAVXLK_MRD", "PAWARW_MRD", "PAWVAB_MRD", "PAXPPX_MRD",
  "2615_IDX", "2631_MRD", "3315_IDX", "3398_MRD", "3840_IDX", "3898_MRD",
  "3984_RLP", "4048_IDX", "4057_RLP", "4232_MRD", "4823_RLP",
  "PAYFXM_IDX", "PAPDNB_IDX", "PAPDNB_MRD", "PAPHWH_IDX", "PAPHWH_MRD",
  "PAXFFK_IDX", "PAXFFK_MRD", "PAYJEG_IDX", "PAYJEG_MRD",
  "PAYIPS_IDX", "PAYIPS_MRD", "PAPIWL_IDX", "PAPIWL_MRD", "PAXYMC_IDX",
  "PAXYMC_MRD", "PAYMKH_IDX", "PAYMKH_MRD", "PAYIYE_IDX", "PAYIYE_MRD",
  "PAYPDR_IDX", "PAYPDR_MRD", "PATGLR_IDX", "PATGLR_MRD", "PATGLR_RLP",
  "PAWEJR_IDX", "PAWEJR_MRD", "PAWEJR_RLP", "PAXEIP_IDX", "PAXEIP_MRD",
  "PAXEIP_RLP", "PAPDUV_IDX", "PAPDUV_MRD", "PAPDUV_RLP", "3928_IDX", "3970_MRD",
  "PAXYRK_IDX", "PAXYRK_MRD", "PAYBCZ_IDX", "PAYBCZ_MRD", "PAYCFJ_IDX", "PAYCFJ_MRD"
)


# ── Projection function ──────────────────────────────────────────────────────
# Projects one patient sample onto the HD reference and returns the query
# Seurat object with transferred labels, pseudotime, and reference UMAP coords.
project_sample <- function(sample0) {
  message("Processing ", sample0)

  # Load re-quantified peak matrix (reference peak set) and prefix barcodes
  mtx <- readRDS(file.path(dir.scATACpro, sample0,
                           "output/filtered_matrix/MACS2/FILTER",
                           "reConstruct_matrix_forProjection/matrix.rds"))
  colnames(mtx) <- paste0(sample0, "_", colnames(mtx))

  # Restrict to annotated barcodes for this sample
  meta.pt <- data.frame(query.annotation[query.annotation$Patient.Sample == sample0, ])
  rownames(meta.pt) <- meta.pt$patient_bc
  mtx <- mtx[, meta.pt$patient_bc]

  # Build query object and run TF-IDF + LSI
  chrom.pt <- CreateChromatinAssay(counts = mtx, sep = c(":", "-"),
                                   annotation = annotation)
  seurat.pt <- CreateSeuratObject(counts = chrom.pt, assay = "ATAC",
                                  meta.data = meta.pt)
  seurat.pt <- RunTFIDF(seurat.pt)
  seurat.pt <- FindTopFeatures(seurat.pt, min.cutoff = "q0")
  seurat.pt <- RunSVD(seurat.pt)
  seurat.pt <- RunUMAP(seurat.pt, reduction = "lsi", dims = LSI_DIMS)
  seurat.pt <- FindNeighbors(seurat.pt, reduction = "lsi", dims = LSI_DIMS)
  seurat.pt <- FindClusters(seurat.pt, verbose = FALSE, algorithm = 3)

  # Anchor-based transfer onto the reference LSI space
  current_anchors <- FindTransferAnchors(
    reference        = seurat.multi,
    reference.assay  = "ATAC",
    query            = seurat.pt,
    query.assay      = "ATAC",
    reference.reduction = REF_REDUCTION,
    dims             = LSI_DIMS
  )
  seurat.pt <- MapQuery(
    anchorset        = current_anchors,
    reference        = seurat.multi,
    query            = seurat.pt,
    refdata          = list(celltype.lsi   = "NewCellTypeUpdate",
                            pseudotime.lsi = "Btraj.PeseudoTime"),
    reference.reduction = REF_REDUCTION,
    reduction.model  = REF_UMAP_MODEL,
    projectumap.args = list(reduction.name = "refUMAPlsi25_ATAC",
                            neighbor.name  = "HDRefbc_UMAPlsi25",
                            return.neighbor = TRUE)
  )

  # Projection UMAP: HD reference (grey) + projected query cells (colored)
  refUMAP0 <- data.frame(seurat.pt@reductions$refUMAPlsi25_ATAC@cell.embeddings)
  refUMAP0 <- cbind(refUMAP0, seurat.pt@meta.data["Level1.anno"])
  names(refUMAP0) <- c("refUMAP_1", "refUMAP_2", "Ctype")

  UMAPforPlot <- rbind(HDUMAP, refUMAP0)
  myColors <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(refUMAP0$Ctype)) + 1)
  names(myColors) <- c(unique(refUMAP0$Ctype), "Healthy")
  myColors["Healthy"] <- "#cccccc"
  UMAPforPlot$Ctype <- factor(UMAPforPlot$Ctype,
                              levels = c("Healthy", sort(unique(refUMAP0$Ctype))))

  pp1 <- ggplot(UMAPforPlot, aes(x = refUMAP_1, y = refUMAP_2, col = Ctype)) +
    geom_point(size = 0.2) +
    scale_colour_manual(name = "annotation", values = myColors)
  ggsave(pp1, device = "eps", width = 9, height = 7,
         filename = file.path(dir.figures,
                              paste0("proj_direct_ATAC_", sample0,
                                     "_toMultiomeRef_UMAPlsi25_SeuratMethod.eps")))

  # Save projected object
  saveRDS(seurat.pt,
          file = file.path(dir.objects,
                           paste0("seurat_", sample0,
                                  "_proj2MultiomeRef_UMAPlsi25_SeuratMethod.rds")))

  seurat.pt@meta.data
}


# ── Run projection over all samples ──────────────────────────────────────────
cellmap <- data.frame()
for (sample0 in sampleID) {
  cellmap <- rbind(cellmap, project_sample(sample0))
}

write.table(cellmap, file = path.output.map, sep = "\t",
            row.names = FALSE, quote = FALSE)

message("Projection complete: ", length(sampleID), " samples.")