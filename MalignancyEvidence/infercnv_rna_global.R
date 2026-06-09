## Running inferCNV on global downsampled object 

## using all normal cell types as reference
## test on blast cells

setwd("/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence")

library(infercnv)
library(Seurat)
library(data.table)

#args = commandArgs(T)
#patientName = args[1]

seurat.rna = readRDS('global_downsample50K_withHealthySimilarity.rds')
meta <- fread('meta.csv')
rownames(meta) <- meta$patient_bc

seurat.rna <- AddMetaData(seurat.rna, meta)

# load annotation data
all_ann = seurat.rna@meta.data
mtx <- seurat.rna@assays$RNA@counts
barcodes <- colnames(mtx)
all_ann = subset(all_ann, select = c(Level1.anno, Subtype))

for(i in 1:nrow(all_ann)){
  if(all_ann$Level1.anno[i]=="Blast"){
    all_ann$Level1.anno[i] <- all_ann$Subtype[i]
  }
}

all_ann = subset(all_ann, select = Level1.anno)

names(all_ann) <- NULL


genelist = read.table('/mnt/isilon/tan_lab/chenc6/Tools/InferCNV/Hg38_Gene_Position.txt', 
                      sep = '\t', stringsAsFactors = F)

genelist = genelist[!duplicated(genelist$V1), ]

rownames(genelist) <- genelist$V1

#reorder chr
chrs = paste0('chr', 1:22)
genes_ann = NULL
for(chr0 in chrs) genes_ann = rbind(genes_ann, genelist[genelist$V2 == chr0, ])

genes_ann$V1 <- NULL

names(genes_ann) <- NULL


infercnv_obj = CreateInfercnvObject(raw_counts_matrix = as.matrix(mtx),
                                    annotations_file = all_ann,
                                    delim = "\t",
                                    gene_order_file = genes_ann,
                                    ref_group_names = c('B', 'Mono', 'Erythroid Progenitor', 'T/NK'))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             min_cells_per_gene = 5,
                             out_dir = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Global', 
                             cluster_by_groups = TRUE,
                             denoise = T,
                             output_format = 'pdf',
                             analysis_mode = 'subclusters', 
                             tumor_subcluster_partition_method = "qnorm",
                             HMM_report_by = "subcluster",
                             tumor_subcluster_pval = 0.05,
                             num_threads = 16,
                             HMM = TRUE)

saveRDS(infercnv_obj, file = 'Global_finalInferCNV.rds')

seurat.rna <- add_to_seurat(seurat_obj = seurat.rna, infercnv_output_path = "InferCNV_Results/Global")

saveRDS(seurat.rna, file = 'Global/Global_downsampledSeuratforInferCNV.rds')

## Plotting

infercnv_obj <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Global/Global_finalInferCNV.rds")

infercnv_obj <- sample_object(infercnv_obj, n_cells=500, on_references = F, on_observations = TRUE)

plot_cnv(
  infercnv_obj,
  out_dir = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Global",
  title = "inferCNV",
  obs_title = "Observations (Cells)",
  ref_title = "References (Cells)",
  cluster_by_groups = TRUE,
  cluster_references = TRUE,
  #k_obs_groups = 3,
  contig_cex = 1,
  x.center = mean(infercnv_obj@expr.data),
  x.range = "auto",
  hclust_method = "ward.D",
  color_safe_pal = FALSE,
  output_filename = "infercnv",
  output_format = "pdf",
  dynamic_resize = 1,
  ref_contig = NULL,
  write_expr_matrix = FALSE
)
