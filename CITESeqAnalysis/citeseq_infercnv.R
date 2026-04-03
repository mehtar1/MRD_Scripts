## using all normal cell types as reference
## test on blast cells

setwd("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq")

library(infercnv)
library(Seurat)
library(data.table)
library(tidyverse)
library(BiocParallel)
`%notin%` <- Negate(`%in%`) 

#args = commandArgs(T)
#patientName = args[1]

merged.obj <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/Round4_mergedFiltered_020226.rds")
merged.obj <- JoinLayers(merged.obj, assay="RNA", layers = c("counts", "data", "scale.data"))
merged.obj <- JoinLayers(merged.obj, assay="ADT")

## 02/03/26 Run without specifying blasts within each patient ------------------

obj.list <- SplitObject(merged.obj, split.by="patient.sample")

obj.list <- bplapply(obj.list, function(x){
  # JoinLayers(x, assay="RNA") %>% NormalizeData() %>% 
  FindVariableFeatures(x) %>% RunPCA(npcs = 40) %>% RunUMAP(dims=1:40)
})

sample.list <- unique(merged.obj$patient.sample)

# for (i in 1:length(sample.list)) {

param <- MulticoreParam(workers = 8)
options(scipen = 100)

options("Seurat.object.assay.version" = "v3")
  
#obj.list <- bplapply(obj.list, BPPARAM = param, function(x) {
  
# for(i in 2:length(obj.list)){ 

i <- commandArgs(trailingOnly = TRUE)
print(i)
  
x <- obj.list[[as.numeric(i)]]

current.sample <- x$patient.sample[1]
print(paste0("Current Sample is ", current.sample))

seurat.rna <- x

# load annotation data
all_ann = seurat.rna@meta.data
mtx <- GetAssayData(seurat.rna, assay="RNA", layer="counts")

# rm(seurat.rna)
barcodes <- colnames(mtx)
all_ann = subset(all_ann, select = c(Level1.anno))
names(all_ann) <- NULL

genelist = read.table('/mnt/isilon/tan_lab/chenc6/Tools/InferCNV/Hg38_Gene_Position.txt',
                      sep = '\t', stringsAsFactors = F)

genelist = genelist[!duplicated(genelist$V1), ]

rownames(genelist) <- genelist$V1

#reorder chr
chrs = paste0('chr', 1:22)
genes_ann = NULL
chrs <- c(chrs, "chrX")
for(chr0 in chrs) genes_ann = rbind(genes_ann, genelist[genelist$V2 == chr0, ])

genes_ann$V1 <- NULL

names(genes_ann) <- NULL

# saveRDS(seurat.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/iAMP21/iAMP21_downsampledSeuratforInferCNV.rds')

infercnv_obj = CreateInfercnvObject(raw_counts_matrix = as.matrix(mtx),
                                    annotations_file = all_ann,
                                    delim = "\t",
                                    gene_order_file = genes_ann,
                                    ref_group_names = c("B","T_NK"),
                                    chr_exclude = c("chrY", "chrM"))

dir.create(paste0('/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/',current.sample))

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir = paste0('/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/',current.sample),
                             cluster_by_groups = TRUE,
                             denoise = T,
                             output_format = 'pdf',
                             analysis_mode = 'subclusters',
                             #tumor_subcluster_partition_method = "qnorm",
                             #HMM_report_by = "subcluster",
                             #tumor_subcluster_pval = 0.05,
                             num_threads = 16,
                             no_plot=T,
                             no_prelim_plot = T,
                             HMM = TRUE)

saveRDS(infercnv_obj, file = paste0('/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/',current.sample,'_finalInferCNV_noRefGroup.rds'))

seurat.rna <- NormalizeData(seurat.rna)
seurat.rna <- add_to_seurat(seurat_obj = seurat.rna, infercnv_output_path = paste0("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/",current.sample))

saveRDS(seurat.rna, file = paste0('/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/',current.sample,'_InferCNV_noRefGroup.rds'))

## Re-plotting (03/25/24)

infercnv_obj <- readRDS(paste0("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/",current.sample,"_finalInferCNV_noRefGroup.rds"))

#infercnv_obj <- sample_object(infercnv_obj, n_cells=500, on_references = F, on_observations = TRUE)

plot_cnv(
  infercnv_obj,
  out_dir = paste0("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/",current.sample),
  title = current.sample,
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
  # }
#}) 

# paste0('/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Scripts/infercnv_results/',current.sample)




