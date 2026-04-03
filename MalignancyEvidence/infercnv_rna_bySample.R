## using all normal cell types as reference
## test on blast cells

setwd("/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay")
library(infercnv)
library(Seurat)
library(data.table)

#args = commandArgs(T)
#patientName = args[1]

iamp21.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/iamp21.rna.1.rds")
ph.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/ph.rna.1.rds")
#kmt2ar.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/kmt2ar.rna.1.rds")
phlike.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike.rna.1.rds")

iamp21.list <- SplitObject(iamp21.rna, split.by = "patientID")
ph.list <- SplitObject(ph.rna, split.by = "patientID")
phlike.list <- SplitObject(phlike.rna, split.by = "patientID")

sample_names <- c("PAYJEG", "PAYMKH", "PAXYRK", "PAREZH", "PARCAH", "PAVALZ", "PATSBP")
# sample_subtypes <- c("Ph-like","Ph-like","Ph-like","Ph+", "Ph+", "iAMP21", "iAMP21")
sample_list <- list(phlike.list$PAYJEG, phlike.list$PAYMKH, phlike.list$PAXYRK, 
                    ph.list$PAREZH, ph.list$PARCAH, iamp21.list$PAVALZ, iamp21.list$PATSBP)

rm(phlike.rna)
rm(ph.rna)
rm(iamp21.rna)

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


for(sample in sample_list){
  annotation <- sample@meta.data
  annotation <- subset(annotation, select = Level1.anno)
  
  #annotating blast by sampleID
  for(i in 1:nrow(annotation)){
    #print(i)
    if(annotation$Level1.anno[i]=="Blast"){
      #ids <- str_split()
      annotation$Level1.anno[i] <- substr((rownames(annotation))[i], start=1, stop=(nchar(rownames(annotation)[i])-19))
    }
  }
  
  mtx <- sample@assays$RNA@counts
  
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix = as.matrix(mtx),
                                      annotations_file = annotation,
                                      delim = "\t",
                                      gene_order_file = genes_ann,
                                      ref_group_names = c('B', 'Mono', 'Erythroid Progenitor', 'T/NK'))
  
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               min_cells_per_gene = 5,
                               out_dir = paste('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Cluster_Filtering/',sample_names[sample]), 
                               cluster_by_groups = FALSE,
                               denoise = T,
                               output_format = 'pdf',
                               analysis_mode = 'subclusters', 
                               tumor_subcluster_partition_method = "qnorm",
                               HMM_report_by = "subcluster",
                               tumor_subcluster_pval = 0.05,
                               num_threads = 16,
                               HMM = TRUE)
  
  sample <- add_to_seurat(seurat_obj = sample, infercnv_output_path = paste('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Cluster_Filtering/',sample_names[sample]))
}

saveRDS(sample_list, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Cluster_Filtering/sample_seurat_list.rds")


#sample.type <- seurat.rna@meta.data[["Subtype"]][1]

# downsample seurat object 
# n <- as.integer(ncol(seurat.rna)/10)
# n = 10000
# seurat.rna.sub <- seurat.rna[, sample(colnames(seurat.rna), size=n, replace=F)]






#all_ann.sub$level.1.anno.viscello = ifelse(grepl(all_ann.sub$level.1.anno.viscello, pattern = 'Blast'),
#                                   'blast', all_ann.sub$level.1.anno.viscello)
#all_ann.sub$level.1.anno.viscello = ifelse(!grepl(all_ann.sub$level.1.anno.viscello, pattern = 'Blast'),
#                                   'normal', all_ann.sub$level.1.anno.viscello)
#mtx.sub = mtx.sub[, all_ann.sub$level.1.anno.viscello != 'NA']
#all_ann.sub = subset(all_ann.sub, level.1.anno.viscello != 'NA')


#names(all_ann.sub) <- NULL






#saveRDS(seurat.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/iAMP21/iAMP21_downsampledSeuratforInferCNV.rds')

#saveRDS(infercnv_obj, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/iAMP21/iAMP21_finalInferCNV.rds')

#seurat.rna.sub <- NormalizeData(seurat.rna.sub)
#seurat.rna.sub <- add_to_seurat(seurat_obj = seurat.rna.sub, infercnv_output_path = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/iAMP21")

#saveRDS(seurat.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/iAMP21/iAMP21_downsampledSeuratforInferCNV.rds')










# 
# infercnv_obj <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/Ph-Positive/run.final.infercnv_obj")
# 
# infercnv_obj@tumor_subclusters = NULL
# plot_cnv(
#   infercnv_obj,
#   out_dir = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/Ph-Positive",
#   title = "inferCNV",
#   obs_title = "Observations (Cells)",
#   ref_title = "References (Cells)",
#   cluster_by_groups = FALSE,
#   cluster_references = TRUE,
#   k_obs_groups = 3,
#   contig_cex = 1,
#   x.center = mean(infercnv_obj@expr.data),
#   x.range = "auto",
#   hclust_method = "ward.D",
#   color_safe_pal = FALSE,
#   output_filename = "infercnv",
#   output_format = "png",
#   png_res = 300,
#   dynamic_resize = 0,
#   ref_contig = NULL,
#   write_expr_matrix = FALSE
# )

