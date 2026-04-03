library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(tidyverse)
library(limma)
library(dplyr)
source('scDataAnalysis_Utilities.R')
source('dsAnalysis_utilities.R')
source('DoubletsFinder.R')
source('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/entropycalc.R')
library(gdata)
library(Matrix.utils)
library(clusterProfiler)
library(enrichplot)
`%notin%` = Negate(`%in%`)

# Setup --------------

phlike.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike.rna.1.rds")
saveRDS(phlike.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike.rna.1.rds")

# Data processing ----------
phlike.rna <- NormalizeData(phlike.rna)
phlike.rna <- FindVariableFeatures(phlike.rna, nfeatures = 4000)
phlike.rna <- RunPCA(phlike.rna, npc = 100, verbose = T, features = VariableFeatures(phlike.rna))
phlike.rna <- RunUMAP(phlike.rna, dims = 1:30)
phlike.rna <- FindNeighbors(phlike.rna, dims = 1:30)
phlike.rna <- FindClusters(phlike.rna, res = 0.4)

# Annotation with marker genes  -------------
#At 0.25 resolution

FeaturePlot(phlike.rna, features=c("nCount_RNA"), order=T, raster=F, max.cutoff = "q95")
phlike.rna@meta.data$Level1.anno <- as.character(phlike.rna@meta.data$RNA_snn_res.0.25)
phlike.rna@meta.data$Level1.anno[phlike.rna@meta.data$Level1.anno %notin% c(0,13,15,17,21,22)] <- 'Blast'
phlike.rna@meta.data$Level1.anno[phlike.rna@meta.data$Level1.anno %in% c(0,13)] <- 'T/NK'
phlike.rna@meta.data$Level1.anno[phlike.rna@meta.data$Level1.anno %in% c(15)] <- 'B'
phlike.rna@meta.data$Level1.anno[phlike.rna@meta.data$Level1.anno %in% c(17,21)] <- 'Erythroid Progenitor'
phlike.rna@meta.data$Level1.anno[phlike.rna@meta.data$Level1.anno %in% c(22)] <- 'Myeloid'

phlike.rna@meta.data %>% filter(Level1.anno != "Blast") %>% ggplot(aes(Level1.anno, RNA_snn_res.2, fill = RNA_snn_res.2)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = RNA_snn_res.2), 
            position = position_dodge(width=1),
            vjust=1.5)

Idents(phlike.rna) <- "Level1.anno"
DimPlot(phlike.rna, group.by="RNA_snn_res.2", raster=F, label=T, repel=T) + NoLegend() + 
  DimPlot(phlike.rna, cells=WhichCells(phlike.rna, idents = "T/NK"), group.by="RNA_snn_res.2", label=T, repel=T)

#DoHeatmap(phlike.rna, cells=WhichCells(phlike.rna, idents="T/NK"), group.by = "RNA_snn_res.2", features=c("CD3D", "GNLY"))

myColors <- colorRampPalette(brewer.pal(6,"Set1"))(6)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/PhlikeRNA_87.pdf", height=10, width=10)
DimPlot(phlike.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024))+NoLegend()
DimPlot(phlike.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024))+NoLegend()
FeaturePlot(phlike.rna, features = c("CD19"), order=T, label=F,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(phlike.rna, features = c("MME"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(phlike.rna, features = c("MS4A1"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(phlike.rna, features = c("CD36"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(phlike.rna, features = c("CD3D"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(phlike.rna, features = c("LYZ"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3)
VlnPlot(phlike.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/PhlikeRNA_AllMarkers.pdf",height=11.350, width=20)
FeaturePlot(phlike.rna, features = c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3, ncol=3)
dev.off()

# Shannon Entropy Score -------------------

phlike.rna <- FindNeighbors(phlike.rna, dims = 1:100, reduction = 'pca')
phlike.rna <- FindClusters(phlike.rna, res = 0.25)
phlike.rna <- FindClusters(phlike.rna, res = 0.5)
phlike.rna <- FindClusters(phlike.rna, res = 0.75)
phlike.rna <- FindClusters(phlike.rna, res = 1)
phlike.rna <- FindClusters(phlike.rna, res = 2)
phlike.rna = get.entropy(phlike.rna,  c("RNA_snn_res.0.25","RNA_snn_res.0.5", "RNA_snn_res.0.75", "RNA_snn_res.1", "RNA_snn_res.2"),join.to.seurat = T)

FeaturePlot(phlike.rna, features=c("average_entropy"), order=T, raster=T)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/ShannonEntropyFeaturePlots_PhLike.pdf", width = 10, height = 10)
FeaturePlot(phlike.rna, features = c("RNA_snn_res.0.25_entropy","RNA_snn_res.0.5_entropy", "RNA_snn_res.0.75_entropy", "RNA_snn_res.1_entropy", "RNA_snn_res.2_entropy", "average_entropy"), 
            raster=T, raster.dpi = c(1024,1024), combine=F)
dev.off()

# Similarity to Healthy Reference -----------

seurat.rna.hd <- readRDS(file="/mnt/isilon/tan_lab/chenc6/MLLr_Project/HTAN_Data_Submission/scRNASeq_Level4/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
seurat.rna.hd <- NormalizeData(seurat.rna.hd)
seurat.rna.hd <- FindVariableFeatures(seurat.rna.hd)
seurat.rna.hd <- ScaleData(seurat.rna.hd)

phlike.rna.sub <- phlike.rna[,sample(colnames(seurat.rna), size = 50000, replace=F)]
phlike.rna.sub <- NormalizeData(phlike.rna.sub)
phlike.rna.sub <- FindVariableFeatures(phlike.rna.sub)
phlike.rna.sub <- ScaleData(phlike.rna.sub)

# Find anchors
anchors <- FindTransferAnchors(reference = seurat.rna.hd, query = phlike.rna.sub, reduction = 'cca')
saveRDS(anchors, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike.healthy.donor.anchors.rds")

# Transfer labels
predictions.Ctype <- TransferData(anchorset = anchors, refdata = seurat.rna.hd$Ctype, k.weight = 50, weight.reduction = 'cca')
phlike.rna.sub <- AddMetaData(object = phlike.rna.sub, metadata = predictions.Ctype)
saveRDS(phlike.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike_downsample50K_withHealthySimilarity.rds')
phlike.rna.sub <- readRDS(file='/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike_downsample50K_withHealthySimilarity.rds')

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/scRNA_SimilarityScore_phlike.pdf", width = 7, height = 7)
phlike.rna.sub@meta.data %>%
  ggplot( aes(x=Level1.anno, y=prediction.score.max, fill=Level1.anno)) +
  geom_boxplot(outlier.size=1) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  ggtitle("Ph-Like scRNA-Seq Similarity Score") +
  xlab("") + ylab("Similarity Score to healthy reference")
dev.off()

# Idents(phlike.rna.sub) <- "Level1.anno"
# DimPlot(phlike.rna.sub, group.by="RNA_snn_res.2", raster=F, label=T, repel=T) + NoLegend() + 
#   DimPlot(phlike.rna.sub, cells=WhichCells(phlike.rna.sub, idents = "T/NK"), group.by="prediction.score.T", label=T, repel=T)

VlnPlot(phlike.rna.sub, group.by = 'Level1.anno', features = c('prediction.score.max'))

FeaturePlot(phlike.rna.sub, features=c("prediction.score.max"), split.by = "Level1.anno", combine=F)
DimPlot(phlike.rna.sub, group.by="predicted.id", split.by="Level1.anno", combine=F, ncol=3)

FeaturePlot(phlike.rna.sub, features=c("prediction.score.max"), order=T, ncol=2) + 
  DimPlot(phlike.rna.sub, group.by="predicted.id") + 
  DimPlot(phlike.rna.sub, group.by="Stage") + 
  DimPlot(phlike.rna.sub, group.by="patientID")

FeaturePlot(phlike.rna.sub, features=c("prediction.score.max","MME"), order=T, ncol=2) 
  
# Percentage of blasts ---------------------------

blasts.metadata <- phlike.rna@meta.data %>% mutate(NormalStatus = ifelse(Level1.anno == "Myeloid" | Level1.anno == "N/K" | Level1.anno == "Erythroid Progenitor" | Level1.anno == "T/NK" | Level1.anno == "B", "Normal",
                                                                   (ifelse(Level1.anno == "Blast", "Blast", NA))))
blasts.metadata <- blasts.metadata %>% 
  group_by(patientID, NormalStatus, Stage) %>% 
  summarize(count=n()) %>% 
  spread(key=NormalStatus, value=count) %>% 
  mutate(BlastLabeledPercent = Blast/(Blast+Normal))

blasts.percentage <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/LeukemicBlastPercentage.csv")

blasts.metadata <- blasts.metadata %>% left_join(blasts.percentage, by=c("patientID"="USI", "Stage")) 
blasts.metadata <- rename(blasts.metadata, "countedLeukemia"="leukemia..")
blasts.metadata <- blasts.metadata %>% mutate(dif_in_percents = countedLeukemia-BlastLabeledPercent) %>% 
  select(!c("Subtype"))
blasts.metadata <- blasts.metadata %>% unite(sampleID, c(patientID,Stage), remove=T)

cor_coef <- paste0("Pearson correlation: r=",round(cor(blasts.metadata$BlastLabeledPercent, blasts.metadata$countedLeukemia, method = "pearson"), digits=3))

require("ggrepel")

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/blast_count_similarity_phlike.pdf", width = 5, height = 5)
ggplot(blasts.metadata, aes(x=BlastLabeledPercent,y=countedLeukemia)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) + 
  labs(title="Ph-Like", x="Annotated Blast Percentage", y="Counted Blast Percentage", caption=cor_coef) + 
  geom_text_repel(label=blasts.metadata$sampleID)
dev.off()
  
t.test(x=blasts.metadata$BlastLabeledPercent,y=blasts.metadata$countedLeukemia, two.tailed=TRUE, var.equal = TRUE)

# InferCNV -------------------------
library(infercnv)
phlike.infercnv <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Ph-Like/Ph-Like_finalInferCNV.rds")
phlike.rna.downsample <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Ph-Like/Ph-Like_downsampledSeuratforInferCNV.rds")

phlike.rna.downsample <- NormalizeData(phlike.rna.downsample)
DimPlot(phlike.rna.downsample, group.by="Level1.anno")

phlike.rna.downsample <- add_to_seurat(seurat_obj = phlike.rna.downsample, infercnv_output_path = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Ph-Like")
DimPlot(phlike.rna.downsample, group.by="top_dupli_4")

#chromosome 13 proves that some of these small cell clusters might be real 
Idents(phlike.rna.downsample) <- "patientID"
FeaturePlot(phlike.rna.downsample, features="proportion_scaled_dupli_chr13", ncol=2) + DimPlot(phlike.rna.downsample, group.by="Level1.anno")
FeaturePlot(phlike.rna.downsample, features="proportion_scaled_cnv_chr19", ncol=2, label=T, repel=T, label.size = 2) + DimPlot(phlike.rna.downsample, group.by="Level1.anno")
FeaturePlot(phlike.rna.downsample, features="proportion_scaled_cnv_chr22", ncol=2, order=T) + DimPlot(phlike.rna.downsample, group.by="Level1.anno")

VlnPlot(phlike.rna.downsample, features=c("proportion_scaled_dupli_chr6"), group.by="Level1.anno")

ggplot(phlike.rna.downsample@meta.data) + geom_boxplot(aes(x=Level1.anno,y=proportion_scaled_dupli_chr6))

# Filtering out small clusters of normal cells ------------------------------------------------------------

PAYJEG_Tclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Cluster_Barcodes_Filtering/phlike_Tcluster_3.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="T/NK" & patientID=="PAYJEG") %>% rownames()
PAYMKH_Tclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Cluster_Barcodes_Filtering/phlike_Tcluster_1.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="T/NK" & patientID=="PAYMKH") %>% rownames()
PAXYRK_Tclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Cluster_Barcodes_Filtering/phlike_Tcluster_2.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="T/NK" & patientID=="PAXYRK") %>% rownames()
PAXYRK_Bclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Cluster_Barcodes_Filtering/phlike_Bcluster_1.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="B" & patientID=="PAXYRK") %>% rownames()
PAYJEG_Bclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Cluster_Barcodes_Filtering/phlike_Bcluster_2.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="B" & patientID=="PAYJEG") %>% rownames()

cells.PAYJEG <- phlike.rna@meta.data %>% dplyr::filter(patientID=="PAYJEG") %>% rownames()
rna.PAYJEG <- subset(phlike.rna, cells = cells.PAYJEG)
rna.PAYJEG <- NormalizeData(rna.PAYJEG)
rna.PAYJEG <- FindVariableFeatures(rna.PAYJEG, nfeatures = 4000)
rna.PAYJEG <- RunPCA(rna.PAYJEG, npc = 30, verbose = T, features = VariableFeatures(rna.PAYJEG))
rna.PAYJEG <- RunUMAP(rna.PAYJEG, dims = 1:30)
rna.PAYJEG <- FindNeighbors(rna.PAYJEG, dims = 1:30)
rna.PAYJEG <- FindClusters(rna.PAYJEG, res = 0.25)

cells.PAYMKH <- phlike.rna@meta.data %>% dplyr::filter(patientID=="PAYMKH") %>% rownames()
rna.PAYMKH <- subset(phlike.rna, cells = cells.PAYMKH)
rna.PAYMKH <- NormalizeData(rna.PAYMKH)
rna.PAYMKH <- FindVariableFeatures(rna.PAYMKH, nfeatures = 4000)
rna.PAYMKH <- RunPCA(rna.PAYMKH, npc = 30, verbose = T, features = VariableFeatures(rna.PAYMKH))
rna.PAYMKH <- RunUMAP(rna.PAYMKH, dims = 1:30)
rna.PAYMKH <- FindNeighbors(rna.PAYMKH, dims = 1:30)
rna.PAYMKH <- FindClusters(rna.PAYMKH, res = 0.25)

cells.PAXYRK <- phlike.rna@meta.data %>% dplyr::filter(patientID=="PAXYRK") %>% rownames()
rna.PAXYRK <- subset(phlike.rna, cells = cells.PAXYRK)
rna.PAXYRK <- NormalizeData(rna.PAXYRK)
rna.PAXYRK <- FindVariableFeatures(rna.PAXYRK, nfeatures = 4000)
rna.PAXYRK <- RunPCA(rna.PAXYRK, npc = 30, verbose = T, features = VariableFeatures(rna.PAXYRK))
rna.PAXYRK <- RunUMAP(rna.PAXYRK, dims = 1:30)
rna.PAXYRK <- FindNeighbors(rna.PAXYRK, dims = 1:30)
rna.PAXYRK <- FindClusters(rna.PAXYRK, res = 0.25)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/ph-like_TumorCellClusterAnnotation.pdf", height=11.350, width=20)
DimPlot(rna.PAYJEG, group.by="RNA_snn_res.0.25") + 
  (DimPlot(rna.PAYJEG, cells.highlight = PAYJEG_Bclust) + ggtitle("PAYJEG_Bclust")) +
  (DimPlot(rna.PAYJEG, cells.highlight = PAYJEG_Tclust) + ggtitle("PAYJEG_Tclust")) +
  FeaturePlot(rna.PAYJEG, features=c("MME"), order=T) + 
  FeaturePlot(rna.PAYJEG, features=c("CD3D"), order=T) + 
  FeaturePlot(rna.PAYJEG, features=c("MS4A1"), order=T) 

DimPlot(rna.PAYMKH, group.by="RNA_snn_res.0.25") + 
  (DimPlot(rna.PAYMKH, cells.highlight = PAYMKH_Tclust) + ggtitle("PAYMKH_Tclust")) +
  FeaturePlot(rna.PAYMKH, features=c("CD3D"), order=T) + 
  FeaturePlot(rna.PAYMKH, features=c("MME"), order=T)

DimPlot(rna.PAXYRK, group.by="RNA_snn_res.0.25") + 
  (DimPlot(rna.PAXYRK, cells.highlight = PAXYRK_Bclust) + ggtitle("PAXYRK_Bclust")) +
  (DimPlot(rna.PAXYRK, cells.highlight = PAXYRK_Tclust) + ggtitle("PAXYRK_Tclust")) +
  FeaturePlot(rna.PAXYRK, features=c("MME"), order=T) + 
  FeaturePlot(rna.PAXYRK, features=c("CD3D"), order=T) + 
  FeaturePlot(rna.PAXYRK, features=c("MS4A1"), order=T)
dev.off()

#PAYJEG_B - non-malignant 
#PAYJEG_T - non-malignant 
#PAYMKH_T - malignant 
#PAXYRK_B - malignant 
#PAXYRK_T - malignant 

# Relabeling small clusters of cells ------------------------------------------------------------
phlike.rna@meta.data$Level1.anno[rownames(phlike.rna@meta.data) %in% PAYMKH_Tclust] <- 'Blast'
phlike.rna@meta.data$Level1.anno[rownames(phlike.rna@meta.data) %in% PAXYRK_Bclust] <- "Blast"
phlike.rna@meta.data$Level1.anno[rownames(phlike.rna@meta.data) %in% PAXYRK_Tclust] <- "Blast"

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Ph-Like/PhlikeRNA_87_relabeled.pdf", height=10, width=10)
DimPlot(phlike.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024))
DimPlot(phlike.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024))+NoLegend()
DimPlot(phlike.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024))+NoLegend()
FeaturePlot(phlike.rna, features = c("CD19"), order=T, label=F,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(phlike.rna, features = c("MME"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(phlike.rna, features = c("MS4A1"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(phlike.rna, features = c("CD36"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(phlike.rna, features = c("CD3D"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(phlike.rna, features = c("LYZ"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3)
VlnPlot(phlike.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

saveRDS(phlike.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike.rna.1.rds")

