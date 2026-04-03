library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(tidyverse)
library(limma)
#source('scDataAnalysis_Utilities.R')
#source('dsAnalysis_utilities.R')
#source('DoubletsFinder.R')
#source('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MalignancyEvidence/Scripts/entropycalc.R')
library(gdata)
# library(Matrix.utils)
library(clusterProfiler)
library(enrichplot)
`%notin%` = Negate(`%in%`)

# Setup --------------

ph.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/ph.rna.1.rds")
saveRDS(ph.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/ph.rna.1.rds")

# Data processing ----------
ph.rna <- NormalizeData(ph.rna)
ph.rna <- FindVariableFeatures(ph.rna, nfeatures = 4000)
ph.rna <- RunPCA(ph.rna, npc = 100, verbose = T, features = VariableFeatures(ph.rna))
ph.rna <- RunUMAP(ph.rna, dims = 1:50, reduction.name="umapPCA50")
ph.rna <- FindNeighbors(ph.rna, dims = 1:50)
ph.rna <- FindClusters(ph.rna, res = 0.50)

DimPlot(ph.rna, group.by="patientID", reduction="umapPCA50", label=T, repel=T)
FeaturePlot(ph.rna, features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), reduction="umapPCA50", ncol=3, label=T, repel=T, order=T)
DotPlot(ph.rna, features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"))

ph.rna@meta.data$seurat_clusters <- ph.rna@meta.data$RNA_snn_res.0.5
Idents(ph.rna) <- "RNA_snn_res.2"
FeaturePlot(ph.rna, features=c("TPM2"), reduction="umapPCA50", label=T, repel=T, order=T, raster=F)
DimPlot(ph.rna, group.by="RNA_snn_res.2", label=T, repel = T)
VlnPlot(ph.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="RNA_snn_res.0.5", pt.size=0, flip=T)

markers.37 <- FindMarkers(ph.rna, ident.1 = 37)
markers.37 %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% arrange(desc(avg_log2FC))
#37 looks like DC, so put in with Monocytes 

#Probable B - 17
#T - 13,7
#Mono - 16
#EP - 14

# Annotation with marker genes  -------------
ph.rna@meta.data$Level1.anno <- as.character(ph.rna@meta.data$RNA_snn_res.0.5)
ph.rna@meta.data$Level1.anno[ph.rna@meta.data$Level1.anno %notin% c(17,13,7,16,14)] <- 'Blast'
ph.rna@meta.data$Level1.anno[ph.rna@meta.data$Level1.anno %in% c(13,7)] <- 'T/NK'
ph.rna@meta.data$Level1.anno[ph.rna@meta.data$Level1.anno %in% c(17)] <- 'B'
ph.rna@meta.data$Level1.anno[ph.rna@meta.data$Level1.anno %in% c(14)] <- 'Erythroid Progenitor'
ph.rna@meta.data$Level1.anno[ph.rna@meta.data$Level1.anno %in% c(16)] <- 'Mono'
#ph.rna@meta.data$Level1.anno[ph.rna@meta.data$Level1.anno %in% c("54")] <- 'DC'

myColors <- colorRampPalette(brewer.pal(6,"Set1"))(6)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/phRNA_87.pdf", height=10, width=10)
DimPlot(ph.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, reduction="umapPCA50", raster.dpi = c(1024,1024))
DimPlot(ph.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA50")
DimPlot(ph.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA50")
DimPlot(ph.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA50")
DimPlot(ph.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024), reduction="umapPCA50")+NoLegend()
DimPlot(ph.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024), reduction="umapPCA50")+NoLegend()
FeaturePlot(ph.rna, features = c("CD19"), reduction="umapPCA50", order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(ph.rna, features = c("MME"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(ph.rna, features = c("MS4A1"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(ph.rna, features = c("CD36"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(ph.rna, features = c("CD3D"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(ph.rna, features = c("LYZ"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024), label.size = 3)
VlnPlot(ph.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/phRNA_AllMarkers.pdf",height=11.350, width=20)
FeaturePlot(ph.rna, features = c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3, ncol=3)
dev.off()

# Shannon Entropy Score -------------------

ph.rna <- FindNeighbors(ph.rna, dims = 1:100, reduction = 'pca')
ph.rna <- FindClusters(ph.rna, res = 0.25)
ph.rna <- FindClusters(ph.rna, res = 0.5)
ph.rna <- FindClusters(ph.rna, res = 0.75)
ph.rna <- FindClusters(ph.rna, res = 1)
ph.rna <- FindClusters(ph.rna, res = 2)
ph.rna = get.entropy(ph.rna,  c("RNA_snn_res.0.25","RNA_snn_res.0.5", "RNA_snn_res.0.75", "RNA_snn_res.1", "RNA_snn_res.2"),join.to.seurat = T)

ph.rna$seurat_clusters <- ph.rna$RNA_snn_res.0.5

FeaturePlot(ph.rna, features=c("average_entropy"), order=T, raster=F)
#VlnPlot(ph.rna, features=c("average_entropy"), group.by = "Level1.anno", split.by = "seurat_clusters", ncol=10)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/ShannonEntropyFeaturePlots_ph.pdf", width = 10, height = 10)
FeaturePlot(ph.rna, features = c("RNA_snn_res.0.25_entropy","RNA_snn_res.0.5_entropy", "RNA_snn_res.0.75_entropy", "RNA_snn_res.1_entropy", "RNA_snn_res.2_entropy", "average_entropy"), 
            raster=T, raster.dpi = c(1024,1024), combine=F)
dev.off()

# Similarity to Healthy Reference -----------

seurat.rna.hd <- readRDS(file="/mnt/isilon/tan_lab/chenc6/MLLr_Project/HTAN_Data_Submission/scRNASeq_Level4/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
seurat.rna.hd <- NormalizeData(seurat.rna.hd)
seurat.rna.hd <- FindVariableFeatures(seurat.rna.hd)
seurat.rna.hd <- ScaleData(seurat.rna.hd)

ph.rna.sub <- ph.rna[,sample(colnames(ph.rna), size = 50000, replace=F)]
ph.rna.sub <- NormalizeData(ph.rna.sub)
ph.rna.sub <- FindVariableFeatures(ph.rna.sub)
ph.rna.sub <- ScaleData(ph.rna.sub)

# Find anchors
anchors <- FindTransferAnchors(reference = seurat.rna.hd, query = ph.rna.sub, reduction = 'cca')
saveRDS(anchors, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/ph.healthy.donor.anchors.rds")

anchors <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/ph.healthy.donor.anchors.rds")

# Transfer labels
predictions.Ctype <- TransferData(anchorset = anchors, refdata = seurat.rna.hd$Ctype, k.weight = 50, weight.reduction = 'cca')
ph.rna.sub <- AddMetaData(object = ph.rna.sub, metadata = predictions.Ctype)
saveRDS(ph.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/ph_downsample50K_withHealthySimilarity.rds')

ph.rna.sub <- readRDS(file='/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/ph_downsample50K_withHealthySimilarity.rds')

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/scRNA_SimilarityScore_phpositive.pdf", width = 7, height = 7)
ph.rna.sub@meta.data %>%
  ggplot( aes(x=Level1.anno, y=prediction.score.max, fill=Level1.anno)) +
  geom_boxplot(outlier.size=1) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  ggtitle("Ph-Positive scRNA-Seq Similarity Score") +
  xlab("") + ylab("Similarity Score to healthy reference")
dev.off()

VlnPlot(ph.rna.sub, group.by = 'Level1.anno', features = c('prediction.score.max'))

FeaturePlot(ph.rna.sub, features=c("prediction.score.MEP"), split.by = "Level1.anno",order=T, ncol=3)
DimPlot(ph.rna.sub, group.by="predicted.id", split.by="Level1.anno", combine=F, ncol=3, label=T, repel=T)

FeaturePlot(ph.rna.sub, features=c("prediction.score.max"), order=T, ncol=2) + 
  DimPlot(ph.rna.sub, group.by="predicted.id") + 
  DimPlot(ph.rna.sub, group.by="Stage") + 
  DimPlot(ph.rna.sub, group.by="patientID")

FeaturePlot(ph.rna.sub, features=c("prediction.score.max","MME"), order=T, ncol=2) 

# Percentage of blasts---------------------------------------

blasts.metadata <- ph.rna@meta.data %>% mutate(NormalStatus = ifelse(Level1.anno == "Mono" | Level1.anno == "N/K" | Level1.anno == "Erythroid Progenitor" | Level1.anno == "T/NK" | Level1.anno == "B", "Normal",
                                                                         (ifelse(Level1.anno == "Blast", "Blast", NA))))
blasts.metadata <- blasts.metadata %>% 
  group_by(patientID, NormalStatus, Stage) %>% 
  summarize(count=n()) %>% 
  spread(key=NormalStatus, value=count) %>% 
  mutate(BlastLabeledPercent = Blast/(Blast+Normal))

blasts.percentage <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/LeukemicBlastPercentage.csv")

blasts.metadata <- blasts.metadata %>% left_join(blasts.percentage, by=c("patientID"="USI", "Stage")) 
blasts.metadata <- rename(blasts.metadata, "countedLeukemia"="leukemia..")
blasts.metadata <- blasts.metadata %>% mutate(dif_in_percents = countedLeukemia-BlastLabeledPercent) %>% 
  select(!c("Subtype"))
blasts.metadata <- blasts.metadata %>% unite(sampleID, c(patientID,Stage), remove=T)

cor_coef <- paste0("Pearson correlation: r=",round(cor(blasts.metadata$BlastLabeledPercent, blasts.metadata$countedLeukemia, method = "pearson"), digits=3))

require("ggrepel")

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/blast_count_similarity_ph.pdf", width = 5, height = 5)
ggplot(blasts.metadata, aes(x=BlastLabeledPercent,y=countedLeukemia)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) + 
  labs(title="Ph-Positive", x="Annotated Blast Percentage", y="Counted Blast Percentage", caption=cor_coef) + 
  geom_text_repel(label=blasts.metadata$sampleID)
dev.off()


t.test(x=blasts.metadata$BlastLabeledPercent,y=blasts.metadata$countedLeukemia, two.tailed=TRUE, var.equal = TRUE)

# InferCNV -------------------------
library(infercnv)
ph.rna.downsample <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/Ph-Positive/Ph-Positive_downsampledSeuratforInferCNV.rds")
DimPlot(ph.rna.downsample, group.by="Level1.anno")
DimPlot(ph.rna.downsample, group.by="top_dupli_4")

FeaturePlot(ph.rna.downsample, features="proportion_scaled_cnv_chr3", ncol=2, order=T) + DimPlot(ph.rna.downsample, group.by="Level1.anno")

FeaturePlot(ph.rna.downsample, features="proportion_scaled_cnv_chr1", ncol=2, order=T) + DimPlot(ph.rna.downsample, group.by="sample.name")

VlnPlot(ph.rna.downsample, features=c("proportion_scaled_dupli_chr6"), group.by="Level1.anno")

ggplot(ph.rna.downsample@meta.data) + geom_boxplot(aes(x=Level1.anno,y=proportion_scaled_dupli_chr6))

# Filtering out small clusters of normal cells ------------------------------------------------------------

PAREZH_Bclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Cluster_Barcodes_Filtering/ph_Bcluster_1.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="B" & patientID=="PAREZH") %>% rownames()
PARCAH_Tclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Cluster_Barcodes_Filtering/ph_Tcluster_1.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="T/NK" & patientID=="PARCAH") %>% rownames()

cells.PAREZH <- ph.rna@meta.data %>% dplyr::filter(patientID=="PAREZH") %>% rownames()
rna.PAREZH <- subset(ph.rna, cells = cells.PAREZH)
rna.PAREZH <- NormalizeData(rna.PAREZH)
rna.PAREZH <- FindVariableFeatures(rna.PAREZH, nfeatures = 4000)
rna.PAREZH <- RunPCA(rna.PAREZH, npc = 30, verbose = T, features = VariableFeatures(rna.PAREZH))
rna.PAREZH <- RunUMAP(rna.PAREZH, dims = 1:30)
rna.PAREZH <- FindNeighbors(rna.PAREZH, dims = 1:30)
rna.PAREZH <- FindClusters(rna.PAREZH, res = 0.25)

cells.PARCAH <- ph.rna@meta.data %>% dplyr::filter(patientID=="PARCAH") %>% rownames()
rna.PARCAH <- subset(ph.rna, cells = cells.PARCAH)
rna.PARCAH <- NormalizeData(rna.PARCAH)
rna.PARCAH <- FindVariableFeatures(rna.PARCAH, nfeatures = 4000)
rna.PARCAH <- RunPCA(rna.PARCAH, npc = 30, verbose = T, features = VariableFeatures(rna.PARCAH))
rna.PARCAH <- RunUMAP(rna.PARCAH, dims = 1:30)
rna.PARCAH <- FindNeighbors(rna.PARCAH, dims = 1:30)
rna.PARCAH <- FindClusters(rna.PARCAH, res = 0.25)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/ph_TumorCellClusterAnnotation.pdf", height=11.350, width=20)
DimPlot(rna.PAREZH, group.by="RNA_snn_res.0.25") + ggtitle("PAREZH_Bclust") +
  DimPlot(rna.PAREZH, cells.highlight = PAREZH_Bclust) +
  DimPlot(rna.PAREZH, group.by="Level1.anno") +
  FeaturePlot(rna.PAREZH, features=c("MME"), order=T)
DimPlot(rna.PARCAH, group.by="RNA_snn_res.0.25") + ggtitle("PARCAH_Tclust") +
  DimPlot(rna.PARCAH, cells.highlight = PARCAH_Tclust) + 
  DimPlot(rna.PARCAH, group.by="Level1.anno") +
  FeaturePlot(rna.PARCAH, features=c("MME"), order=T) + 
  FeaturePlot(rna.PARCAH, features=c("CD3D"), order=T)
dev.off()

# Relabeling small clusters of cells ------------------------------------------------------------
ph.rna@meta.data$Level1.anno[rownames(ph.rna@meta.data) %in% PAREZH_Bclust] <- "Blast"
ph.rna@meta.data$Level1.anno[ph.rna@meta.data$RNA_snn_res.2 %in% c(23,37,47,45)] <- 'Mono'

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/Ph-Positive/phRNA_87_relabeled.pdf", height=10, width=10)
DimPlot(ph.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, reduction="umapPCA50", raster.dpi = c(1024,1024))
DimPlot(ph.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA50")
DimPlot(ph.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA50")
DimPlot(ph.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA50")
DimPlot(ph.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024), reduction="umapPCA50")+NoLegend()
DimPlot(ph.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024), reduction="umapPCA50")+NoLegend()
FeaturePlot(ph.rna, features = c("CD19"), reduction="umapPCA50", order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(ph.rna, features = c("MME"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(ph.rna, features = c("MS4A1"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(ph.rna, features = c("CD36"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(ph.rna, features = c("CD3D"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(ph.rna, features = c("LYZ"), reduction="umapPCA50",order=T, label=T, raster=T, raster.dpi = c(1024,1024), label.size = 3)
VlnPlot(ph.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

saveRDS(ph.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/ph.rna.1.rds")

