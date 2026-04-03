library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(tidyverse)
library(limma)
source('scDataAnalysis_Utilities.R')
source('dsAnalysis_utilities.R')
source('DoubletsFinder.R')
source('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/entropycalc.R')
library(gdata)
library(Matrix.utils)
library(clusterProfiler)
library(enrichplot)
`%notin%` = Negate(`%in%`)

# Setup --------------

iamp21.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/iamp21.rna.1.rds")
saveRDS(iamp21.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/iamp21.rna.1.rds")

# Data processing ----------
iamp21.rna <- NormalizeData(iamp21.rna)
iamp21.rna <- FindVariableFeatures(iamp21.rna, nfeatures = 4000)
iamp21.rna <- RunPCA(iamp21.rna, npc = 100, verbose = T, features = VariableFeatures(iamp21.rna))
iamp21.rna <- RunUMAP(iamp21.rna, dims = 1:50, reduction.name = "umapPCA50")
iamp21.rna <- FindNeighbors(iamp21.rna, dims = 1:50)
iamp21.rna <- FindClusters(iamp21.rna, res = 0.50)

FeaturePlot(iamp21.rna, features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), reduction="umapPCA100", ncol=3, label=T, repel=T, order=T)

# Annotation with marker genes  -------------
DimPlot(iamp21.rna, group.by="RNA_snn_res.2", reduction="umapPCA100", label=T, repel=T)

iamp21.rna@meta.data$Level1.anno <- as.character(iamp21.rna@meta.data$RNA_snn_res.1)
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$Level1.anno %notin% c(1,3,42,10,18,22,30,24)] <- 'Blast'
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$Level1.anno %in% c(1,3,42)] <- 'T/NK'
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$Level1.anno %in% c(10)] <- 'B'
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$Level1.anno %in% c(18,22,30)] <- 'Erythroid Progenitor'
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$Level1.anno %in% c(24)] <- 'Mono'
#iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$Level1.anno %in% c("54")] <- 'DC'

myColors <- colorRampPalette(brewer.pal(6,"Set1"))(6)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/iamp21RNA_87.pdf", height=10, width=10)
DimPlot(iamp21.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024), reduction = "umapPCA100")+NoLegend()
DimPlot(iamp21.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024), reduction = "umapPCA100")+NoLegend()
FeaturePlot(iamp21.rna, features = c("CD19"), order=T, label=F,  raster=T, raster.dpi = c(1024,1024), label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("MME"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("MS4A1"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("CD36"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3, reduction= "umapPCA100")
FeaturePlot(iamp21.rna, features = c("CD3D"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("LYZ"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3, reduction = "umapPCA100")
VlnPlot(iamp21.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/iamp21_AllMarkers.pdf",height=11.350, width=20)
FeaturePlot(iamp21.rna, features = c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3, ncol=3)
dev.off()

# Shannon Entropy Score -------------------

iamp21.rna <- FindNeighbors(iamp21.rna, dims = 1:100, reduction = 'pca')
iamp21.rna <- FindClusters(iamp21.rna, res = 0.25)
iamp21.rna <- FindClusters(iamp21.rna, res = 0.5)
iamp21.rna <- FindClusters(iamp21.rna, res = 0.75)
iamp21.rna <- FindClusters(iamp21.rna, res = 1)
iamp21.rna <- FindClusters(iamp21.rna, res = 2)
iamp21.rna = get.entropy(iamp21.rna,  c("RNA_snn_res.0.25","RNA_snn_res.0.5", "RNA_snn_res.0.75", "RNA_snn_res.1", "RNA_snn_res.2"),join.to.seurat = T)

FeaturePlot(iamp21.rna, features=c("average_entropy"), order=T, raster=F)
VlnPlot(iamp21.rna, features=c("average_entropy"), group.by = "Level1.anno", split.by = "seurat_clusters", ncol=10)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/ShannonEntropyFeaturePlots_iAMP21.pdf", width = 10, height = 10)
FeaturePlot(iamp21.rna, features = c("RNA_snn_res.0.25_entropy","RNA_snn_res.0.5_entropy", "RNA_snn_res.0.75_entropy", "RNA_snn_res.1_entropy", "RNA_snn_res.2_entropy", "average_entropy"), 
            raster=T, raster.dpi = c(1024,1024), combine=F, reduction="umapPCA100")
dev.off()

# Similarity to Healthy Reference -----------

seurat.rna.hd <- readRDS(file="/mnt/isilon/tan_lab/chenc6/MLLr_Project/HTAN_Data_Submission/scRNASeq_Level4/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
seurat.rna.hd <- NormalizeData(seurat.rna.hd)
seurat.rna.hd <- FindVariableFeatures(seurat.rna.hd)
seurat.rna.hd <- ScaleData(seurat.rna.hd)

iamp21.rna.sub <- iamp21.rna[,sample(colnames(iamp21.rna), size = 50000, replace=F)]
iamp21.rna.sub <- NormalizeData(iamp21.rna.sub)
iamp21.rna.sub <- FindVariableFeatures(iamp21.rna.sub)
iamp21.rna.sub <- ScaleData(iamp21.rna.sub)

# Find anchors
anchors <- FindTransferAnchors(reference = seurat.rna.hd, query = iamp21.rna.sub, reduction = 'cca')
saveRDS(anchors, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/iamp21.healthy.donor.anchors.rds")

anchors <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/iamp21.healthy.donor.anchors.rds")

# Transfer labels
predictions.Ctype <- TransferData(anchorset = anchors, refdata = seurat.rna.hd$Ctype, k.weight = 50, weight.reduction = 'cca')
iamp21.rna.sub <- AddMetaData(object = iamp21.rna.sub, metadata = predictions.Ctype)
saveRDS(iamp21.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/iamp21_downsample50K_withHealthySimilarity.rds')

iamp21.rna.sub <- readRDS(file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/iamp21_downsample50K_withHealthySimilarity.rds')

VlnPlot(iamp21.rna.sub, group.by = 'Level1.anno', features = c('prediction.score.max'))

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/scRNA_SimilarityScore_iAMP21.pdf", width = 7, height = 7)
iamp21.rna.sub@meta.data %>%
  ggplot( aes(x=Level1.anno, y=prediction.score.max, fill=Level1.anno)) +
  geom_boxplot(outlier.size=1) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  ggtitle("iAMP21 scRNA-Seq Similarity Score") +
  xlab("") + ylab("Similarity Score to healthy reference")
dev.off()

FeaturePlot(iamp21.rna.sub, features=c("prediction.score.Pre.B"), split.by = "Level1.anno", combine=F, order=T)
DimPlot(iamp21.rna.sub, group.by="predicted.id", split.by="Level1.anno", combine=F, ncol=3)

FeaturePlot(iamp21.rna.sub, features=c("prediction.score.max"), order=T, ncol=2) + 
  DimPlot(iamp21.rna.sub, group.by="predicted.id") + 
  DimPlot(iamp21.rna.sub, group.by="Stage") + 
  DimPlot(iamp21.rna.sub, group.by="patientID")

FeaturePlot(iamp21.rna.sub, features=c("prediction.score.max","MME"), order=T, ncol=2) 

# Percentage of blasts #####

blasts.metadata <- iamp21.rna@meta.data %>% mutate(NormalStatus = ifelse(Level1.anno == "Mono" | Level1.anno == "N/K" | Level1.anno == "Erythroid Progenitor" | Level1.anno == "T/NK" | Level1.anno == "B", "Normal",
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

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/blast_count_similarity_iamp21.pdf", width = 5, height = 5)
ggplot(blasts.metadata, aes(x=BlastLabeledPercent,y=countedLeukemia)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) + 
  labs(title="iAMP21", x="Annotated Blast Percentage", y="Counted Blast Percentage", caption=cor_coef) + 
  geom_text_repel(label=blasts.metadata$sampleID)
dev.off()

t.test(x=blasts.metadata$BlastLabeledPercent,y=blasts.metadata$countedLeukemia, two.tailed=TRUE, var.equal = TRUE)

# InferCNV -------------------------
library(infercnv)
iamp21.rna.downsample <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/InferCNV_Results/iAMP21/iAMP21_downsampledSeuratforInferCNV.rds")

DimPlot(iamp21.rna.downsample, group.by="Level1.anno", reduction="umapPCA100")

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/iamp21_infercnv_chr21.pdf", height=11.350, width=20)
FeaturePlot(iamp21.rna.downsample, features="proportion_scaled_cnv_chr21", ncol=2, reduction="umapPCA100") + 
  DimPlot(iamp21.rna.downsample, group.by="Level1.anno", reduction="umapPCA100", cols = myColors)
dev.off()

VlnPlot(iamp21.rna.downsample, features=c("proportion_scaled_dupli_chr6"), group.by="Level1.anno")

ggplot(iamp21.rna.downsample@meta.data) + geom_boxplot(aes(x=Level1.anno,y=proportion_scaled_dupli_chr6))

# Filtering out small clusters of normal cells ------------------------------------------------------------

PAVALZ_Bclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Cluster_Barcodes_Filtering/iamp21_Bcluster_1.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="B" & patientID=="PAVALZ") %>% rownames()
PATSBP_Tclust <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Cluster_Barcodes_Filtering/iamp21_Tcluster_1.csv", row.names=1) %>%
  dplyr::filter(Level1.anno=="T/NK" & patientID=="PATSBP") %>% rownames()

cells.PAVALZ <- iamp21.rna@meta.data %>% dplyr::filter(patientID=="PAVALZ") %>% rownames()
rna.PAVALZ <- subset(iamp21.rna, cells = cells.PAVALZ)
rna.PAVALZ <- NormalizeData(rna.PAVALZ)
rna.PAVALZ <- FindVariableFeatures(rna.PAVALZ, nfeatures = 4000)
rna.PAVALZ <- RunPCA(rna.PAVALZ, npc = 30, verbose = T, features = VariableFeatures(rna.PAVALZ))
rna.PAVALZ <- RunUMAP(rna.PAVALZ, dims = 1:30)
rna.PAVALZ <- FindNeighbors(rna.PAVALZ, dims = 1:30)
rna.PAVALZ <- FindClusters(rna.PAVALZ, res = 0.25)

cells.PATSBP <- iamp21.rna@meta.data %>% dplyr::filter(patientID=="PATSBP") %>% rownames()
rna.PATSBP <- subset(iamp21.rna, cells = cells.PATSBP)
rna.PATSBP <- NormalizeData(rna.PATSBP)
rna.PATSBP <- FindVariableFeatures(rna.PATSBP, nfeatures = 4000)
rna.PATSBP <- RunPCA(rna.PATSBP, npc = 30, verbose = T, features = VariableFeatures(rna.PATSBP))
rna.PATSBP <- RunUMAP(rna.PATSBP, dims = 1:30)
rna.PATSBP <- FindNeighbors(rna.PATSBP, dims = 1:30)
rna.PATSBP <- FindClusters(rna.PATSBP, res = 0.25)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/iamp21_TumorCellClusterAnnotation.pdf", height=11.350, width=20)
DimPlot(rna.PAVALZ, group.by="RNA_snn_res.0.25") + ggtitle("PAVALZ_Bclust") + 
  DimPlot(rna.PAVALZ, cells.highlight = PAVALZ_Bclust) +
  FeaturePlot(rna.PAVALZ, features=c("MS4A1"), order=T) + 
  FeaturePlot(rna.PAVALZ, features=c("MME"), order=T) + 
  DimPlot(rna.PAVALZ, group.by="Level1.anno")
DimPlot(rna.PATSBP, group.by="RNA_snn_res.0.25") + ggtitle("PATSBP_Tclust") + 
  DimPlot(rna.PATSBP, cells.highlight = PATSBP_Tclust) + 
  FeaturePlot(rna.PATSBP, features=c("CD3D"), order=T) + 
  FeaturePlot(rna.PATSBP, features=c("MME"), order=T)
dev.off()

# Relabeling small clusters of cells ------------------------------------------------------------
iamp21.rna@meta.data$Level1.anno[rownames(iamp21.rna@meta.data) %in% PAVALZ_Bclust] <- "Blast"
iamp21.rna@meta.data$Level1.anno[rownames(iamp21.rna@meta.data) %in% PATSBP_Tclust] <- "Blast"
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$RNA_snn_res.2 %in% c(37)] <- 'Mono'
iamp21.rna@meta.data$Level1.anno[iamp21.rna@meta.data$RNA_snn_res.2 %in% c(57)] <- 'B'


pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/iAMP21/iamp21RNA_87_relabeled.pdf", height=10, width=10)
DimPlot(iamp21.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction = "umapPCA100")
DimPlot(iamp21.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024), reduction = "umapPCA100")+NoLegend()
DimPlot(iamp21.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024), reduction = "umapPCA100")+NoLegend()
FeaturePlot(iamp21.rna, features = c("CD19"), order=T, label=F,  raster=T, raster.dpi = c(1024,1024), label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("MME"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("MS4A1"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("CD36"), order=T, label=F, raster=T, raster.dpi = c(1024,1024),label.size = 3, reduction= "umapPCA100")
FeaturePlot(iamp21.rna, features = c("CD3D"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3, reduction = "umapPCA100")
FeaturePlot(iamp21.rna, features = c("LYZ"), order=T, label=F, raster=T, raster.dpi = c(1024,1024), label.size = 3, reduction = "umapPCA100")
VlnPlot(iamp21.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

saveRDS(iamp21.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/iamp21.rna.1.rds")


