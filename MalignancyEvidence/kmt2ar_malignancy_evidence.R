library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(tidyverse)
library(limma)
# source('scDataAnalysis_Utilities.R')
# source('dsAnalysis_utilities.R')
# source('DoubletsFinder.R')
# source('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/entropycalc.R')
library(gdata)
# library(Matrix.utils)
library(clusterProfiler)
library(enrichplot)
library(infercnv)
`%notin%` = Negate(`%in%`)

# Setup --------------

kmt2ar.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/kmt2ar.rna.1.rds")
saveRDS(kmt2ar.rna, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/kmt2ar.rna.1.rds")

# Data processing ----------
kmt2ar.rna <- NormalizeData(kmt2ar.rna)
kmt2ar.rna <- FindVariableFeatures(kmt2ar.rna, nfeatures = 4000)
kmt2ar.rna <- RunPCA(kmt2ar.rna, npc = 100, verbose = T, features = VariableFeatures(kmt2ar.rna))
kmt2ar.rna <- RunUMAP(kmt2ar.rna, dims = 1:75, reduction.name="umapPCA75")
kmt2ar.rna <- FindNeighbors(kmt2ar.rna, dims = 1:75)
kmt2ar.rna <- FindClusters(kmt2ar.rna, res = 0.50)

DimPlot(kmt2ar.rna, group.by="RNA_snn_res.2", reduction="umapPCA75", label=T, repel=T, raster=F) + NoLegend() +
  (FeaturePlot(kmt2ar.rna, features="nCount_RNA", reduction="umapPCA75", label=F, repel=T, raster=F) + NoLegend())

VlnPlot(kmt2ar.rna, features="nCount_RNA", group.by="RNA_snn_res.2") + NoLegend() #49 and 55 have the highest counts 

FeaturePlot(kmt2ar.rna, features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), reduction="umapPCA75", ncol=3, label=T, repel=T, order=T)
DotPlot(kmt2ar.rna, features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"))

kmt2ar.rna@meta.data$seurat_clusters <- kmt2ar.rna@meta.data$RNA_snn_res.2
Idents(kmt2ar.rna) <- "seurat_clusters"
FeaturePlot(kmt2ar.rna, features=c("TFRC"), reduction="umapPCA75", label=T, repel=T, order=T, raster=F)

ggplot(kmt2ar.rna@meta.data) +
  geom_density(aes(x=nCount_RNA)) + 
  facet_wrap(~patientID)

# Annotation with marker genes  -------------

DimPlot(kmt2ar.rna, group.by="RNA_snn_res.2")

kmt2ar.rna@meta.data$Level1.anno <- as.character(kmt2ar.rna@meta.data$RNA_snn_res.2)
kmt2ar.rna@meta.data$Level1.anno[kmt2ar.rna@meta.data$Level1.anno %notin% c(15,12,22,43,38,29,44,48,67,27,28,19,16,34,49,31,32,51)] <- 'Blast'
kmt2ar.rna@meta.data$Level1.anno[kmt2ar.rna@meta.data$Level1.anno %in% c(15,12,22,43,38,29,44,48,67,27)] <- 'T/NK'
kmt2ar.rna@meta.data$Level1.anno[kmt2ar.rna@meta.data$Level1.anno %in% c(28,19)] <- 'B'
kmt2ar.rna@meta.data$Level1.anno[kmt2ar.rna@meta.data$Level1.anno %in% c(16,34,49)] <- 'Erythroid Progenitor'
kmt2ar.rna@meta.data$Level1.anno[kmt2ar.rna@meta.data$Level1.anno %in% c(31,32,51)] <- 'Mono'
#kmt2ar.rna@meta.data$Level1.anno[kmt2ar.rna@meta.data$Level1.anno %in% c("54")] <- 'DC'

myColors <- colorRampPalette(brewer.pal(6,"Set1"))(6)

#CD9

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/KMT2Ar/kmt2arRNA_87.pdf", height=10, width=10)
DimPlot(kmt2ar.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, reduction="umapPCA75", raster.dpi = c(1024,1024))
DimPlot(kmt2ar.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA75")
DimPlot(kmt2ar.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA75")
DimPlot(kmt2ar.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA75")
DimPlot(kmt2ar.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024), reduction="umapPCA75")+NoLegend()
DimPlot(kmt2ar.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024), reduction="umapPCA75")+NoLegend()
FeaturePlot(kmt2ar.rna, features = c("CD34"), reduction="umapPCA75", order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("CD63"), reduction="umapPCA75", order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("CD19"), reduction="umapPCA75", order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("MME"), reduction="umapPCA75",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("MS4A1"), reduction="umapPCA75",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("CD36"), reduction="umapPCA75",order=T, label=T, raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("CD3D"), reduction="umapPCA75",order=T, label=T, raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(kmt2ar.rna, features = c("LYZ"), reduction="umapPCA75",order=T, label=T, raster=T, raster.dpi = c(1024,1024), label.size = 3)
VlnPlot(kmt2ar.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
dev.off()

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/KMT2Ar/kmt2arRNA_AllMarkers.pdf",height=11.350, width=20)
FeaturePlot(kmt2ar.rna, features = c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"), order=T, label=T,  raster=T, raster.dpi = c(1024,1024), label.size = 3, ncol=3)
dev.off()

# Shannon Entropy Score -------------------

kmt2ar.rna <- FindNeighbors(kmt2ar.rna, dims = 1:100, reduction = 'pca')
kmt2ar.rna <- FindClusters(kmt2ar.rna, res = 0.25)
kmt2ar.rna <- FindClusters(kmt2ar.rna, res = 0.5)
kmt2ar.rna <- FindClusters(kmt2ar.rna, res = 0.75)
kmt2ar.rna <- FindClusters(kmt2ar.rna, res = 1)
kmt2ar.rna <- FindClusters(kmt2ar.rna, res = 2)
kmt2ar.rna = get.entropy(kmt2ar.rna,  c("RNA_snn_res.0.25","RNA_snn_res.0.5", "RNA_snn_res.0.75", "RNA_snn_res.1", "RNA_snn_res.2"),join.to.seurat = T)

kmt2ar.rna$seurat_clusters <- kmt2ar.rna$RNA_snn_res.0.5

FeaturePlot(kmt2ar.rna, features=c("average_entropy"), order=T, raster=F)
VlnPlot(kmt2ar.rna, features=c("average_entropy"), group.by = "Level1.anno", split.by = "seurat_clusters", ncol=10)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/KMT2Ar/ShannonEntropyFeaturePlots_kmt2ar.pdf", width = 10, height = 10)
FeaturePlot(kmt2ar.rna, features = c("RNA_snn_res.0.25_entropy","RNA_snn_res.0.5_entropy", "RNA_snn_res.0.75_entropy", "RNA_snn_res.1_entropy", "RNA_snn_res.2_entropy", "average_entropy"), 
            raster=T, raster.dpi = c(1024,1024), combine=F, reduction="umapPCA75")
dev.off()

# Similarity to Healthy Reference -----------

seurat.rna.hd <- readRDS(file="/mnt/isilon/tan_lab/chenc6/MLLr_Project/HTAN_Data_Submission/scRNASeq_Level4/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
seurat.rna.hd <- NormalizeData(seurat.rna.hd)
seurat.rna.hd <- FindVariableFeatures(seurat.rna.hd)
seurat.rna.hd <- ScaleData(seurat.rna.hd)

kmt2ar.rna.sub <- kmt2ar.rna[,sample(colnames(kmt2ar.rna), size = 50000, replace=F)]
kmt2ar.rna.sub <- NormalizeData(kmt2ar.rna.sub)
kmt2ar.rna.sub <- FindVariableFeatures(kmt2ar.rna.sub)
kmt2ar.rna.sub <- ScaleData(kmt2ar.rna.sub)

# Find anchors
anchors <- FindTransferAnchors(reference = seurat.rna.hd, query = kmt2ar.rna.sub, reduction = 'cca')
saveRDS(anchors, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/kmt2ar.healthy.donor.anchors.rds")

anchors <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/kmt2ar.healthy.donor.anchors.rds")

# Transfer labels
predictions.Ctype <- TransferData(anchorset = anchors, refdata = seurat.rna.hd$Ctype, k.weight = 50, weight.reduction = 'cca')
kmt2ar.rna.sub <- AddMetaData(object = kmt2ar.rna.sub, metadata = predictions.Ctype)
saveRDS(kmt2ar.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/kmt2ar_downsample50K_withHealthySimilarity.rds')

kmt2ar.rna.sub <- readRDS(file='/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Scripts/R_Objects/kmt2ar_downsample50K_withHealthySimilarity.rds')

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/KMT2Ar/scRNA_SimilarityScore_kmt2ar.pdf", width = 7, height = 7)
kmt2ar.rna.sub@meta.data %>%
  ggplot( aes(x=Level1.anno, y=prediction.score.max, fill=Level1.anno)) +
  geom_boxplot(outlier.size=1) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  ggtitle("KMT2A2 scRNA-Seq Similarity Score") +
  xlab("") + ylab("Similarity Score to healthy reference")
dev.off()

VlnPlot(kmt2ar.rna.sub, group.by = 'Level1.anno', features = c('prediction.score.max'))

FeaturePlot(kmt2ar.rna.sub, features=c("prediction.score.MEP"), split.by = "Level1.anno", ncol=2)
DimPlot(kmt2ar.rna.sub, group.by="predicted.id", split.by="Level1.anno", combine=F, ncol=3, label=T, repel=T)

FeaturePlot(kmt2ar.rna.sub, features=c("prediction.score.max"), order=T, ncol=2) + 
  DimPlot(kmt2ar.rna.sub, group.by="predicted.id") + 
  DimPlot(kmt2ar.rna.sub, group.by="Stage") + 
  DimPlot(kmt2ar.rna.sub, group.by="patientID")

FeaturePlot(kmt2ar.rna.sub, features=c("prediction.score.max","MME"), order=T, ncol=2) 

# Percentage of blasts -------

blasts.metadata <- kmt2ar.rna@meta.data %>% mutate(NormalStatus = ifelse(Level1.anno == "Mono" | Level1.anno == "N/K" | Level1.anno == "Erythroid Progenitor" | Level1.anno == "T/NK" | Level1.anno == "B", "Normal",
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

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_JointAssay/Plots/KMT2Ar/blast_count_similarity_kmt2ar.pdf", width = 5, height = 5)
ggplot(blasts.metadata, aes(x=BlastLabeledPercent,y=countedLeukemia)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE) + 
  labs(title="KMT2Ar", x="Annotated Blast Percentage", y="Counted Blast Percentage", caption=cor_coef) + 
  geom_text_repel(label=blasts.metadata$sampleID)
dev.off()


t.test(x=blasts.metadata$BlastLabeledPercent,y=blasts.metadata$countedLeukemia, two.tailed=TRUE, var.equal = TRUE)

# InferCNV -------------------------
kmt2ar.infercnv <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/KMT2Ar/KMT2Ar_finalInferCNV.rds")
kmt2ar.rna.downsample <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/KMT2Ar/KMT2Ar_downsampledSeuratforInferCNV.rds")

kmt2ar.rna.downsample <- NormalizeData(kmt2ar.rna.downsample)
kmt2ar.rna.downsample$Level1.anno <- NULL
metadata_to_add <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_Projection/Tables/metadata_blast_aggregateCtype_110723.csv", row.names=1)
rownames(metadata_to_add) <- metadata_to_add$patient_bc
kmt2ar.rna.downsample <- AddMetaData(kmt2ar.rna.downsample, metadata_to_add)

DimPlot(kmt2ar.rna.downsample, group.by="Level1.anno")

kmt2ar.rna.downsample <- add_to_seurat(seurat_obj = kmt2ar.rna.downsample, infercnv_output_path = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/InferCNV_Results/KMT2Ar")
DimPlot(kmt2ar.rna.downsample, group.by="top_dupli_4")

FeaturePlot(kmt2ar.rna.downsample, features="proportion_cnv_chr6")

VlnPlot(kmt2ar.rna.downsample, features=c("proportion_scaled_dupli_chr4"), group.by="Level1.anno")

ggplot(kmt2ar.rna.downsample@meta.data) + geom_boxplot(aes(x=Level1.anno,y=proportion_scaled_dupli_chr6))

