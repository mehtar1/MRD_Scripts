library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(RColorBrewer)
library(R.utils)
library(ggplot2)
library(tidyverse)
library(limma)
library(ggrastr)
library(ggpubr)
# source('scDataAnalysis_Utilities.R')
# source('dsAnalysis_utilities.R')
# source('DoubletsFinder.R')
source('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/entropycalc.R')
library(gdata)
# library(Matrix.utils)
library(clusterProfiler)
library(enrichplot)
`%notin%` = Negate(`%in%`)

iamp21.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/iamp21.rna.1.rds")
ph.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/ph.rna.1.rds")
kmt2ar.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/kmt2ar.rna.1.rds")
phlike.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/phlike.rna.1.rds")

myColors <- colorRampPalette(brewer.pal(6,"Set1"))(6)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/Level1_annotation_bySubtype_relabeled.pdf", height=11.350, width=20)
(DimPlot(phlike.rna, group.by="Level1.anno", ncol=2, cols = myColors) + ggtitle("Ph-Like")) + 
  (DimPlot(ph.rna, group.by="Level1.anno", cols = myColors, reduction="umapPCA50") + ggtitle("Ph-Positive")) +
  (DimPlot(kmt2ar.rna, group.by="Level1.anno", cols = myColors, reduction="umapPCA75") + ggtitle("KMT2Ar")) + 
  (DimPlot(iamp21.rna, group.by="Level1.anno", cols = myColors, reduction="umapPCA100") + ggtitle("iAMP21"))
dev.off()

(FeaturePlot(phlike.rna, features=c("CD34"), ncol=2, order=T) + ggtitle("Ph-Like")) + 
  (FeaturePlot(ph.rna, features=c("CD34"), reduction="umapPCA50", order=T) + ggtitle("Ph-Positive")) +
  (FeaturePlot(kmt2ar.rna, features=c("CD34"), reduction="umapPCA75", order=T) + ggtitle("KMT2Ar")) + 
  (FeaturePlot(iamp21.rna, features=c("CD34"), reduction="umapPCA100", order=T) + ggtitle("iAMP21"))

(VlnPlot(phlike.rna, features=c("CD81", "CD34"), group.by="Level1.anno", stack = T, flip = T) + ggtitle("Ph-Like")) + 
  (VlnPlot(ph.rna, features=c("CD81", "CD34"), group.by="Level1.anno", stack = T, flip = T) + ggtitle("Ph-Positive")) + 
  (VlnPlot(kmt2ar.rna, features=c("CD81", "CD34"), group.by="Level1.anno", stack = T, flip = T) + ggtitle("KMT2Ar")) + 
  (VlnPlot(iamp21.rna, features=c("CD81", "CD34"), group.by="Level1.anno", stack=T, flip = T) + ggtitle("iAMP21")) 

p1 <- (DimPlot(phlike.rna, group.by="Level1.anno", ncol=2, cols = myColors, raster=F) + ggtitle("Ph-Like")) + 
  (DimPlot(ph.rna, group.by="Level1.anno", cols = myColors, reduction="umapPCA50", raster=F) + ggtitle("Ph-Positive")) +
  (DimPlot(kmt2ar.rna, group.by="Level1.anno", cols = myColors, reduction="umapPCA75", raster=F) + ggtitle("KMT2Ar")) + 
  (DimPlot(iamp21.rna, group.by="Level1.anno", cols = myColors, reduction="umapPCA100", raster=F) + ggtitle("iAMP21"))

ggsave(p1, filename = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/Figures/Figure1/Subtype_UMAPs_Level1_080323.eps',
       device = 'eps', width = 12, height = 8)

  
#Possible markers:
#CD22 --> weak blast and B cell marker in Ph-like and iAMP21, absent in other subtypes for blasts 
#MME --> strong blast marker in all subtypes except KMT2Ar
#CD9 --> weak blast marker for all except KMT2Ar, strongest in iAMP21
#PTPRC/CD45 --> strongest blast marker in KMT2Ar, but expressed in all cell types for all subtypes 
#TDT/DNTT --> weak/medium blast specific marker in all subtypes, strongest in ph-positive/iamp21, but also expressed in monocytes for ph-pos
#CD81 --> blast marker in all subtypes, but present in iAMP21 b cells as well 
#CD34 --> blast marker, highest in iAMP21/ph-positive blasts 
#Next, try on a global.rna object

#global.rna object ------------------------------------------------

global.rna.1 <- merge(iamp21.rna, y = list(ph.rna, kmt2ar.rna, phlike.rna))
global.rna <- NormalizeData(global.rna)
global.rna <- FindVariableFeatures(global.rna, nfeatures = 4000) 
global.rna <- ScaleData(global.rna)
global.rna <- RunPCA(global.rna, npc = 100, verbose = T, features = VariableFeatures(global.rna)) 
global.rna <- RunUMAP(global.rna, dims = 1:100, reduction.name="umapPCA100") 
global.rna <- RunUMAP(global.rna, dims = 1:50, reduction.name="umapPCA50") 
global.rna <- FindNeighbors(global.rna, dims = 1:100)
global.rna <- FindClusters(res = 1.0)

saveRDS(global.rna.1, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/merged.postLabeling.011124.rds")

# global joint object 

global.joint <- readRDS('/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/jointModality_downsample0.2_withWNN.rds')

meta <- fread('/mnt/isilon/tan_lab/chenc6/HTAN_Leukemia_MRD/JointAssay/Scripts/metadata/Projection_MultiomeHD_WenbaoMethod/HTAN_MRD_AllSamples_RNAproj2healthyCells_to_MultiomeRef_UMAPpca20_Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels_Aug21_2023.csv')
rownames(meta) <- meta$patient_bc

global.joint <- AddMetaData(global.joint, meta)

global.joint <- FindMultiModalNeighbors(global.joint, reduction.list=c("pca", "lsi"), 
                                        dims.list = list(1:50, 2:50))
global.joint <- FindClusters(global.joint, res = 0.5, graph.name = "wsnn", algorithm = 3)
global.joint <- FindClusters(global.joint, res = 1, graph.name = "wsnn", algorithm = 3)
global.joint <- FindClusters(global.joint, res = 1.5, graph.name = "wsnn", algorithm = 3)
global.joint <- FindClusters(global.joint, res = 3, graph.name = "wsnn", algorithm = 3)
global.joint = get.entropy(global.joint,  c("wsnn_res.0.5","wsnn_res.1", "wsnn_res.3"),join.to.seurat = T)

global.joint$wsnn_res.0.5 <- as.numeric(global.joint$wsnn_res.0.5)
global.joint$wsnn_res.1 <- as.numeric(global.joint$wsnn_res.1)
global.joint$wsnn_res.3 <- as.numeric(global.joint$wsnn_res.3)

entropy.joint <- FeaturePlot(global.joint, features = c("wsnn_res.0.5_entropy","wsnn_res.1_entropy", "wsnn_res.3_entropy"), 
                             raster=T, reduction="umap_wsnn", ncol = 3) 

ggsave(entropy.joint, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/Figures/Figure1/GlobalJoint_entropy_032524.eps", 
       device = "pdf", height=6, width=12, bg='transparent')

### Fusion status 
fusion.status <- read_xlsx("/mnt/isilon/tan_lab/mehtar1/MRD_RM/References/SampleProcessing_HTAN_MRD.xlsx")
fusion.status <- fusion.status[-1,c(3,7)]
fusion.status <- distinct(fusion.status)
colnames(fusion.status) <- c("Fusion", "PT.number")

meta <- left_join(global.joint@meta.data, fusion.status, by="PT.number")
rownames(meta) <- meta$cell.name
global.joint <- AddMetaData(global.joint, meta)

fusion.joint <- DimPlot(global.joint, group.by = "Fusion", 
                             raster=T, reduction="umap_wsnn") 

ggsave(fusion.joint, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/Figures/Figure1/GlobalJoint_FusionStatus_040524.eps", 
       device = "pdf", height=6, width=8, bg='transparent')


























DimPlot(global.rna, group.by = "Level1.anno", reduction="umapPCA100", raster=F) + 
  (DimPlot(global.rna, group.by = "RNA_snn_res.1", raster=F) + NoLegend())
DimPlot(global.rna, group.by = "Level1.anno", reduction="umapPCA50", raster=F)

blast_list <- global.rna@meta.data %>% filter(Level1.anno == "Blast") %>% rownames() 
merged_onlyBlast <- subset(global.rna, cells = blast_list)
merged_onlyBlast <- NormalizeData(merged_onlyBlast)

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/GlobalRNA_PCA100_postLabelingMerged.pdf", height=10, width=10)
DimPlot(global.rna, group.by = 'Level1.anno', label = F, cols = myColors, raster=T, raster.dpi = c(1024,1024), reduction="umapPCA100")
DimPlot(global.rna, group.by = 'Stage', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA100")
DimPlot(global.rna, group.by = 'Subtype', cols = myColors[1:4], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA100")
DimPlot(global.rna, group.by = 'MRDGroup', cols = myColors[1:3], raster=T, raster.dpi = c(1024,1024), reduction="umapPCA100")
DimPlot(global.rna, group.by = "patientID", label=T, repel=T, raster.dpi=c(1024,1024), reduction="umapPCA100")+NoLegend()
DimPlot(global.rna, group.by = "patientID", label=F, raster.dpi=c(1024,1024), reduction="umapPCA100")+NoLegend()
FeaturePlot(global.rna, features = c("CD19"), order=T, label=T, reduction="umapPCA100", raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(global.rna, features = c("MME"), order=T, label=T, reduction="umapPCA100", raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(global.rna, features = c("MS4A1"), order=T, label=T, reduction="umapPCA100", raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(global.rna, features = c("CD36"), order=T, label=T, reduction="umapPCA100", raster=T, raster.dpi = c(1024,1024),label.size = 3)
FeaturePlot(global.rna, features = c("CD3D"), order=T, label=T, reduction="umapPCA100", raster=T, raster.dpi = c(1024,1024), label.size = 3)
FeaturePlot(global.rna, features = c("LYZ"), order=T, label=T, reduction="umapPCA100", raster=T, raster.dpi = c(1024,1024), label.size = 3)
VlnPlot(global.rna,features=c("CD19", "MME", "MS4A1", "CD36", "CD3D", "LYZ"),stack = T, group.by="Level1.anno", pt.size=0, flip=T)
VlnPlot(merged_onlyBlast,features=c("CD22", "MME", "CD9", "PTPRC", "DNTT", "CD81", "CD34"),stack = T, group.by="Subtype", pt.size=0, flip=T)
dev.off()

#CD22 --> weak blast and B cell marker in Ph-like and iAMP21, absent in other subtypes for blasts 
#MME --> strong blast marker in all subtypes except KMT2Ar
#CD9 --> weak blast marker for all except KMT2Ar, strongest in iAMP21
#PTPRC/CD45 --> strongest blast marker in KMT2Ar, but expressed in all cell types for all subtypes 
#TDT/DNTT --> weak/medium blast specific marker in all subtypes, strongest in ph-positive/iamp21, but also expressed in monocytes for ph-pos
#CD81 --> blast marker in all subtypes, but present in iAMP21 b cells as well 
#CD34 --> blast marker, highest in iAMP21/ph-positive blasts 

## Joint Embedding --------------------------------------------------------------------------

global.atac.ds <- readRDS("/mnt/isilon/tan_lab/chenc6/HTAN_Leukemia_MRD/JointAssay/Scripts/SeuratObj/seurat_ATAC_FullDataset_Downsample0.2.rds")
global.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/merged.postLabeling.rds")

ds.cell.names <- rownames(global.atac.ds@meta.data)

global.rna.ds <- subset(global.rna, cells = ds.cell.names)

joint.ds <- CreateSeuratObject(
  counts = global.rna.ds@assays$RNA@data,
  assay = "RNA"
)

joint.ds[["RNA"]] <-global.rna.ds[["RNA"]]
joint.ds[["ATAC"]] <- global.atac.ds[["ATAC"]]
joint.ds@reductions <- append(global.rna.ds@reductions,global.atac.ds@reductions)

joint.ds <- FindMultiModalNeighbors(
  object = joint.ds,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0
)

joint.ds <- RunUMAP(
  object = joint.ds,
  nn.name = "weighted.nn",
  reduction.name = "umap_wsnn",
  assay = "RNA",
  verbose = TRUE
)

metadata_forJoint <- read_csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/Integration_and_WNN/Scripts/Split_by_Patient/Tables/metadata_allcells_updated.csv")
metadata_forJoint <- metadata_forJoint %>% select(Stage, Level1.anno, S_method, W_method, MRD, cell.name, S.Score, G2M.Score, PT.number)
rownames(metadata_forJoint) <- metadata_forJoint$cell.name
joint.ds <- AddMetaData(joint.ds, metadata = metadata_forJoint)

priming.meta <- read.delim("/mnt/isilon/tan_lab/chenc6/HTAN_Leukemia_MRD/JointAssay/Scripts/metadata/HTAN_MRD_AllSamples_ATACproj2healthyCells_to_MultiomeRef_UMAPlsi25_Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels_Aug21_2023.csv")
rownames(priming.meta) <- priming.meta$patient_bc
joint.ds <- AddMetaData(joint.ds, metadata = priming.meta)

joint.ds$joint.embedding.Level1 <- joint.ds$Level1.anno

joint.ds@meta.data <- joint.ds@meta.data %>% 
  mutate(joint.embedding.Level1 = case_when(Level1.anno=='Blast' & Stage=="IDX" ~ 'IDX Blast',
                                            Level1.anno=='Blast' & Stage!="IDX" ~ 'MRD/RLP Blast', 
                                            Level1.anno!="Blast" ~ Level1.anno))

saveRDS(joint.ds, file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/jointModality_downsample0.2_withWNN.rds")

p3 <- DimPlot(joint.ds, reduction = "umap_wsnn", group.by="joint.embedding.Level1") 

p4 <- DimPlot(joint.ds, group.by="Patient.Sample", reduction = "umap_wsnn", order=T, shuffle=T, label=T, repel=T, label.size = 2) + NoLegend()

p6 <- DimPlot(joint.ds, group.by="Projected.Ctype.Wmethod", reduction = "umap_wsnn", order=T, shuffle=T)

p7 <- DimPlot(joint.ds, group.by="predicted.celltype.lsi", reduction = "umap_wsnn", order=T, shuffle=T)


ggsave(p3, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/JointEmbedding_global.eps", 
       device = "eps", height=10, width=11)

ggsave(p4, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/JointEmbedding_global_PatientSamples.eps", 
       device = "eps", height=10, width=11)

ggsave(p6, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/JointEmbedding_global_RNAProjection.eps", 
       device = "eps", height=10, width=11)

ggsave(p7, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/JointEmbedding_global_ATACProjection.eps", 
       device = "eps", height=10, width=11)


# Shannon Entropy Score -------------------

global.rna <- FindNeighbors(global.rna, dims = 1:100, reduction = 'pca')
global.rna <- FindClusters(global.rna, res = 0.5)
global.rna <- FindClusters(global.rna, res = 1)
global.rna <- FindClusters(global.rna, res = 1.5)
global.rna <- FindClusters(global.rna, res = 3)
global.rna = get.entropy(global.rna,  c("RNA_snn_res.0.5","RNA_snn_res.1", "RNA_snn_res.3"), join.to.seurat = T)

global.rna$RNA_snn_res.0.5_entropy <- as.numeric(global.rna$RNA_snn_res.0.5_entropy)
global.rna$RNA_snn_res.1_entropy <- as.numeric(global.rna$RNA_snn_res.1_entropy)
global.rna$RNA_snn_res.3_entropy <- as.numeric(global.rna$RNA_snn_res.3_entropy)


pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/ShannonEntropyFeaturePlots_global.pdf", width = 8, height = 15)
FeaturePlot(global.rna, features = c("RNA_snn_res.0.5_entropy","RNA_snn_res.1_entropy", "RNA_snn_res.3_entropy"), 
            raster=T, raster.dpi=c(1024,1024), reduction="umapPCA100")
dev.off()


# Dotplot Markers  -------------------------------------------------------------

# T/NK cells --> CD3D
# Myeloid cells --> CD36
# Healthy B cells --> MS4A1
# EP cells --> CD36

global.rna <- readRDS(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/merged.postLabeling.011124.rds")
metadata_to_add <- read.csv(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_Projection/Tables/metadata_blast_aggregateCtype_110723.csv", row.names=1)
rownames(metadata_to_add) <- metadata_to_add$patient_bc
#Relabel metadata 
global.metadata <- global.rna@meta.data %>% select(Level1.anno, Subtype, Stage)
global.metadata$'Cell.Type' <- global.metadata$Level1.anno

global.metadata <- global.metadata %>%
  mutate('Cell.Type' = case_when(Level1.anno=="Blast" & Stage=="IDX" & Subtype=='Ph-like' ~ 'Ph-Like IDX',
                                 Level1.anno=="Blast" & Stage=="IDX" & Subtype=='Ph+' ~ 'Ph-Positive IDX',
                                 Level1.anno=="Blast" & Stage=="IDX" & Subtype=='KMT2Ar' ~ 'KMT2Ar IDX',
                                 Level1.anno=="Blast" & Stage=="IDX" & Subtype=='iAMP21' ~ 'iAMP21 IDX',
                                 Level1.anno=="Blast" & Stage=="MRD" & Subtype=='Ph-like' ~ 'Ph-Like MRD',
                                 Level1.anno=="Blast" & Stage=="MRD" & Subtype=='Ph+' ~ 'Ph-Positive MRD',
                                 Level1.anno=="Blast" & Stage=="MRD" & Subtype=='KMT2Ar' ~ 'KMT2Ar MRD',
                                 Level1.anno=="Blast" & Stage=="MRD" & Subtype=='iAMP21' ~ 'iAMP21 MRD',
                                 Level1.anno=="Blast" & Stage=="RLP" & Subtype=='Ph-like' ~ 'Ph-Like RLP',
                                 Level1.anno=="Blast" & Stage=="RLP" & Subtype=='Ph+' ~ 'Ph-Positive RLP',
                                 Level1.anno=="Blast" & Stage=="RLP" & Subtype=='KMT2Ar' ~ 'KMT2Ar RLP',
                                 Level1.anno=="Blast" & Stage=="RLP" & Subtype=='iAMP21' ~ 'iAMP21 RLP',
                                 Level1.anno=="B" ~ 'B', 
                                 Level1.anno=="Erythroid Progenitor" ~ 'Erythroid Progenitor', 
                                 Level1.anno=="Mono" ~ 'Mono',
                                 Level1.anno=="T/NK" ~ 'T/NK')) 

global.rna <- AddMetaData(global.rna, metadata=global.metadata$'Cell.Type', col.name= 'Cell.Type')
global.rna <- AddMetaData(global.rna, metadata=metadata_to_add)

metadata_to_add <- global.rna@meta.data %>% mutate(Blast_Lineage = case_when(
  Projected.Ctype.Wmethod %in% c("cDC", "DC-progenitor", "GMP", "MEP", "Monocyte", "pDC") & Level1.anno=="Blast" ~ "Myeloid-like",
  Projected.Ctype.Wmethod %in% c("HSC/MPP") & Level1.anno=="Blast" ~ "HSC-like",
  Projected.Ctype.Wmethod %in% c("NK", "T-cell") & Level1.anno=="Blast" ~ NA,
  Projected.Ctype.Wmethod %in% c("CLP", "Immature-B", "LMPP", "Mature-B-1", 
                                 "Mature-B-2", "Pre-B", "Pre-pro-B", "Pro-B") & Level1.anno=="Blast" ~ "B-lineage", 
  .default = Level1.anno
))

global.rna <- AddMetaData(global.rna, metadata=metadata_to_add)

table(global.rna$Level1.anno, global.rna$Cell.Type, useNA = "always")

global.rna$Cell.Type <- factor(global.rna$Cell.Type, 
                            levels=c("Ph-Like IDX", "Ph-Like MRD", 
                                     "Ph-Positive IDX","Ph-Positive MRD", "Ph-Positive RLP", 
                                     "KMT2Ar IDX","KMT2Ar MRD", "KMT2Ar RLP", 
                                     "iAMP21 IDX", "iAMP21 MRD", "iAMP21 RLP",  
                                     "B", "T/NK", "Mono", "Erythroid Progenitor"))

global.rna <- NormalizeData(global.rna)
# global.rna <- ScaleData(global.rna)

to.plot <- c("MME", "CD19", "CD34", "CD38", 
             "MS4A1", "CD3D", "LYZ", "CD36")

p2 <- DotPlot(global.rna, features = to.plot, scale = TRUE, col.min = -1, 
        col.max = 1,scale.min = 0, scale.max = 80, group.by="Cell.Type") + 
  RotatedAxis() + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

ggsave(p2, filename = "/mnt/isilon/tan_lab/mehtar1/MRD_RM/Figures/Figure1/All87Samples_DotPlot.eps", 
       device = "eps", height=5, width=7)

global.rna$Blast_Lineage <- factor(global.rna$Blast_Lineage, 
                               levels=c("T/NK","Erythroid Progenitor", "Mono","B","Myeloid-like", "HSC-like","B-lineage"))
Idents(global.rna) <- "Cell.Type"
p3 <- DotPlot(global.rna, features = c("CD19", "MS4A1", "CD22", "CD24", "CD79A", 
                                       "CD7","PTPRC", "ITGAM", "CD33", "CD34"), 
              scale = TRUE, col.min = -1, 
              col.max = 1,scale.min = 0, scale.max = 80, 
              group.by="Blast_Lineage") + 
  RotatedAxis() + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

p4 <- DotPlot(global.rna, features = to.plot, scale = TRUE, col.min = -1, 
              col.max = 1,scale.min = 0, scale.max = 80, group.by="Blast_Lineage") + 
  RotatedAxis() + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


# DotPlot(global.rna, features = to.plot, scale = TRUE, col.min = -1, 
#         col.max = 1,scale.min = 0, scale.max = 80, group.by="Level1.anno") + RotatedAxis()+ 
#   scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
# 
# RidgePlot(global.rna, features=to.plot, group.by = "Cell.Type")

pdf(file="/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/BlastMarkers_global.pdf", width = 10, height = 10)
FeaturePlot(global.rna, features=c("MME", "CD19"), order=T, raster=F)
# FeaturePlot(global.rna, features="CD34", order=T, raster=F)
# FeaturePlot(global.rna, features="CD38", order=T, raster=F)
FeaturePlot(global.rna, features="MS4A1", order=T, raster=F) + DimPlot(global.rna, group.by="Level1.anno",raster=F)
FeaturePlot(global.rna, features="CD36", order=T, raster=F) + DimPlot(global.rna, group.by="Level1.anno",raster=F)
FeaturePlot(global.rna, features="CD3D", order=T, raster=F) + DimPlot(global.rna, group.by="Level1.anno",raster=F)
dev.off()

Idents(global.rna) <- "Cell.Type"
# markers.to.plot <- FindAllMarkers(global.rna, only.pos = T, logfc.threshold = 0.5)
# write.csv(markers.to.plot, "/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Tables/markers.global.dotplot.csv")

### Similarity to Healthy Reference ---------------------------------------------

seurat.rna.hd <- readRDS(file="/mnt/isilon/tan_lab/chenc6/MLLr_Project/HTAN_Data_Submission/scRNASeq_Level4/seurat_pool_logNorm_gini_FiveHD_10Xv3_downsample10000HSPC_Oct25_2021.rds")
seurat.rna.hd <- NormalizeData(seurat.rna.hd)
seurat.rna.hd <- FindVariableFeatures(seurat.rna.hd)
seurat.rna.hd <- ScaleData(seurat.rna.hd)

global.rna.sub <- global.rna[,sample(colnames(global.rna), size = 50000, replace=F)]
global.rna.sub <- NormalizeData(global.rna.sub)
global.rna.sub <- FindVariableFeatures(global.rna.sub)
global.rna.sub <- ScaleData(global.rna.sub)

# Find anchors
anchors <- FindTransferAnchors(reference = seurat.rna.hd, query = global.rna.sub, reduction = 'pcaproject')
# Transfer labels
predictions.Ctype <- TransferData(anchorset = anchors, refdata = seurat.rna.hd$Ctype, k.weight = 50, weight.reduction = 'pcaproject')
global.rna.sub <- AddMetaData(object = global.rna.sub, metadata = predictions.Ctype)
saveRDS(global.rna.sub, file = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/global_downsample50K_withHealthySimilarity.rds')
global.rna.sub <- readRDS(file='/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/global_downsample50K_withHealthySimilarity.rds')

p1 <- global.rna.sub@meta.data %>%
  ggplot( aes(x=Level1.anno, y=prediction.score.max, fill=Level1.anno)) +
  geom_boxplot(outlier.size=0) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  ggtitle("Global scRNA-Seq Similarity Score") +
  xlab("") + ylab("Similarity Score to healthy reference") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p1

ggsave(p1, filename = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/Figures/Supplemental/QC/Global_SimilaritytoHealthyReference.eps',
       device = 'eps', width = 5, height = 5)

p2 <- global.rna.sub@meta.data %>%
  ggplot( aes(x=Level1.anno, y=prediction.score.max, fill=Level1.anno)) +
  geom_boxplot(outlier.size=0, outlier.shape = NA) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
  ggtitle("Global scRNA-Seq Similarity Score") +
  xlab("") + ylab("Similarity Score to healthy reference") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p2

ggsave(p1, filename = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/Figures/Supplemental/QC/Global_SimilaritytoHealthyReference_NoOutlier.eps',
       device = 'eps', width = 5, height = 5)

# Entropy calculation between stages ------------------------------------------

global.blasts <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/RNA_DEG/Scripts/R_Objects/global.rna.blasts.onlyidxmrd.011124.rds")

entropy.sample <- get.entropy(global.blasts,  c("Patient.Sample"), join.to.seurat = F)[[1]]

entropy.sample.meta <- entropy.sample %>% left_join(
  global.blasts@meta.data %>% 
    dplyr::select("Patient.Sample", "Stage", "Subtype", "PT.number", "MRD", "MRDGroup"), by="Patient.Sample") %>% unique()

entropy.sample.meta <- entropy.sample.meta %>% filter(PT.number %notin% c("PT5", "PT17"))

p1 <- ggpaired(entropy.sample.meta, x="Stage", y="Patient.Sample_entropy", facet.by = "Subtype", id="PT.number", 
         line.color = "gray", line.size = 0.4, color="Stage", palette = "npg", ncol = 4) + 
  stat_compare_means(paired = T)
  ggsave(p1, filename = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/ChangeInEntropyByStage_SplitBySubtype_102224.eps',
         device = 'eps', width = 8, height = 5)
  
p2 <- ggpaired(entropy.sample.meta, x="Stage", y="Patient.Sample_entropy", id="PT.number", 
         line.color = "gray", line.size = 0.4, color="Stage", palette = "npg") + 
  stat_compare_means(paired = T)
  ggsave(p2, filename = '/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Plots/Global/ChangeInEntropyByStage_Global_102224.eps',
         device = 'eps', width = 5, height = 5)

p3 <- ggpaired(entropy.sample.meta, x="Stage", y="Patient.Sample_entropy", facet.by = "MRDGroup", id="PT.number", 
           line.color = "gray", line.size = 0.4, color="Stage", palette = "npg", ncol = 4) + 
    stat_compare_means(paired = T)
  
p4 <- ggboxplot(entropy.sample.meta %>% dplyr::filter(Stage=="IDX"), x="MRDGroup", y="Patient.Sample_entropy")
p5 <- ggboxplot(entropy.sample.meta %>% dplyr::filter(Stage=="MRD"), x="MRDGroup", y="Patient.Sample_entropy") + stat_compare_means()

# Test integration by Stage -----

global.ds <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_MalignancyEvidence/Scripts/R_Objects/jointModality_downsample0.2_withWNN.rds")
metadata_to_add <- read.csv("/mnt/isilon/tan_lab/mehtar1/MRD_RM/MRD_Projection/Tables/UpdatedProjectionResults_AllCells_100125.csv", 
                                               row.names = 1)

global.ds <- AddMetaData(global.ds, metadata_to_add)
global.ds@reductions[["umap_wsnn"]] <- NULL
global.ds@reductions[["umapPCA100"]] <- NULL

global.ds <- UpdateSeuratObject(global.ds)

DimPlot(global.ds, group.by="Stage")

library(future)
plan(multicore, workers = 8)
options(future.globals.maxSize = 10 * 1024^3)

# List of Seurat objects, already filtered & normalized or with raw counts
objs <- SplitObject(global.ds, "Stage")

# 1. Normalize + find variable features (log-normalization example)
objs <- lapply(objs, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  x
})

# 2. Select integration features
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)

# Optional: SCT workflow instead of log-normal
# objs <- PrepSCTIntegration(object.list = objs, anchor.features = features)

# 3. Find anchors (CCA or RPCA)
objs <- lapply(objs, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  x
})

anchors <- FindIntegrationAnchors(
  object.list = objs,
  anchor.features = features,
  reduction = "rpca",         # use "cca" for CCA
  dims = 1:30
)

# 4. Integrate data
options(future.globals.maxSize = 8 * 1024^3)
integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  new.assay.name = "integrated_RNA"
)

DefaultAssay(integrated) <- "integrated_RNA"

# 5. Downstream analysis
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated, npcs = 30)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)







