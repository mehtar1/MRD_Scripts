# Goal: QC and analysis of CITE-Seq data 
# Date: 11/19/25

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(patchwork)
#BiocManager::install("dsb")
library(dsb)
library(ComplexHeatmap)
library(data.table)
library(RColorBrewer)
library(RSpectra)
library(uwot)

options(future.globals.maxSize =  100 * 1000 * 1024^2)


## Functions -----

adt_heatmap <- function(seurat_obj,
                        prots = rownames(seurat_obj@assays$ADT@features)[1:37],
                        groupByVar,
                        clusterIdents = T,
                        title=NULL){
  
  require(ComplexHeatmap)
  
  s <- seurat_obj
  s <- JoinLayers(s, assay="ADT")
  s <- NormalizeData(s, assay="ADT", normalization.method="CLR")
  cite_data = GetAssayData(s,slot = 'data',assay = 'ADT') %>% t()
  
  # combine into dataframe
  d = cbind(s@meta.data, cite_data)

  # calculate the median protein expression per cluster
  dat_plot = d %>%
    dplyr::group_by(.data[[groupByVar]]) %>%
    dplyr::summarize_at(.vars = c(prots), .funs = median) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(groupByVar)
  
  prot_col = circlize::colorRamp2(breaks = seq(0.25,2, by = 0.0625),
                                  colors = viridis::viridis(n = 29, option = "B"))
  
  p1 = Heatmap(t(dat_plot)[prots, ],
               name = "protein",
               col = prot_col,
               use_raster = T,
               cluster_columns = clusterIdents, 
               row_names_gp = gpar(color = "black", fontsize = 10),
               column_title = title, column_title_gp = gpar(fontsize = 12, col = "black"))
  
  print(p1)
  return(p1)
  
}

cluster.plot <- function(table, x.var, y.var, return.table=FALSE, color.scale=NULL){
  require(dplyr)
  cluster.freq.sample <- table %>%
    dplyr::select({{x.var}}, {{y.var}}) %>% group_by({{x.var}}) %>%
    mutate(total=n()) %>% group_by({{x.var}}, {{y.var}}) %>% mutate(freq=n()/total) %>%
    dplyr::select({{x.var}}, {{y.var}}, freq) %>% distinct()
  
  totals <- table %>% dplyr::select({{x.var}}) %>% 
    group_by({{x.var}}) %>% mutate(total_n=n()) %>% unique()
  
  if(is.null(color.scale)){
    p1 <- ggplot(cluster.freq.sample) +
      geom_col(aes(x={{x.var}}, y=freq, fill={{y.var}}), color="gray") +
      geom_text(data = totals, aes(x = {{x.var}}, y = 1, label = total_n), vjust = -0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 
      
  }else{
    p1 <- ggplot(cluster.freq.sample) +
      geom_col(aes(x={{x.var}}, y=freq, fill={{y.var}}), color="gray") +
      geom_text(data = totals, aes(x = {{x.var}}, y = 1, label = total_n), vjust = -0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
      
    p1 <- p1 + scale_fill_manual(values=color.scale)
  }

  if(return.table){
    return(list(p1, cluster.freq.sample))
  } else{
    return(p1)
  }
}

# Data loading ----

hd.ref <- readRDS('/mnt/isilon/tan_lab/xuj5/HealthyReference/other_r_objects/thymus_bm/2021_7_26_full.ref.gini.genes.cc.removed.rds')

merged.all <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/5Pairs_IntegratedAnnotated_OnlyRepeatSamples_111925.rds")

infercnv.6099R <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/6099R_InferCNV.rds")

infercnv.6061 <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/6061_InferCNV.rds")

infercnv.6146 <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/6146_InferCNV.rds")

infercnv.6092 <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/6092_InferCNV.rds")

# inferCNV for annotation ------ 

infercnv.list <- list(infercnv.6099R, infercnv.6061, infercnv.6146, infercnv.6092)

# Generate plots for each sample 

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/inferCNVResults_5Pairs_IDXMRDPairs_111925.pdf",
    height=7, width=8)

# Merged representation by sample ----

# Pair 4: 6099R (IDX) and 6146 (MRD)
pair4.s <- merge(infercnv.6099R, infercnv.6146)
pair4.s <- pair4.s %>% ScaleData() %>% FindVariableFeatures() %>%
  RunPCA(ndims=40) %>% FindNeighbors() %>%
  FindClusters(resolution = 1.0, cluster.name="merged_res_1.0") %>% RunUMAP(dims = 1:40)

DimPlot(pair4.s, group.by="Sample", ncol=1) + DimPlot(pair4.s, group.by="has_dupli_chr21")
DimPlot(pair4.s, group.by="Sample", ncol=1) + DimPlot(pair4.s, group.by="celltype.binned")
#DimPlot(pair4.s, group.by="Sample", ncol=1) + (DimPlot(pair4.s, group.by="Level1.anno") + ggtitle("Original Anno"))
DimPlot(pair4.s, group.by="celltype.binned", ncol=1) + (DimPlot(pair4.s, group.by="merged_res_1.0", label=T, repel=T) + NoLegend())
FeaturePlot(pair4.s, features=c("proportion_dupli_chr21", "proportion_dupli_chrX", "proportion_loss_chrX"), order=T, ncol=2) & NoAxes()
VlnPlot(pair4.s, features=c("proportion_dupli_chr21", "proportion_dupli_chrX", "proportion_loss_chrX"),
        stack=T, flip = T, group.by="merged_res_1.0") + NoLegend()

pair4.s@meta.data <- mutate(pair4.s@meta.data, Level1.anno = 
                              case_when(merged_res_1.0 %in% c(13) ~ "Myeloid",
                                        merged_res_1.0 %in% c(6, 9, 14, 15) ~ "T_NK",
                                        merged_res_1.0 %in% c(7,8) ~ "B",
                                        .default = "Blast"))
DimPlot(pair4.s, group.by="Sample", ncol=1) + (DimPlot(pair4.s, group.by="Level1.anno") + ggtitle("New Anno"))

FeaturePlot(pair4.s, features=c("MME", "CLEC12A", "MS4A1", "CD3E")) & NoAxes()
FeaturePlot(pair4.s, features=c("CD19", "CD33", "MS4A1", "CD3E"), order=T) & NoAxes()


pair4.s@meta.data <- mutate(pair4.s@meta.data, Blast_Lineage = case_when(
  Level1.anno=="Blast" & celltype.short %in% hspc.ctypes ~ "HSPC-like", 
  Level1.anno=="Blast" & celltype.short %in% my.ctypes ~ "Myeloid-like", 
  Level1.anno=="Blast" & celltype.short %in% b.ctypes ~ "B-like", 
  Level1.anno=="Blast" & celltype.short %in% t.ctypes ~ "T-cell-like", 
  .default = Level1.anno
))

pair4.s$Blast_Lineage <- factor(pair4.s$Blast_Lineage, 
                                   levels=c("B-like", "Myeloid-like", "HSPC-like", "T-cell-like", "B", "T_NK", "Myeloid"))
lineage.cols <- c("#377eb8","#e41a1c", "#4daf4a", "#984ea3", "#cbd5e8", "#fdcdac", "#f4cae4")

DimPlot(pair4.s, reduction="umap", group.by="Stage", ncol=1) + 
  DimPlot(pair4.s, reduction="umap", group.by="Blast_Lineage", cols = lineage.cols)

## Reannotate cells by cluster 

# Pair 3: 6092 (IDX) and 6061 (MRD)
pair3.s <- merge(infercnv.6092, infercnv.6061)
pair3.s <- pair3.s %>% ScaleData() %>% FindVariableFeatures() %>%
  RunPCA(ndims=40) %>% FindNeighbors() %>%
  FindClusters(resolution = 1.0, cluster.name="merged_res_1.0") %>% RunUMAP(dims = 1:40)

DimPlot(pair3.s, group.by="Sample", ncol=1) + DimPlot(pair3.s, group.by="celltype.short")
DimPlot(pair3.s, group.by="Sample", ncol=1) + DimPlot(pair3.s, group.by="celltype.binned")
DimPlot(pair3.s, group.by="Sample", ncol=1) + (DimPlot(pair3.s, group.by="Level1.anno") + ggtitle("Original Anno"))

DimPlot(pair3.s, group.by="celltype.binned", ncol=1) + (DimPlot(pair3.s, group.by="merged_res_1.0", label=T, repel=T) + NoLegend())

FeaturePlot(pair3.s, features=c("nCount_RNA"))

FeaturePlot(pair3.s, features=c("proportion_dupli_chr4", "proportion_dupli_chr8", "proportion_dupli_chr21", "proportion_dupli_chr22"), 
            ncol=2) & NoAxes()
VlnPlot(pair3.s, features=c("proportion_dupli_chr4", "proportion_dupli_chr8", "proportion_dupli_chr21", "proportion_dupli_chr22"),
        stack=T, flip = T, group.by="merged_res_1.0") + NoLegend()
pair3.s@meta.data <- mutate(pair3.s@meta.data, Level1.anno = 
                              case_when(merged_res_1.0 %in% c(10,16,17) ~ "Myeloid",
                                        merged_res_1.0 %in% c(11,15) ~ "T_NK",
                                        merged_res_1.0 %in% c(18) ~ "B",
                                        .default = "Blast"))
DimPlot(pair3.s, group.by="Sample", ncol=1) + (DimPlot(pair3.s, group.by="Level1.anno") + ggtitle("New Anno"))

FeaturePlot(pair3.s, features=c("CD19", "CD33", "MS4A1", "CD3E"), order=T) & NoAxes()

pair3.s@meta.data <- mutate(pair3.s@meta.data, Blast_Lineage = case_when(
  Level1.anno=="Blast" & celltype.short %in% hspc.ctypes ~ "HSPC-like", 
  Level1.anno=="Blast" & celltype.short %in% my.ctypes ~ "Myeloid-like", 
  Level1.anno=="Blast" & celltype.short %in% b.ctypes ~ "B-like", 
  Level1.anno=="Blast" & celltype.short %in% t.ctypes ~ "T-cell-like", 
  .default = Level1.anno
))

pair3.s$Blast_Lineage <- factor(pair3.s$Blast_Lineage, 
                                levels=c("B-like", "Myeloid-like", "HSPC-like", "T-cell-like", "B", "T_NK", "Myeloid"))

DimPlot(pair3.s, reduction="umap", group.by="Stage", ncol=1) + 
  DimPlot(pair3.s, reduction="umap", group.by="Blast_Lineage", cols = lineage.cols)

dev.off()

# Add in metadata for CNVs and new Level1.anno from each object to the merged object ----

metadata.subset <- c("has_dupli_chr21", "proportion_cnv_chr21", "proportion_dupli_chr21",
                     "has_dupli_chr22", "proportion_cnv_chr22", "proportion_dupli_chr22",
                     "has_dupli_chrX", "proportion_cnv_chrX", "proportion_dupli_chrX",
                     "has_loss_chrX", "proportion_loss_chrX",
                     "has_dupli_chr4", "proportion_cnv_chr4", "proportion_dupli_chr4",
                     "Level1.anno")

pairs.infercnv.list <- list(pair3.s, pair4.s)

for(i in 1:length(pairs.infercnv.list)){
  current.obj <- pairs.infercnv.list[[i]]
  meta <- current.obj@meta.data %>% dplyr::select(all_of(metadata.subset))
  merged.all <- AddMetaData(merged.all, meta)
}

# Look at CNVs in merged and integrated object ---

#merged.all <- RunUMAP(merged.all, dims=1:40)

merged.all@meta.data <- mutate(merged.all@meta.data, Blast_Lineage = case_when(
  Level1.anno=="Blast" & celltype.short %in% hspc.ctypes ~ "HSPC-like", 
  Level1.anno=="Blast" & celltype.short %in% my.ctypes ~ "Myeloid-like", 
  Level1.anno=="Blast" & celltype.short %in% b.ctypes ~ "B-like", 
  Level1.anno=="Blast" & celltype.short %in% t.ctypes ~ "T-cell-like", 
  .default = Level1.anno
))

merged.all$Blast_Lineage <- factor(merged.all$Blast_Lineage, 
                                   levels=c("B-like", "Myeloid-like", "HSPC-like", "T-cell-like", "B", "T_NK", "Myeloid"))

merged.all <- RunUMAP(merged.all, dims = 1:40)

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/inferCNV_5Pairs_SamplesInFullCohort.pdf",
    height=5, width=7)

DimPlot(merged.all, reduction="umap.rpca.pc40", group.by="Level1.anno")
DimPlot(merged.all, reduction="umap", group.by="Level1.anno")
DimPlot(merged.all, reduction="umap", group.by="Blast_Lineage", cols=lineage.cols)
DimPlot(merged.all, reduction="umap", group.by="Sample")
FeaturePlot(merged.all, reduction="umap", features=c("nCount_RNA"))
DimPlot(merged.all, reduction="umap.rpca.pc40", group.by="has_dupli_chr21")
DimPlot(merged.all, reduction="umap.rpca.pc40", group.by="has_dupli_chrX")

DimPlot(merged.all, reduction="umap.rpca.pc40", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6099R" & has_dupli_chr21 == T)),
        sizes.highlight = 0.5) +
  ggtitle("6099R, Duplication in chr21, Integrated")

DimPlot(merged.all, reduction="umap.rpca.pc40", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6099R" & has_dupli_chrX == T)),
        sizes.highlight = 0.5) +
  ggtitle("6099R, Duplication in chrX, Integrated")

DimPlot(merged.all, reduction="umap.rpca.pc40", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6146" & has_dupli_chr21 == T)),
        sizes.highlight = 0.5) +
  ggtitle("6146, Duplication in chr21, Integrated")

DimPlot(merged.all, reduction="umap.rpca.pc40", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6146" & has_dupli_chrX == T)),
        sizes.highlight = 0.5) +
  ggtitle("6146, Duplication in chrX, Integrated")

DimPlot(merged.all, reduction="umap", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6099R" & has_dupli_chr21 == T)),
        sizes.highlight = 0.5) +
  ggtitle("6099R, Duplication in chr21")

DimPlot(merged.all, reduction="umap", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6099R" & has_dupli_chrX == T)),
        sizes.highlight = 0.5) +
  ggtitle("6099R, Duplication in chrX")

DimPlot(merged.all, reduction="umap", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6146" & has_dupli_chr21 == T)),
        sizes.highlight = 0.5) +
  ggtitle("6146, Duplication in chr21")

DimPlot(merged.all, reduction="umap", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6146" & has_dupli_chrX == T)),
        sizes.highlight = 0.5) +
  ggtitle("6146, Duplication in chrX")

DimPlot(merged.all, reduction="umap", 
        cells.highlight = colnames(subset(merged.all, subset = Sample == "6146" & has_dupli_chrX == T)),
        sizes.highlight = 0.5) +
  ggtitle("6146, Duplication in chrX")

dev.off()

# Final annotation and ADT maps ----

## Removing all metadata columns from doublet finder 
merged.all@meta.data <- merged.all@meta.data[, !grepl(paste0("^", "DF"), colnames(merged.all@meta.data))]
merged.all@meta.data <- merged.all@meta.data[, !grepl(paste0("^", "pANN"), colnames(merged.all@meta.data))]
merged.all$Level1.anno.Stage <- NULL
merged.all$patient.sample.Stage <- NULL

merged.all@meta.data <- mutate(merged.all@meta.data, Level1.anno.Sample =
                                    case_when(Level1.anno == "Blast" ~ paste0(Level1.anno, "_", Sample),
                                              .default = as.character(Level1.anno)))

merged.all$celltype.short <- factor(merged.all$celltype.short,
                                                levels=c("CLP", "Pre-pro-B","Pro-B", "Pre-B", "Immature-B", "Mature-B",
                                                         "HSPC", "LMPP",
                                                         "GMP", "MEP", "DC-Progenitor", "pDC", "cDC","Mono",
                                                         "Pro-T", "Pre-T","CyclingDP", "Naive_T", "alpha-beta", "alpha-beta(mature)", "EffectorT",
                                                         "gamma-delta","NK"))

hspc.ctypes <- c("HSPC", "LMPP")
my.ctypes <- c("GMP", "MEP", "DC-Progenitor", "pDC", "cDC","Mono")
b.ctypes <- c("CLP", "Pre-pro-B","Pro-B", "Pre-B", "Immature-B", "Mature-B")
t.ctypes <- c("Pro-T", "Pre-T","CyclingDP", "Naive_T", "alpha-beta", "alpha-beta(mature)", "EffectorT",
              "gamma-delta","NK")

merged.all$Level1.anno.Sample <- factor(merged.all$Level1.anno.Sample, 
                                    levels=c("Blast_4526R", "Blast_5481", "Blast_5995", "Blast_6024", 
                                             "Blast_6092", "Blast_6061", "Blast_6099R", "Blast_6146", 
                                             "Blast_5850", "Blast_5907", "B", "Myeloid", "T_NK"))

blasts.all <- subset(merged.all, subset = Level1.anno == "Blast")

blast.freq.dt.colors <- c("#619ccf","#75a8d5","#89b5db","#9cc1e1","#b0cee7","#c4daed",
                          "#00ba38", "#60d483",
                          "#f8766d","#f9877f","#fa9892","#fba9a4","#fbbbb6","#fcccc8", 
                          "#703ac5","#8050cc", "#9066d2", "#a07cd9", "#b092df", "#bfa7e5", "#cfbdec", "#dfd3f2", "#efe9f9")

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/FrequencyPlots_PostInferCNVAnno_112025.pdf",
    height=6, width=9)
cluster.plot(merged.all@meta.data, Level1.anno.Sample, celltype.short, color.scale = blast.freq.dt.colors) + ggtitle("All Cells")
cluster.plot(blasts.all@meta.data, Level1.anno.Sample, celltype.short, color.scale = blast.freq.dt.colors) + ggtitle("Blasts")
dev.off()



h1 <- adt_heatmap(merged.all, groupByVar = "Level1.anno.Sample", clusterIdents = F, title="Median ADT, Blasts and Healthy Samples") 

h2 <- adt_heatmap(merged.all, groupByVar = "Blast_Lineage", clusterIdents = F, title="Median ADT, Blast Lineages and Healthy Samples") 

merged.all@meta.data <- mutate(merged.all@meta.data, Blast_Stage = case_when(
  Level1.anno=="Blast" & Stage=="IDX" ~ "Blast_IDX", 
  Level1.anno=="Blast" & Stage=="MRD"  ~ "Blast_MRD", 
  .default = Level1.anno
))

h3 <- adt_heatmap(merged.all, groupByVar = "Blast_Stage", clusterIdents = F, title="Median ADT, Blast Stages and Healthy Samples") 

infercnv.6061@meta.data <- mutate(infercnv.6061@meta.data, Blast_Lineage = case_when(
  Level1.anno=="Blast" & celltype.short %in% hspc.ctypes ~ "HSPC-like", 
  Level1.anno=="Blast" & celltype.short %in% my.ctypes ~ "Myeloid-like", 
  Level1.anno=="Blast" & celltype.short %in% b.ctypes ~ "B-like", 
  Level1.anno=="Blast" & celltype.short %in% t.ctypes ~ "T-cell-like", 
  .default = Level1.anno
))

infercnv.6061$Blast_Lineage <- factor(infercnv.6061$Blast_Lineage, 
                                   levels=c("B-like", "Myeloid-like", "HSPC-like", "T-cell-like", "B", "T_NK", "Myeloid"))


h4 <- adt_heatmap(infercnv.6061, groupByVar = "Blast_Lineage",
                  prots = c("CD45", "CD10", "CD19.1", "CD33.1", "CD34.1", "CD64",
                            "CD20", "CD3"),
                  clusterIdents = F)

infercnv.6146@meta.data <- mutate(infercnv.6146@meta.data, Blast_Lineage = case_when(
  Level1.anno=="Blast" & celltype.short %in% hspc.ctypes ~ "HSPC-like", 
  Level1.anno=="Blast" & celltype.short %in% my.ctypes ~ "Myeloid-like", 
  Level1.anno=="Blast" & celltype.short %in% b.ctypes ~ "B-like", 
  Level1.anno=="Blast" & celltype.short %in% t.ctypes ~ "T-cell-like", 
  .default = Level1.anno
))

infercnv.6146$Blast_Lineage <- factor(infercnv.6146$Blast_Lineage, 
                                      levels=c("B-like", "Myeloid-like", "HSPC-like", "T-cell-like", "B", "T_NK", "Myeloid"))

h5 <- adt_heatmap(infercnv.6146, groupByVar = "Blast_Lineage",
                  prots = c("CD45", "CD10", "CD19.1", "CD33.1","HLA-DR","CD38.1", "CD117", "CD72.1", 
                            "CD20", "CD3"),
                  clusterIdents = F)
h6<- adt_heatmap(pair3.s, groupByVar = "Blast_Lineage",
                 prots = c("CD20","CD24.1","CD22.1","CD33.1", "CD36.1",
                           "CD371", "CD117", "CD298", "CD45" ,"CD56", "CD11b",
                           "CD123", "CD64", 
                           "CD34.1", "CD72.1", "CD19.1"),
                 clusterIdents = F)

h7 <- adt_heatmap(pair4.s, groupByVar = "Blast_Lineage",
                  prots=c("CD45","CD3","CD38.1", "CD305", "CD117",
                          "CD36.1", "CD371", "CD20", "CD33.1", 
                          "CD34.1", "CD298", "CD13", 
                          "CD22.1","CD72.1","CD19.1", "CD10"),
               clusterIdents = F)

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/ADTPlots_PostInferCNVAnno_112025.pdf",
    height=6, width=6)
h1
h2
h3
h4
h5
h6
h7
dev.off()

saveRDS(merged.all, "/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/5Pairs_Integrated_inferCNVAnno_112025.rds")

merged.all <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/5Pairs_Integrated_inferCNVAnno_112025.rds")

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/ADTPlots_PostInferCNVAnno_Reduced_122525.pdf",
    height=6, width=6)
h1 <- adt_heatmap(merged.all, 
                  groupByVar = "Blast_Lineage", 
                  clusterIdents = F, title="Median ADT, Blast Lineages and Healthy Samples",
                  prots=c("CD33.1", "CD371", "CD36.1", "CD305",
                          "CD3", "CD20", "CD10", "CD24.1", "CD19.1", "CD45",
                          "CD72.1", "CD123", "CD64", "CD83.1", "CD13", "CD117", 
                          "CD38.1", "CD34.1")) 

merged.all$Blast_Stage <- factor(merged.all$Blast_Stage, levels=c("Blast_IDX", "Blast_MRD", 
                                                                  "B", "T_NK", "Myeloid"))

h2 <- adt_heatmap(merged.all, 
                  groupByVar = "Blast_Stage", 
                  clusterIdents = F, title="Median ADT, Blast Stages and Healthy Samples",
                  prots=c("CD33.1", "CD371", "CD36.1", "CD305",
                          "CD3", "CD20", "CD10", "CD24.1", "CD19.1", "CD45",
                          "CD72.1", "CD123", "CD64", "CD83.1", "CD13", "CD117", 
                          "CD38.1", "CD34.1")) 

merged.all$Blast_Lineage_Stage <- paste0(merged.all$Blast_Lineage, "_", merged.all$Stage)

merged.all$Blast_Lineage_Stage <- factor(merged.all$Blast_Lineage_Stage, 
                                         levels=c("B-like_IDX","Myeloid-like_IDX","HSPC-like_IDX","T-cell-like_IDX","B_IDX","T_NK_IDX","Myeloid_IDX",
                                                  "B-like_MRD","Myeloid-like_MRD","HSPC-like_MRD","T-cell-like_MRD","B_MRD","T_NK_MRD","Myeloid_MRD"))

h3 <- adt_heatmap(merged.all, 
                  groupByVar = "Blast_Lineage_Stage", 
                  clusterIdents = F, title="Median ADT, Blast Lineages by Stage and Healthy Samples",
                  prots=c("CD33.1", "CD371", "CD36.1", "CD305",
                          "CD3", "CD20", "CD10", "CD24.1", "CD19.1", "CD45",
                          "CD72.1", "CD123", "CD64", "CD83.1", "CD13", "CD117", 
                          "CD38.1", "CD34.1"))

merged.idx <- subset(merged.all, subset = Stage == "IDX")

h4 <- adt_heatmap(merged.idx, 
                  groupByVar = "Blast_Lineage_Stage", 
                  clusterIdents = F, title="Median ADT, Blast Lineages IDX and Healthy Samples",
                  prots=c("CD33.1", "CD371", "CD36.1", "CD305",
                          "CD3", "CD20", "CD10", "CD24.1", "CD19.1", "CD45",
                          "CD72.1", "CD123", "CD64", "CD83.1", "CD13", "CD117", 
                          "CD38.1", "CD34.1"))

merged.mrd <- subset(merged.all, subset = Stage == "MRD")

h5 <- adt_heatmap(merged.mrd, 
                  groupByVar = "Blast_Lineage_Stage", 
                  clusterIdents = F, title="Median ADT, Blast Lineages MRD and Healthy Samples",
                  prots=c("CD33.1", "CD371", "CD36.1", "CD305",
                          "CD3", "CD20", "CD10", "CD24.1", "CD19.1", "CD45",
                          "CD72.1", "CD123", "CD64", "CD83.1", "CD13", "CD117", 
                          "CD38.1", "CD34.1"))

dev.off()

# Plot individual ADT Maps ---

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/ADTPlots_PostInferCNVAnno_BySample_112125.pdf",
    height=8, width=5)
#merged.samples <- SplitObject(merged.all, "Sample")

for(i in 1:length(merged.samples)){
  adt_heatmap(merged.samples[[i]], 
              groupByVar = "Level1.anno.Sample", 
              clusterIdents = F, 
              title=paste0("Median ADT, ", names(merged.samples)[i], " Level1 anno"))
  
  adt_heatmap(merged.samples[[i]], 
              groupByVar = "Blast_Lineage", 
              clusterIdents = F, 
              title=paste0("Median ADT, ", names(merged.samples)[i], " Lineage"))
}

dev.off()


## Plotting supplementary UMAPs ---

merged.all <- readRDS("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Objects/5Pairs_Integrated_inferCNVAnno_112025.rds")

merged.all@meta.data <- mutate(merged.all@meta.data, Subtype = case_when(
  patient.sample %in% c("PT3", "PT5") ~ "KMT2Ar", 
  .default = "Ph-like"
))

merged.all@meta.data <- mutate(merged.all@meta.data, Patient.Stage.Healthy = case_when(
  Level1.anno == "Blast" ~ Patient.Stage.All, 
  .default = Level1.anno
))


pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/inferCNV_5Pairs_SamplesInFullCohort_AdditionalUMAPs.pdf",
    height=5, width=7)

DimPlot(merged.all, group.by="Stage")
DimPlot(merged.all, group.by="patient.sample")
DimPlot(merged.all, group.by="Subtype")
DimPlot(merged.all, group.by="Patient.Stage.All")

merged.all$Patient.Stage.Healthy <- factor(merged.all$Patient.Stage.Healthy, 
                                           levels=c("B", "Myeloid", "T_NK", "PT1_IDX", "PT1_MRD",
                                                    "PT2_IDX", "PT2_MRD","PT3_IDX", "PT3_MRD",
                                                    "PT4_IDX", "PT4_MRD", "PT5_IDX", "PT5_MRD"))

DimPlot(merged.all, group.by="Patient.Stage.Healthy",
        cols=c("#944939","#7F7A3A","#1B4A4A",
               "#f25e60","#e41a1c",
               "#6ea8d6","#377eb8",
               "#7fd977","#4daf4a",
               "#c285cd","#984ea3",
               "#ffb54d","#ff7f00"))

dev.off()

merged.blasts <- subset(merged.all, subset = Level1.anno == "Blast")
merged.blasts$Blast_Lineage_Stage <- paste0(merged.blasts$Blast_Lineage, "_", merged.blasts$Stage)

RidgePlot(merged.blasts, features=c("CD305"), assay="ADT", group.by = "Blast_Stage")
p1 <- RidgePlot(merged.blasts, features=c("CD305"), assay="ADT", group.by = "Blast_Lineage_Stage")

ggsave("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/CD305_RidgePlot_020526.pdf",
       p1)

## Volcano Plots for marker differences 

merged.all <- JoinLayers(merged.all, assay="ADT")
merged.all$Blast_Lineage_Stage <- paste0(merged.all$Blast_Lineage, "_", merged.all$Stage)

myl.vs.b.idx <- FindMarkers(merged.all, ident.1 = "Myeloid-like_IDX", ident.2 = "B-like_IDX",
                            assay = "ADT", group.by = "Blast_Lineage_Stage",
                            logfc.threshold = 0)

myl.vs.b.mrd <- FindMarkers(merged.all, ident.1 = "Myeloid-like_MRD", ident.2 = "B-like_MRD",
                            assay = "ADT", group.by = "Blast_Lineage_Stage",
                            logfc.threshold = 0)

mrd.vs.idx.b <- FindMarkers(merged.all, ident.1 = "B-like_MRD", ident.2 = "B-like_IDX",
                            assay = "ADT", group.by = "Blast_Lineage_Stage",
                            logfc.threshold = 0)

mrd.vs.idx.myl <- FindMarkers(merged.all, ident.1 = "Myeloid-like_MRD", ident.2 = "Myeloid-like_IDX",
                              assay = "ADT", group.by = "Blast_Lineage_Stage",
                              logfc.threshold = 0)

library(EnhancedVolcano)

EnhancedVolcano(
  myl.vs.b.idx,
  lab = rownames(myl.vs.b.idx),
  x   = "avg_log2FC",   # or "log2FC" in newer Seurat
  y   = "p_val_adj",    # adjusted p values
  title    = "Myeloid vs. B (IDX)",
  xlab     = "log2 fold change",
  ylab     = "-log10 adjusted p-value",
  pCutoff  = 0.05,
  FCcutoff = 0.25
)

markers.list <- list(myl.vs.b.idx, myl.vs.b.mrd, mrd.vs.idx.b, mrd.vs.idx.myl)
titles.list <- c("Myeloid vs. B IDX", "Myeloid vs. B MRD", "MRD vs. IDX B", "MRD vs. IDX Myl")

plot.list <- list()

pdf("/mnt/isilon/tan_lab/mehtar1/MRD_RM/CITE-Seq/Plots/Sample_Plots_Round3/MarkerQuantification_VolcanoPlots_020626.pdf",
    height=12, width=10)
for(i in 1:length(markers.list)){
  
  # markers: output of FindMarkers with columns 'avg_log2FC' and 'p_val_adj'
  # define significance cutoffs
  p_cut  <- 0.05
  fc_cut <- 0.25
  
  markers <- markers.list[[i]]
  
  # build a named colour vector: red if significant & logFC>0, blue if significant & logFC<0, grey otherwise
  col_vec <- ifelse(
    markers$p_val_adj < p_cut & markers$avg_log2FC >  fc_cut, "red",
    ifelse(
      markers$p_val_adj < p_cut & markers$avg_log2FC < -fc_cut, "blue",
      "grey"
    )
  )
  
  # Initialize all as NS
  names(col_vec) <- "NS"
  
  # Mark upregulated (red) as Up
  names(col_vec)[col_vec == "red"]  <- "Up"
  
  # Mark downregulated (blue) as Down
  names(col_vec)[col_vec == "blue"] <- "Down"
  # names(col_vec) <- rownames(markers)
  
  p1 <- EnhancedVolcano(
    markers,
    lab        = rownames(markers),
    x          = "avg_log2FC",
    y          = "p_val_adj",
    title      = titles.list[i],
    xlab       = "log2 fold change",
    ylab       = "-log10 adjusted p-value",
    pCutoff    = p_cut,
    FCcutoff   = fc_cut,
    colCustom  = col_vec,      # use your custom colours
    colAlpha   = 0.4,            # make colours solid (optional)
    legendLabels = c("NS", "NS", "Up", "Down"),  # optional, for legend clarity
    # turn off gridlines
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    # draw a box around the whole panel
    border       = "full",   # "partial" would only draw axes
    borderWidth  = 1.0,
    borderColour = "black"  
    )
  
  plot.list[[i]] <- p1
  
}

print(ggarrange(plotlist=plot.list, ncol=2, nrow = 2))

dev.off()

