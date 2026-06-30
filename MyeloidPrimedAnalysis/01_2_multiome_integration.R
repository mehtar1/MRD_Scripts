library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(tidyverse)
library(harmony)
library(chromVAR)
# Link gene and peaks
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(scATAC.Peak.annotation)
library(pheatmap)


##### 1. Patient.Samples in k_4_Myeloid C1&C2

df_joint <- read.csv("/PATH/HTAN_MRD_AllSamples_ATACproj2healthyCells_to_MultiomeRef_UMAPlsi25_Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels_Aug21_2023_k4cluster_092024.csv",sep=",")
table(df_joint$k_4_Myeloid)

df_k_4_Myeloid_C1_2 <- df_joint[which(df_joint$k_4_Myeloid %in% c(1,2)),]
df_k_4_Myeloid_C1_2_sample <- unique(df_k_4_Myeloid_C1_2$Patient.Sample)



##### 2. load C1&C2 samples and select blast cells

sample_list = unique(df_k_4_Myeloid_C1_2$Patient.Sample)
Object_list = list()
#  
for (sample in sample_list){
       print(sample)   
    # get the file directory
      input_Objfiledir = str_c("/PATH/JointAssay_SeuratObj/",sample)
    # read the Seurat object
      Joint.object <- readRDS(str_c(input_Objfiledir,"_JointSeuratObj_doubletRemoved.rds"))
      DefaultAssay(Joint.object) <- "ATAC"

    # remove fragment information from assay
      Fragments(Joint.object) <- NULL  
    # Constructing the Fragment class and adding it to the Seurat object
      input_Fragfiledir = "/PATH/fragments/persample/new/"
      frag.path <- str_c(input_Fragfiledir,sample,".bed.gz")
      fragments <- CreateFragmentObject(
        path = frag.path,
        cells = colnames(Joint.object), 
        validate.fragments = TRUE)
      Fragments(Joint.object) <- fragments
      print(sample)

    # put the Seurat object into the list
      Object_list[[sample]] = Joint.object
    }
#
length(Object_list)

# Merge the Seurat object
BALL.Object_1 <- merge(x = Object_list[[1]], y = Object_list[2:25])
BALL.Object_2 <- merge(x = Object_list[[26]], y = Object_list[27:length(Object_list)])


# add the metadata
# BALL_1
df_BALL_1_meta <- BALL.Object_1@meta.data
identical(rownames(df_BALL_1_meta), df_BALL_1_meta$rn)

intersect_col <- intersect(colnames(df_joint), colnames(df_BALL_1_meta))
df_BALL_1_meta <- df_BALL_1_meta[,!colnames(df_BALL_1_meta) %in% intersect_col]

DF_meta_1 <- merge(x = df_BALL_1_meta, y = df_joint, by.x = "rn", by.y = "patient_bc", all.x = TRUE)
identical(rownames(BALL.Object_1@meta.data), DF_meta_1$rn)

rownames(DF_meta_1) <- DF_meta_1$rn
DF_meta_1 <- DF_meta_1[rownames(BALL.Object_1@meta.data),]
identical(rownames(BALL.Object_1@meta.data), DF_meta_1$rn)

BALL.Object_1@meta.data <- DF_meta_1


# BALL_2
df_BALL_2_meta <- BALL.Object_2@meta.data
identical(rownames(df_BALL_2_meta), df_BALL_2_meta$rn)

intersect_col <- intersect(colnames(df_joint), colnames(df_BALL_2_meta))
df_BALL_2_meta <- df_BALL_2_meta[,!colnames(df_BALL_2_meta) %in% intersect_col]

DF_meta_2 <- merge(x = df_BALL_2_meta, y = df_joint, by.x = "rn", by.y = "patient_bc", all.x = TRUE)
identical(rownames(BALL.Object_2@meta.data), DF_meta_2$rn)

rownames(DF_meta_2) <- DF_meta_2$rn
DF_meta_2 <- DF_meta_2[rownames(BALL.Object_2@meta.data),]
identical(rownames(BALL.Object_2@meta.data), DF_meta_2$rn)

BALL.Object_2@meta.data <- DF_meta_2