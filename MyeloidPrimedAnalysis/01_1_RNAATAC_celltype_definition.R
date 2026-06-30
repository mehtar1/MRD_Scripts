library(ggplot2)
library(tidyverse)

DF <- read.csv("/PATH/HTAN_MRD_AllSamples_ATACproj2healthyCells_to_MultiomeRef_UMAPlsi25_Seurat_AllDonors_DownsampleBtraj15000_JointAnalysis_RNAGini_ATACDAP_RefineLabels_Aug21_2023_k4cluster_092024.csv", sep = ",")

### 1. select samples 
DF <- DF[!is.na(DF$k_4_Myeloid),]
DF <- DF[which(DF$k_4_Myeloid %in% c("1","2")),]

### 2. subset blast cell and filter out T/NK cells
DF <- DF[which(DF$Level1.anno %in% c("Blast")),] # 238380

### 3. define Bcell & Myeloid lineage committed cells
## 3.1.1 define RNA combined cell type
DF$Projected.Ctype.Wmethod_bin <- DF$Projected.Ctype.Wmethod
DF$Projected.Ctype.Wmethod_bin <- ifelse(
    DF$Projected.Ctype.Wmethod_bin %in% c("HSC/MPP","LMPP"), "HSPC", 
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("CLP"), "CLP", 
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("Pre-pro-B"), "Pre-pro-B", 
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("Pro-B","Pre-B"), "Pro-B-Pre-B",
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("Immature-B","Mature-B-1","Mature-B-2"), "Immature-Mature-B",
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("GMP","DC-progenitor","MEP"), "Myeloid-Progenitor",
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("Monocyte","cDC","pDC"), "Mature-Myeloid",
    ifelse(DF$Projected.Ctype.Wmethod_bin %in% c("T-cell"), "T-cell","NK"))))))))
table(DF$Projected.Ctype.Wmethod_bin)

# 3.1.2 define RNA Lineage
DF$Projected.Ctype.Wmethod_Lineage <- DF$Projected.Ctype.Wmethod
DF$Projected.Ctype.Wmethod_Lineage <- ifelse(
    DF$Projected.Ctype.Wmethod_Lineage %in% c("HSC/MPP","LMPP"), "HSPC", 
    ifelse(DF$Projected.Ctype.Wmethod_Lineage %in% c("CLP", "Pre-pro-B", "Pro-B","Pre-B", "Immature-B","Mature-B-1","Mature-B-2"), "B-Lineage-committed", 
    ifelse(DF$Projected.Ctype.Wmethod_Lineage %in% c("GMP","DC-progenitor","Monocyte","cDC","pDC","MEP"), "M-Lineage-committed","TNK")))
# filter TNK
DF <- DF %>% filter(Projected.Ctype.Wmethod_Lineage != "TNK")
table(DF$Projected.Ctype.Wmethod_Lineage)
DF$Projected.Ctype.Wmethod_Lineage <- factor(DF$Projected.Ctype.Wmethod_Lineage, levels = c("HSPC","B-Lineage-committed","M-Lineage-committed"))

# 3.2.1 define ATAC Lineage
DF$predicted.celltype.lsi_Lineage <- DF$predicted.celltype.lsi
DF$predicted.celltype.lsi_Lineage <- ifelse(
    DF$predicted.celltype.lsi_Lineage %in% c("HSC/MPP","LMPP"), "HSPC", 
    ifelse(DF$predicted.celltype.lsi_Lineage %in% c("CLP","Pre-pro-B","Pro-B","Pre-B","Immature-B","Mature-B-1","Mature-B-2"), "B-Lineage-committed",
    ifelse(DF$predicted.celltype.lsi_Lineage %in% c("GMP","DC-progenitor","Monocyte","cDC","pDC","MEP"), "M-Lineage-committed","TNK")))
table(DF$predicted.celltype.lsi_Lineage)
#  filter out TNK
DF <- DF[DF$predicted.celltype.lsi_Lineage != "TNK",]

# 3.2.2 add ATAC define to RNA define
DF$RNA.Lineage_ATAC.Lineage <- paste(DF$Projected.Ctype.Wmethod_Lineage, DF$predicted.celltype.lsi_Lineage, sep = "_")
DF$RNA.Lineage_ATAC.Lineage_bin <- DF$RNA.Lineage_ATAC.Lineage
DF$RNA.Lineage_ATAC.Lineage_bin <- ifelse(
    DF$RNA.Lineage_ATAC.Lineage_bin %in% c("B-Lineage-committed_B-Lineage-committed","B-Lineage-committed_HSPC"), "B-Lineage-committed_B-Lineage-committed",
    ifelse(DF$RNA.Lineage_ATAC.Lineage_bin %in% c("B-Lineage-committed_M-Lineage-committed"), "B-Lineage-committed_M-Lineage-committed",
    ifelse(DF$Projected.Ctype.Wmethod_Lineage %in% c("M-Lineage-committed"), "M-Lineage-committed",
    ifelse(DF$RNA.Lineage_ATAC.Lineage %in% c("HSPC_HSPC"), "HSPC_HSPC",
    ifelse(DF$RNA.Lineage_ATAC.Lineage %in% c("HSPC_B-Lineage-committed"), "HSPC_B-Lineage-committed",
    ifelse(DF$RNA.Lineage_ATAC.Lineage %in% c("HSPC_M-Lineage-committed"), "HSPC_M-Lineage-committed", "others"))))))
table(DF$RNA.Lineage_ATAC.Lineage_bin)

write.csv(DF, "/PATH/BALL_C1C2_BLAST_25PTs_55Samples_metadata_20250515.csv")