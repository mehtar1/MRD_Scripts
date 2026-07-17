library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(Matrix)

# =============================================================================
# 1. load Seurat obj
# =============================================================================
BALL.Obj <- readRDS("/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts.rds")




# =============================================================================
# 2. prepare differential data
# =============================================================================

## diff Motifs
Up_Motif <- c("JUN","JUNB","JUND","FOS","FOSB","FOSL2","BACH1","BACH2","CEBPA","CEBPG","ATF4","SPI1","BCL11A","MAFK","MAFG","MAFF","NFE2L2","NFE2L3","NFIL3","ELF1","ETV6","ETS2","HIVEP2","SMARCC1","CREB5","ATF7","PRDM1","STAT2")

################################################################
data("human_pwms_v2")
DF_chromvar <- BALL.Obj@assays$chromvar@data
DF_chromvar <- t(DF_chromvar) %>% as.data.frame()
colnames(DF_chromvar) <- sapply(strsplit(colnames(DF_chromvar),"-"), `[`, 3) %>% unlist()
identical(rownames(DF_chromvar),rownames(BALL.Obj@meta.data))
DF_chromvar$RNA.Lineage_ATAC.Lineage_bin <- BALL.Obj@meta.data$RNA.Lineage_ATAC.Lineage_bin
DF_chromvar_avg <- DF_chromvar[,c(Up_Motif,"RNA.Lineage_ATAC.Lineage_bin")]
DF_chromvar_avg <- reshape2::melt(DF_chromvar_avg, id.vars = "RNA.Lineage_ATAC.Lineage_bin")
DF_chromvar_avg <- DF_chromvar_avg %>% dplyr::group_by(RNA.Lineage_ATAC.Lineage_bin, variable) %>% dplyr::mutate(mean_value = mean(value))
DF_chromvar_avg <- DF_chromvar_avg[,c("RNA.Lineage_ATAC.Lineage_bin","variable","mean_value")] %>% unique() %>% as.data.frame()
DF_chromvar_avg <- reshape2::dcast(DF_chromvar_avg, variable ~ RNA.Lineage_ATAC.Lineage_bin, value.var = "mean_value")

# B_M - B_B 
DF_chromvar_avg$avg_diff <- DF_chromvar_avg[,c("B-Lineage-committed_M-Lineage-committed")] - DF_chromvar_avg[,c("B-Lineage-committed_B-Lineage-committed")] 
DF_chromvar_avg <- DF_chromvar_avg[,c("variable","avg_diff")]
DF_chromvar_avg$cluster <- "B_B_Down"
DF_chromvar_avg <- DF_chromvar_avg %>% dplyr::rename("gene"="variable")
saveRDS(DF_chromvar_avg, file = '/PATH/IDX_MyeloidPrim_BMvsBB_diff_chromVAR_dscore_by_atac_cell_states_allTFsKept_CISBP_AllUpMotif.rds')

# M - B_B 
DF_chromvar_avg$avg_diff <- DF_chromvar_avg[,c("M-Lineage-committed")] - DF_chromvar_avg[,c("B-Lineage-committed_B-Lineage-committed")] 
DF_chromvar_avg <- DF_chromvar_avg[,c("variable","avg_diff")]
DF_chromvar_avg$cluster <- "B_B_Down"
DF_chromvar_avg <- DF_chromvar_avg %>% dplyr::rename("gene"="variable")
saveRDS(DF_chromvar_avg, file = '/PATH/IDX_MyeloidPrim_MvsBB_diff_chromVAR_dscore_by_atac_cell_states_allTFsKept_CISBP_AllUpMotif.rds')

# directly input the diff.TFs
diff.tf_BMvsBB = readRDS(file = '/PATH/IDX_MyeloidPrim_BMvsBB_diff_chromVAR_dscore_by_atac_cell_states_allTFsKept_CISBP_AllUpMotif.rds')
diff.tf_MvsBB = readRDS(file = '/PATH/IDX_MyeloidPrim_MvsBB_diff_chromVAR_dscore_by_atac_cell_states_allTFsKept_CISBP_AllUpMotif.rds')
diff.tf <- rbind(diff.tf_BMvsBB, diff.tf_MvsBB)


############################################################
## degs
BMvsBB.DEG <- read.csv("/PATH/DEG_IDX_MyePrim_BMvsBB.csv",row.names="X")
BMvsBB.DEG <- BMvsBB.DEG[BMvsBB.DEG$p_val_adj < 0.05,]
BMvsBB.DEG <- BMvsBB.DEG[abs(BMvsBB.DEG$avg_log2FC) > 0.25,]
BMvsBB.DEG_Up <- BMvsBB.DEG[which(BMvsBB.DEG$FC_change == "Up"),]
BMvsBB.DEG_Up <- BMvsBB.DEG_Up[BMvsBB.DEG_Up$pct.1 > 0.1,]
BMvsBB.DEG_Down <- BMvsBB.DEG[which(BMvsBB.DEG$FC_change == "Down"),]
BMvsBB.DEG_Down <- BMvsBB.DEG_Down[BMvsBB.DEG_Down$pct.2 > 0.1,]
BMvsBB.DEG <- rbind(BMvsBB.DEG_Up, BMvsBB.DEG_Down)
BMvsBB.DEG$gene <- rownames(BMvsBB.DEG)


MvsBB.DEG <- read.csv("/PATH/DEG_IDX_MyePrim_MvsBB.csv",row.names="X")
MvsBB.DEG <- MvsBB.DEG[MvsBB.DEG$p_val_adj < 0.05,]
MvsBB.DEG <- MvsBB.DEG[abs(MvsBB.DEG$avg_log2FC) > 0.25,]
MvsBB.DEG_Up <- MvsBB.DEG[which(MvsBB.DEG$FC_change == "Up"),]
MvsBB.DEG_Up <- MvsBB.DEG_Up[MvsBB.DEG_Up$pct.1 > 0.1,]
MvsBB.DEG_Down <- MvsBB.DEG[which(MvsBB.DEG$FC_change == "Down"),]
MvsBB.DEG_Down <- MvsBB.DEG_Down[MvsBB.DEG_Down$pct.2 > 0.1,]
MvsBB.DEG <- rbind(MvsBB.DEG_Up, MvsBB.DEG_Down)
MvsBB.DEG$gene <- rownames(MvsBB.DEG)

DEG <- rbind(BMvsBB.DEG, MvsBB.DEG)


############################################################
## daps or not
# DAP
BMvsBB.DAP <- read.csv(file = "/PATH/DAP_IDX_MyePrim_BMvsBB.csv",row.names="X")
BMvsBB.DAP <- BMvsBB.DAP[BMvsBB.DAP$p_val_adj < 0.05,]
BMvsBB.DAP <- BMvsBB.DAP[abs(BMvsBB.DAP$avg_log2FC) > 0.25,]
BMvsBB.DAP_Up <- BMvsBB.DAP[which(BMvsBB.DAP$FC_change == "Up"),]
BMvsBB.DAP_Up <- BMvsBB.DAP_Up[BMvsBB.DAP_Up$pct.1 > 0.1,]
BMvsBB.DAP_Down <- BMvsBB.DAP[which(BMvsBB.DAP$FC_change == "Down"),]
BMvsBB.DAP_Down <- BMvsBB.DAP_Down[BMvsBB.DAP_Down$pct.2 > 0.1,]
BMvsBB.DAP <- rbind(BMvsBB.DAP_Up, BMvsBB.DAP_Down)
BMvsBB.DAP$Peak <- rownames(BMvsBB.DAP)
table(BMvsBB.DAP$FC_change)

MvsBB.DAP <- read.csv(file = "/PATH/DAP_IDX_MyePrim_MvsBB.csv",row.names="X")
MvsBB.DAP <- MvsBB.DAP[MvsBB.DAP$p_val_adj < 0.05,]
MvsBB.DAP <- MvsBB.DAP[abs(MvsBB.DAP$avg_log2FC) > 0.25,]
MvsBB.DAP_Up <- MvsBB.DAP[which(MvsBB.DAP$FC_change == "Up"),]
MvsBB.DAP_Up <- MvsBB.DAP_Up[MvsBB.DAP_Up$pct.1 > 0.1,]
MvsBB.DAP_Down <- MvsBB.DAP[which(MvsBB.DAP$FC_change == "Down"),]
MvsBB.DAP_Down <- MvsBB.DAP_Down[MvsBB.DAP_Down$pct.2 > 0.1,]
MvsBB.DAP <- rbind(MvsBB.DAP_Up, MvsBB.DAP_Down)
MvsBB.DAP$Peak <- rownames(MvsBB.DAP)
table(MvsBB.DAP$FC_change)


DAP <- rbind(BMvsBB.DAP, MvsBB.DAP)

# 2. filter weak peaks
weak.pks = readRDS('/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_0.2weak_peaks.rds')
open.frac.pks = readRDS('/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_peak_openFrac_byState.rds')


############################################################
## filter/add TF by gene expression
expr.frac.state = readRDS('/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_frac_expr_state.rds')
expr.median.state = readRDS('/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_median_expr_state.rds')
expr.avg.state = readRDS('/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_avg_expr_state.rds')
# tfs.deg = readRDS('MetaData/tfs_deg_state.rds')

##motif scan
motif_pk_match = readRDS('/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_motif_match_pv5E-05_CISBP.rds')
library(SingleCellExperiment)
motif_pk_match_count = assay(motif_pk_match)
colnames(motif_pk_match_count) <- motif_pk_match@colData$name
rownames(motif_pk_match_count) <- paste0(motif_pk_match@rowRanges@seqnames, '-', motif_pk_match@rowRanges@ranges)
motif_pk_match_count = as.matrix(motif_pk_match_count)




# =============================================================================
# 3. controlling parameters
# =============================================================================
filtered_links = NULL
filterWeakPeak = FALSE
useDAP = TRUE ## filter peaks by daps/foldChange or by peak open fraction(>0.05)
useAllDEGs = TRUE ## use all DEGs or filter
useEnhancerPeakOnly = TRUE ## use promoter enhancer interaction only or include PPIs
peakFracCutoff = 0.2 ## used to filter peaks per state (default 0.1)
tfexprfracCutoff = 0.1

if(useEnhancerPeakOnly){
  predicted_links = fread('/PATH/new_regrRes4_gene_peak_ep_links_metacell_IDX_MyeloidPrim_BB_BM_Mye.tsv')
}

# if(filterWeakPeak) predicted_links = predicted_links[peak_name %notin% weak.pks]
DAP.Up <- DAP[which(DAP$FC_change == "Up"),]
if(useDAP) predicted_links = predicted_links[peak_name %in% unique(DAP.Up$Peak)]




# =============================================================================
# 4. construct TRN per state
# =============================================================================

edges_list = nodes_list = list()
## provide or read your interested gene and tf lists

## Genes
DEG.Up <- DEG[which(DEG$FC_change == "Up"),]
DEG.Down <- DEG[which(DEG$FC_change == "Down"),]

interested_genes = unique(DEG.Up$gene)
open.frac.pks0 = open.frac.pks[, c("B-Lineage-committed_M-Lineage-committed","M-Lineage-committed")]

## as TF strength
avg.zscore.state = subset(diff.tf_BMvsBB, select = c('gene', 'avg_diff', 'cluster'))
avg.zscore.state0 = avg.zscore.state[which(avg.zscore.state$cluster == "B_B_Down"),]
avg.zscore.state0 = as.data.table(avg.zscore.state0)
setkey(avg.zscore.state0, gene)


interested_tfs = unique(Up_Motif)
interested_tfs = interested_tfs[interested_tfs %notin% unique(DEG.Down$gene)] 


## filtering peaks
filtered.peaks = unique(DAP.Up$Peak) # 13122
# if(!useDAP) filtered.peaks = names(which(open.frac.pks0 > peakFracCutoff)) # 0.2 : 18193 peaks


## further filtering TF by expression
expr.frac.state0 <- expr.frac.state[, c("B-Lineage-committed_M-Lineage-committed","M-Lineage-committed")]
expr.median.state0 <- expr.median.state[, c("B-Lineage-committed_M-Lineage-committed","M-Lineage-committed")]
expr.avg.state0 <- expr.avg.state[, c("B-Lineage-committed_M-Lineage-committed","M-Lineage-committed")]

tf.expr.frac = expr.frac.state0[interested_tfs,]
# tf.expr.frac = tf.expr.frac[!is.na(tf.expr.frac)]
valid_rows <- apply(tf.expr.frac, 1, function(row) all(row > tfexprfracCutoff))
tf.expr.frac <- tf.expr.frac[valid_rows,]
interested_tfs = rownames(tf.expr.frac)

# predicted links
predicted_links0 = predicted_links[gene_name %in% interested_genes]
sele.tf.mat = motif_pk_match_count[, interested_tfs]
sele.peaks = names(which(rowSums(sele.tf.mat > 0) > 0))
sele.peaks = intersect(sele.peaks, filtered.peaks) #  DAP : 11177 peaks
predicted_links0 = predicted_links0[peak_name %in% sele.peaks] # DAP : 4664 links


## save filtered links
predicted_links0 = predicted_links0[Estimate > 0.2 & fdr < 0.05] # 1048 links
predicted_links0$cell_state = "B_M&Myeloid"
filtered_links = rbind(filtered_links, predicted_links0) # DAP: 1048 links


## split links by tf
tf_regulons = list()
cell_state0 = "B-Lineage-committed_M-Lineage-committed"
open.frac.pks0 = open.frac.pks[, cell_state0]
for(TF0 in interested_tfs){
    peaks0 = names(which(sele.tf.mat[, TF0] > 0))
    ep0 = predicted_links0[peak_name %in% peaks0]
    ep0 = subset(ep0, select = c('gene_name', 'Estimate',  'fdr', 'peak_name'))
    ep0$peak_strength = open.frac.pks0[ep0$peak_name]
    ep0 = ep0[peak_strength > 0.05]
    ep0[, 'score' := -log10(fdr)]
    ep0$TF = TF0
    ep0[, 'edge_strength' := sum(Estimate * peak_strength), by = gene_name]
    ep0[, 'Estimate' := sum(Estimate), by = gene_name]
    ep0[, 'score' := mean(score), by = gene_name]
    ep0[, 'peaks' := paste0(peak_name, collapse = ','), by = gene_name]
    ep0[, 'peak_strengths' := paste0(peak_strength, collapse = ','), by = gene_name]
    ep0[, 'npeaks' := .N, by = gene_name]
    ep0 = subset(ep0, select = c('TF', 'gene_name', 'edge_strength', 'peak_strengths','Estimate', 'score', 'peaks', 'npeaks'))
    ## node inf
    ep0$TF_strength = avg.zscore.state0[ep0$TF]$avg_diff
    ep0$gene_strength = expr.avg.state0[ep0$gene_name]

    tf_regulons[[TF0]] = ep0[!duplicated(ep0)]
}

length(tf_regulons)


tf_regulons.comb = do.call('rbind', tf_regulons)
tf_regulons.comb = tf_regulons.comb[order(-edge_strength)]
message(paste0('B_M', ': ', nrow(tf_regulons.comb), ' regulons'))
  
nodes1 = subset(tf_regulons.comb, select = c('TF', 'TF_strength'))
nodes2 = subset(tf_regulons.comb, select = c('gene_name', 'gene_strength'))
names(nodes1) = names(nodes2) = c('gene', 'node_strength')
nodes1$node_strength = nodes1$node_strength/max(nodes1$node_strength)
nodes2$node_strength = nodes2$node_strength/max(nodes2$node_strength)
nodes1$node_type = 'TF'
nodes2$node_type = 'gene'
nodes2 = nodes2[gene %notin% nodes1$gene]
nodes = rbind(nodes1, nodes2) %>% unique()


cell_state0 = "B_M_Myeloid"
outdir = "/PATH/"
Peak_sele = "DAP"
#
tf_regulons.comb$cell_state = cell_state0
nodes$cell_state = cell_state0
 
tf_regulon_fname = paste0(outdir,'edges_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP.tsv')
write.table(tf_regulons.comb, file = tf_regulon_fname, sep = '\t', row.names = F, quote = F)

# 
nodes_fname = paste0(outdir,'nodes_', cell_state0, "_", Peak_sele,  '_pv5E-05_CISBP.tsv')
write.table(nodes, file = nodes_fname, sep = '\t', row.names = F, quote = F)

#
edges_list[[cell_state0]] = tf_regulons.comb
nodes_list[[cell_state0]] = nodes
# plot tf ranking
pdata = data.table('TF' = names(sort(table(tf_regulons.comb$TF))),
                    'nTargets' = as.numeric(sort(table(tf_regulons.comb$TF))))
pdata$TF = factor(pdata$TF, levels = names(sort(table(tf_regulons.comb$TF))))
p1 <- ggplot(data = pdata, aes(x = TF, y = nTargets, fill = nTargets)) + 
      geom_bar(stat = 'identity') +
      coord_flip() +
      theme_classic() + 
      scale_fill_viridis_c(option = "magma", begin = 0.2, direction = -1) +
      ggtitle(cell_state0) + xlab('')
tf_rank_fname = paste0(outdir,'ConnectRank_TF_regulons_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP.pdf')
# output
ggsave(p1, filename = tf_rank_fname, device = 'pdf', width = 5, height = 2 * ceiling(nrow(pdata)/20))


filtered_ep_file = paste0(outdir,'filtered_regrRes4_gene_peak_links_metacell_', cell_state0,"_", Peak_sele, '_ep_pv5E-05_CISBP.tsv')
fwrite(filtered_links, file = filtered_ep_file, sep = '\t')

edge_file = paste0(outdir,'TRN_edges_', cell_state0, "_", Peak_sele,'_pv5E-05_CISBP.xlsx')
node_file = paste0(outdir,'TRN_node_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP.xlsx')
writexl::write_xlsx(edges_list,path = edge_file)
writexl::write_xlsx(nodes_list,path = node_file)

edges_combined = do.call('rbind', edges_list)
nodes_combined = do.call('rbind', nodes_list)
edges_combined_file = paste0(outdir,'TRN_edges_', cell_state0, "_", Peak_sele, '_combined_pv5E-05_CISBP.xlsx')
nodes_combined_file = paste0(outdir,'TRN_nodes_', cell_state0, "_", Peak_sele, '_combined_pv5E-05_CISBP.xlsx')
writexl::write_xlsx(edges_combined, path = edges_combined_file)
writexl::write_xlsx(nodes_combined, path = nodes_combined_file)


# subset TF-Target linked for targets GO enrichment
df_edge <- read.csv(paste0(outdir,'edges_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP.tsv'), sep = "\t")
df_edge <- df_edge[,c("TF","gene_name")]
write.csv(df_edge, file = paste0(outdir,'edges_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP_Netw_forGO.csv'), row.names = F)

#
df <- read.csv(paste0(outdir,'edges_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP_Netw_forGO.csv'))
df$gene <- df$gene_name
df_TF <- pivot_wider(df, names_from = TF, values_from = gene_name)
df_TF <- df_TF %>% dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "", .))) %>% dplyr::select(-gene)
write.csv(df_TF, file = paste0(outdir,'edges_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP_Netw_forGO.csv'), row.names = F)


## prepare files for Cytoscape
# 1. target->type
df_edge <- read.csv(paste0(outdir,'edges_', cell_state0, "_", Peak_sele, '_pv5E-05_CISBP.tsv'), sep = "\t")
df1 <- df_edge[,c("TF","cell_state")]
df1$type <- "TF"
colnames(df1) <- c("gene","cell_state","type")
df2 <- df_edge[,c("gene_name","cell_state")]
df2$type <- "Target"
colnames(df2) <- c("gene","cell_state","type")
df3 <- rbind(df1,df2) %>% dplyr::distinct() %>% dplyr::select(-cell_state)
write.csv(df3, file = paste0(outdir,'edges_', cell_state0, '_', Peak_sele, '_pv5E-05_CISBP_Netw_TF_Target_type.csv'), row.names = F)

# 2. TF->Target->avg.LogFC
DEG <- read.csv("/PATH/DEG_IDX_MyePrim_MyeloidvsBM_NoFilter.csv", row.names = "X")
DEG$gene <- rownames(DEG)
df_merge <- merge(df_edge, DEG, by.x = "gene_name", by.y = "gene", all.x = T) %>% dplyr::distinct()
df_merge$FC_Type <- paste(df_merge$Versus, df_merge$FC, sep = "_")
df_avg <- df_merge[,c("TF","gene_name","avg_log2FC","FC_Type","p_val_adj")]
colnames(df_avg) <- c("TF","Target","Target.MvsBM.avg.LogFC","FC_Type","p_val_adj")
write.csv(df_avg, file =  paste0(outdir,'edges_', cell_state0, '_', Peak_sele, '_pv5E-05_CISBP_Netw_TF_Target_type_avgLogFC.csv'), row.names = F)