library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(tidyverse)
library(harmony)
library(pheatmap)
library(gridExtra)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(Matrix)


BALL.Obj <- readRDS("/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts.rds")
# RUN PCA
DefaultAssay(BALL.Obj) <- "RNA"
BALL.Obj <- NormalizeData(BALL.Obj) %>% FindVariableFeatures(., nfeatures = 2000) %>% ScaleData(.) %>% RunPCA(.,npc = 30, verbose = F)
saveRDS(BALL.Obj, "/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts.rds")


# =============================================================================
# 1. call metacell on joint obj
# =============================================================================

DefaultAssay(BALL.Obj) <- "RNA"
BALL.Obj <- SetupForWGCNA(
  BALL.Obj,
  gene_select = "fraction", # the gene selection approach
  #fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "Mid12PTsAndBot6PTs_BB_BM_Mye_hdWGCNA", # the name of the hdWGCNA experiment
  # features = genes4anchors
)

# construct metacells  in each group
BALL.Obj <- MetacellsByGroups(
  seurat_obj = BALL.Obj,
  group.by = c("RNA.Lineage_ATAC.Lineage_bin"), # specify the columns in seurat_obj@meta.data to group by
  k = 20, # nearest-neighbors parameter
  max_shared = 3, # maximum number of shared cells between two metacells
  min_cells = 100,
  reduction = 'pca',
  ident.group = 'RNA.Lineage_ATAC.Lineage_bin' # set the Idents of the metacell seurat object
)


# normalize metacell expression matrix:
BALL.Obj <- NormalizeMetacells(BALL.Obj)
plan("sequential") # It is not recommended to use parallelization when doing the WGCNA analysis.
metacell_obj <- GetMetacellObject(BALL.Obj)

if(ncol(metacell_obj) < 15){
  coembed <- MetacellsByGroups(
    seurat_obj = BALL.Obj,
    group.by = c("RNA.Lineage_ATAC.Lineage_bin"), # specify the columns in seurat_obj@meta.data to group by
    k = 20, # nearest-neighbors parameter
    max_shared = 3, # maximum number of shared cells between two metacells
    min_cells = 100, # ignore cluster with less than min_cells
    reduction = 'pca',
    ident.group = 'RNA.Lineage_ATAC.Lineage_bin'
  )
  # normalize metacell expression matrix:
  BALL.Obj <- NormalizeMetacells(BALL.Obj)
  metacell_obj <- GetMetacellObject(BALL.Obj)
  
}

DefaultAssay(BALL.Obj) <- "RNA"
BALL.Obj <- ScaleMetacells(BALL.Obj, features=VariableFeatures(BALL.Obj))
message(paste('#of metacells:', ncol(metacell_obj)))

if(ncol(metacell_obj) > 50){
  BALL.Obj <- RunPCAMetacells(BALL.Obj, features=VariableFeatures(BALL.Obj), npcs = 10)
  #coembed <- RunHarmonyMetacells(coembed, group.by.vars='Sample')
  BALL.Obj <- RunUMAPMetacells(BALL.Obj, reduction='pca', dims=1:10)
  
  #p5 <- DimPlotMetacells(BALL.Obj, group.by = 'RNA.Lineage_ATAC.Lineage_bin') + ggtitle("cluster_metacell")
  #p6 <- DimPlot(coembed, group.by = 'cell_state') + ggtitle('Coembedded')
  #ggsave(p0 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3), 
  #       width = 15, height = 12, 
  #       filename = paste0('data/intermediate/coembedBySample/figures/new_coembed_', sampleID0, '.pdf'), 
  #       device = 'pdf')
  
}else{
  #ggsave(p0 + p2 + p3 + p4 + plot_layout(ncol = 2, heights = c(2, 2)), 
  #       width = 15, height = 12, 
  #       filename = paste0('data/intermediate/coembedBySample/figures/new_coembed_', sampleID0, '.pdf'), 
  #       device = 'pdf')
}



message('Call metacells Done! Saving metadata matrices...')

## construct and save metacell rna and atac matrices 
merged_cells = metacell_obj$cells_merged
# nrna_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'RNA_'))
# natac_cells = sapply(merged_cells, function(x) stringr::str_count(x, pattern = 'ATAC_'))

write(merged_cells, file = '/PATH/metacells_names_IDX_MyeloidPrim_BB_BM_Mye.txt')


# sele.metacells = merged_cells[nrna_cells >= 5 & natac_cells >= 5] ## ignore metacells which are mostly comprised of single modality
# nrna_cells = nrna_cells[names(sele.metacells)]
# natac_cells = natac_cells[names(sele.metacells)]


## manually
mask <- sapply(names(merged_cells), function(x) colnames(BALL.Obj) %in%
                 unlist(strsplit(merged_cells[x], ',', fixed = T)))
rownames(mask) = colnames(BALL.Obj)

# rna.mask <- sapply(names(merged_cells), function(x) mask[, x]/merged_cells[x])
rna.mtx <- BALL.Obj[['RNA']]$data %*% mask # remember atac cells in RNA assay have 0 account

# atac.mask <- sapply(names(natac_cells), function(x) mask[, x]/natac_cells[x])
atac.mtx <- BALL.Obj[['ATAC']]$data %*% mask

# colnames(rna.mtx) = colnames(atac.mtx) = paste0(sampleID0, '_', names(sele.metacells))

saveRDS(rna.mtx, file = paste0('/PATH/rna_metacell_mtx_IDX_MyeloidPrim_BB_BM_Mye.rds'))
saveRDS(atac.mtx, file = paste0('/PATH/atac_metacell_mtx_IDX_MyeloidPrim_BB_BM_Mye.rds'))
message(paste('#of final metacells:', ncol(rna.mtx)))

message('All Done!')




# =============================================================================
# 2. Regress gene peak links by metacell
# =============================================================================

#### 2.1 load self-defined functions ####
`%notin%` = Negate(`%in%`)
library(parallel)
library(Signac)
library(Seurat)
library(data.table)
library(magrittr)

source('/PATH/regr_gene_peak_links_metacell_malignant.R')



#### 2.2 load matrix and prepare for downstream regression analysis ####
metacell.rna.mtx <- readRDS('/PATH/rna_metacell_mtx_IDX_MyeloidPrim_BB_BM_Mye.rds')
metacell.atac.mtx <- readRDS('/PATH/atac_metacell_mtx_IDX_MyeloidPrim_BB_BM_Mye.rds')

## filter the features of metacells ####
## filter peaks that accessible in less than 3% of all cell type
peaks.mean.ctype <- sapply(unique(BALL.Obj$RNA.Lineage_ATAC.Lineage_bin), function(x){rowMeans(BALL.Obj@assays$ATAC@counts[, BALL.Obj$RNA.Lineage_ATAC.Lineage_bin == x] > 0)})
rmax = apply(peaks.mean.ctype, 1,  max)
summary(rmax)
filtered.peaks = names(which(rmax > 0.03)) # 109259 peaks
cand.peaks = rownames(metacell.atac.mtx)
all(filtered.peaks %in% cand.peaks)
metacell.atac.mtx = metacell.atac.mtx[cand.peaks %in% filtered.peaks, ]

## filter genes by expression in less than 3% of all cell type
mtx = BALL.Obj[['RNA']]$counts
mtx.ctype = sapply(unique(BALL.Obj$RNA.Lineage_ATAC.Lineage_bin), function(x) rowMeans(BALL.Obj[['RNA']]$counts[, BALL.Obj$RNA.Lineage_ATAC.Lineage_bin == x] >0) )
mfreq_gene = qlcMatrix::rowMax(mtx.ctype)
sele.genes = intersect(rownames(mtx), rownames(mtx)[which(mfreq_gene > 0.03)]) # 13341 genes


## gene-peak-affinity map ####
## construct gene peak affinity binary matrix
#--new version, use longest transcripts
annotation = readRDS("/PATH/GRanges_hg38.rds")
tss_ann <- Signac:::CollapseToLongestTranscript(annotation)
tss_ann <- as.data.frame(resize(tss_ann, width = 1, fix = 'start'))
names(tss_ann)[1] = 'chr'
tss_ann$Tss = tss_ann$start
gene_ann = tss_ann = data.table(tss_ann)
  
gene2peak.map <- get_gene2peak_map(gene.list = data.table('gene_name' = sele.genes), peak_names = filtered.peaks, gene_ann = gene_ann, distal_dist = 5e5)

peaks.used = colnames(gene2peak.map)
sele.genes = sele.genes[sele.genes %in% rownames(gene2peak.map)]

metacell.rna.mtx = metacell.rna.mtx[sele.genes, ]
metacell.atac.mtx = metacell.atac.mtx[peaks.used, ]

save(metacell.rna.mtx, metacell.atac.mtx, gene2peak.map,file = '/PATH/rna_atac_metacell_matrices_IDX_MyeloidPrim_BB_BM_Mye.RData')



#### 2.3 do regression gene by gene ####
load('/PATH/rna_atac_metacell_matrices_IDX_MyeloidPrim_BB_BM_Mye.RData')
ls()
sele.genes = rownames(metacell.rna.mtx)
peaks.used = rownames(metacell.atac.mtx)

metacells.mtx = rbind(metacell.rna.mtx, metacell.atac.mtx)
metacells.df = data.frame(t(metacells.mtx))
pnames = names(metacells.df)[-c(1:length(sele.genes))]
pnames = gsub('.', '-', pnames, fixed = T)
names(metacells.df)[c(1:length(sele.genes))] = sele.genes
names(metacells.df)[-c(1:length(sele.genes))] = pnames

## do it in parallel
runParallel = TRUE 
stime = Sys.time()

if(runParallel){
  regress_function <- function(gene0) {
    peaks = names(which(gene2peak.map[gene0, ] == 1))
    data <- subset(metacells.df, select = c(gene0, peaks))
    names(data)[1] = 'y'
    model <- lm(y ~ -1 + ., data=data)
    return(coef(summary(model)))
  }
  # Use the mclapply function to perform the regressions in parallel
  ncores = detectCores() - 1
  regr_list <- mclapply(sele.genes,  FUN = regress_function, 
                        mc.cores = ncores)
  message(paste('#cores:', ncores))
  names(regr_list) = sele.genes
  
}else{
  curate_dataList = list()
  for(gene0 in sele.genes){
    peaks = names(which(gene2peak.map[gene0, ] == 1))
    data <- subset(metacells.df, select = c(gene0, peaks))
    names(data)[1] = 'y'
    curate_dataList[[gene0]] = data
  }
  
  regr_list = list()
  for(x in sort(sele.genes)){
    rdata <- curate_dataList[[x]]
    regr_list[[x]] <- coef(summary(lm(y ~ . -1, data = rdata)))
  }
}
etime = Sys.time()
etime-stime /PATH/rna_atac_metacell_matrices_IDX_MyeloidPrim_BB_BM_Mye.RData
saveRDS(regr_list, "/PATH/new_regrRes4ep_prediction_metacell_IDX_MyeloidPrim_BB_BM_Mye.rds")

## summarize/filter loops ####
regr_list = readRDS("/PATH/new_regrRes4ep_prediction_metacell_IDX_MyeloidPrim_BB_BM_Mye.rds")
regr.sum <- lapply(names(regr_list), function(t){
  x = regr_list[[t]]
  x = x[, c(1, 4), drop = F]
  x = data.frame(x)
  x = data.table(x, keep.rownames = T)
  x$gene_name = t
  return(x)
})

regr.sum = do.call('rbind', regr.sum)
names(regr.sum)[c(1, 3)] = c('peak_name', 'P_value')
regr.sum[, 'p_val_bonf' := pmin(1, P_value*nrow(regr.sum))]
regr.sum$fdr = p.adjust(regr.sum$P_value, method = 'fdr')

regr.sum$peak_name = sapply(regr.sum$peak_name, function(x) gsub('`', '', x, fixed = T))
predicted_links = regr.sum[P_value < 0.01 & abs(Estimate) > 0.05 & grepl(peak_name, pattern = '^chr')]

fwrite(predicted_links, file = '/PATH/new_regrRes4ep_prediction_metacell_IDX_MyeloidPrim_BB_BM_Mye.tsv',sep = '\t')

useEnhancerPeakOnly = T

if(useEnhancerPeakOnly){
  ##older version
  if(F){
    # gene_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/gene_ann_hg38.txt')
    # gene_ann[, 'Tss' := ifelse(strand == '+', start, end)]
    
    # gene_ann_sele = gene_ann[gene_name %in% predicted_links$gene_name, ]
    
   
    # ## avoid peaks to overlap with promoters -- this step is optional 
    # #-- all transcripts
    # tss_ann = fread('/mnt/isilon/tan_lab/yuw1/R_work_dir/MLLr/MetaData/transcript_ann_hg38.txt')
    # tss_ann = tss_ann[gene_biotype %in% c('protein_coding', 'lincRNA', 'miRNA')]
    # tss_ann[, 'Tss' := ifelse(strand == '+', start, end)]
    
  }else{
    
    #--new version, use longest transcripts
    annotation = readRDS('/PATH/01.ref/GRanges_hg38.rds')
    tss_ann <- Signac:::CollapseToLongestTranscript(annotation)
    tss_ann <- as.data.frame(resize(tss_ann, width = 1, fix = 'start'))
    names(tss_ann)[1] = 'chr'
    tss_ann$Tss = tss_ann$start
    gene_ann = tss_ann = data.table(tss_ann)
    
    gene_ann_sele = gene_ann[gene_name %in% predicted_links$gene_name, ]
    
  }
  
  setkey(gene_ann_sele, gene_name)
  predicted_links[, 'chr' := gene_ann_sele[J(predicted_links$gene_name)]$chr]
  predicted_links[, 'Tss' := gene_ann_sele[J(predicted_links$gene_name)]$Tss]
  predicted_links[, 'promoter_pos' := paste(chr, Tss - 1000, Tss + 1000, sep = '-')]
  
  # any peak within promoter region
  peak.ann = annPeak2Gene(unique(predicted_links$peak_name), tss_ann, 1000)
  setkey(peak.ann, peak_name)
  peaks.nprom = peak.ann[nchar(gene_name) == 0]$peak_name
  peaks.prom = peak.ann[nchar(gene_name) > 0]$peak_name
  
  # predicted_links.ep = predicted_links[peak_name %in% peak.ann$peaks.nprom]
  predicted_links.ep = predicted_links[peak_name %in% peaks.nprom]
  
  ## assign nearest peak to promoter (within 1kb)
  gene_list <- subset(predicted_links.ep, select = c(gene_name, chr, Tss)) %>%
    .[!duplicated(.)]
  gene_list[, 'start' := Tss - 1000]
  gene_list[, 'end' := Tss + 1000]
  gene_list2peak = geneOverlapPeak(gene_list, peak_names = peaks.used,
                                   mid_dist = 1000)
  setkey(gene_list2peak, gene_name)
  predicted_links.ep[, 'promoter_peak' := gene_list2peak[J(predicted_links.ep$gene_name)]$peak_name]
  
  predicted_links.ep = predicted_links.ep[promoter_peak != 'Not_Found'] # filter genes whose promters are not open
  
  predicted_links.ep = subset(predicted_links.ep, select = c(gene_name, promoter_pos, promoter_peak, peak_name, Tss, P_value, p_val_bonf, fdr, Estimate))
    
  predicted_links.ep[, 'start' := as.integer(unlist(strsplit(peak_name, '-'))[2]), by = peak_name]
  predicted_links.ep[, 'end' := as.integer(unlist(strsplit(peak_name, '-'))[3]), by = peak_name]
  
  predicted_links.ep[, 'ep_dist' := abs(Tss - start/2 - end/2)]
  predicted_links.ep = subset(predicted_links.ep, select = c(gene_name, promoter_pos, promoter_peak, peak_name, Tss, ep_dist, P_value, p_val_bonf, fdr, Estimate))
  #names(predicted_links.ep)[4] = 'enhancer_peak'
  fwrite(predicted_links.ep, file = '/PATH/new_regrRes4_gene_peak_ep_links_metacell_IDX_MyeloidPrim_BB_BM_Mye.tsv', sep = '\t')
}



# =============================================================================
# 3. predict tf regulon per cell state
# =============================================================================
library(data.table)
library(Seurat)
library(magrittr)
library(ggplot2)
library(matrixStats)
library(viridis)
library(openxlsx)
library(stringr)
library(enrichR)
`%notin%` = Negate(`%in%`)
library(motifmatchr)
library(GenomicRanges)


#### 3.1.1 motif scan --- enrich Motif at open peak regions ####
library(chromVAR)
library(chromVARmotifs)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(GenomicRanges)


BALL.Obj <- readRDS("/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts.rds")

DefaultAssay(BALL.Obj) <- 'ATAC'

peaks = rownames(BALL.Obj)
peaks = sapply(peaks, function(x) unlist(strsplit(x, ','))[1])
names(peaks) = NULL

pks = tidyr::separate(data = data.table('peak_name' = peaks),
                      col = 'peak_name', into = c('chr', 'start', 'end'),
                      remove = F)
pks$start = as.integer(pks$start)
pks$end = as.integer(pks$end)
pks <- as.data.table(pks)
setkey(pks, chr, start)
# Make a set of peaks
peaks <- GenomicRanges::GRanges(seqnames = pks$chr,
                 ranges = IRanges::IRanges(start = pks$start,
                                  end = pks$end))

motif_ix <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 1e-3)
                        
motif_ix_str <- matchMotifs(human_pwms_v2, peaks, 
                        genome = BSgenome.Hsapiens.UCSC.hg38,
                        p.cutoff = 5e-5)

if(F){
        motif_name = motif_ix@colData$name
        pk_inf = motif_ix@rowRanges
        motif_pk_match <- motif_ix@assays@data$motifMatches
        rownames(motif_pk_match) = paste0(pk_inf@seqnames, '-', pk_inf@ranges)
        names(motif_name) = NULL
        colnames(motif_pk_match) = motif_name      
}

saveRDS(motif_ix, file = '/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_motif_match_pvE-03_CISBP.rds')
saveRDS(motif_ix_str, file = '/PATH/IDX_MyeloidPrim_BB_BM_Mye_Blasts_motif_match_pv5E-05_CISBP.rds')