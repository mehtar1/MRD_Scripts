#c alculate entropy per sample


get.entropy = function(seurat.obj, cluster_columns, join.to.seurat = F){
  
  entropy.list = list()
  
  for (i in (1:length(cluster_columns))){
    
    cluster_column = cluster_columns[i]
    
    print(paste("calculating entropy for", cluster_column))
    
    number.unique =  seurat.obj@meta.data[[cluster_column]] %>% unique() %>% length()
    
    print(paste(number.unique, "unique levels in", cluster_column))
    
    print("calculating entropy")

    p.log.p.table = seurat.obj@meta.data %>% 
      group_by(.data[[cluster_column]], Projected.Ctype.Wmethod, .drop = T) %>% #group by predicted cell type & sample name
      dplyr::summarize(n = length(Projected.Ctype.Wmethod)) %>% #sum
      dplyr::left_join(seurat.obj@meta.data %>% 
                  group_by(.data[[cluster_column]]) %>% 
                  dplyr::summarize(total = length(.data[[cluster_column]])))%>%  #calculate total cells / sample
      dplyr::mutate(p = n/total, plogp = p*log(p, base = 2))
    
    entropy = p.log.p.table %>% 
      group_by(.data[[cluster_column]]) %>% 
      dplyr::summarize(entropy = -sum(plogp) / log(number.unique,2)) #calculate entropy by negative summing p log p across each cluster
    
    colnames(entropy) = c(cluster_column, paste0(cluster_column, "_entropy"))
    entropy.list[[cluster_column]] = entropy
    #print(head(entropy))
    
    #print(entropy)
    
    seurat.obj$cell.name =rownames(seurat.obj@meta.data)
    #print(head(seurat.obj@meta.data))
    temp <- seurat.obj@meta.data %>% dplyr::left_join(entropy) %>% dplyr::select(paste0(cluster_column, "_entropy")) 
    rownames(temp) <- rownames(seurat.obj@meta.data)
    print(temp)
    print(rownames(temp))
    seurat.obj <- AddMetaData(seurat.obj, metadata = temp, col.name = paste0(cluster_column, "_entropy"))
    
  }
  
  entropy.names <- paste0(cluster_columns, "_entropy")
  entropy.individual <- seurat.obj@meta.data %>% as.data.frame() %>% dplyr::select(all_of(entropy.names))
  print(entropy.individual)
  # #print(colnames(seurat.obj@meta.data))
  entropy.table <- rowMeans(entropy.individual, na.rm = T) %>% as.data.frame()
  colnames(entropy.table) <- "average_entropy"
  entropy.table <- entropy.table %>% rownames_to_column(var = "cell_name_entropy")

  
  if(join.to.seurat == T){
    return(seurat.obj)
  }else{

    return(entropy.list)
  }
  
  
}
