
# custome wrapper around Find_DEGs

FindAll_DEGs = function(seurat_object,group_var,idents_name = "seurat_clusters",...){
  
  all_marker_lists = list()
  all_clusters = unique(seurat_object@meta.data[,idents_name])
  
  for(i in 1:length(all_clusters)){
    current_cluster = all_clusters[i]
    message(paste0("Running ",i," out of ",length(all_clusters)," clusters..."))
    currentmarkers =  Find_DEGs(seurat_object=seurat_object,subset_ident= current_cluster,idents_name = idents_name,group_var=group_var,...)
    all_marker_lists[[current_cluster]]=currentmarkers
  }
  res=do.call(rbind,all_marker_lists)
  return(res)
}

# run findmarkers with between two groups

Find_DEGs = function(seurat_object,subset_ident,idents_name = "seurat_clusters",group_var,group1 = NULL,group2 =NULL,pval_filter = 5e-3,def_assay="RNA",...){
  
  # req
  require(Seurat)
  
  # set defaults
  DefaultAssay(seurat_object) <- def_assay
  Idents(seurat_object) <- idents_name
  result=data.frame()
  
  # if group1 not specified take first from group_var
  if(is.null(group1)){
    group1 = unique(seurat_object@meta.data[,group_var])[1]
  }
  
  # run function with error handling
  result <- tryCatch({
    conditionGenes = FindMarkers(seurat_object, ident.1 = group1, ident.2 = group2, group.by = group_var, subset.ident = subset_ident,...)
    # name comp
    first_string = paste0(group1,collapse = ".")
    if(is.null(group2)){
      unique_group_values =  unique(seurat_object@meta.data[,group_var])
      if(length(unique_group_values)==2){
        second_string = unique_group_values[unique_group_values!=group1]
      }else{
        second_string = "all"
      }
    }else{
      second_string = paste0(group2,collapse = ".")
    }
    conditionGenes$Identity = subset_ident
    conditionGenes$comparison = paste0(first_string,"_vs_",second_string)
    
    # reformat
    conditionGenes = conditionGenes %>% 
      dplyr::mutate(gene = rownames(conditionGenes)) %>%
      dplyr::select(gene,Identity,comparison,adj_pval = p_val_adj,p_val,contains("pct"),contains("fold"),contains("FC")) %>%
      dplyr::filter(adj_pval<pval_filter)
    
    return(conditionGenes)
  },
  error=function(cond) {
    message(paste("Find_DEGs caused an error, will skip the current comparison:"))
    message(cond)
    print("")
    return(data.frame())
  },
  warning=function(cond) {
    message(paste("Find_DEGs caused a warning:"))
    message(cond)
    print("")
    return(data.frame())
  })
  
  return(result)
}