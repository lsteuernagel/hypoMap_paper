##########
### get_expression_stats
##########

#' Get gene pct and average expression for cells
#' TODO: does nearly the same as the above function! Refactor?
#' @param object object
#' @param cells.1 which cells
#' @param features genes
#' @param thresh.min min expressio
#' @return dataframe

### function to get_expression_stats
get_expression_stats <- function(object, cells.1,features = NULL,  thresh.min = 0) {
  features <- features %||% rownames(x = object)
  # Calculate percent expressed
  #thresh.min <- 0
  pct.1 <- round(
    x = rowSums(x = object@assays$RNA@data[features, cells.1] > thresh.min) /length(x = cells.1),
    digits = 3
  )
  mean.1 <- rowMeans(object@assays$RNA@data[features, cells.1]) #mean.fxn(object[features, cells.1, drop = FALSE])
  expression_stats <- data.frame(features,mean.1, pct.1)
  colnames(expression_stats) <- c("gene", "mean.1", "pct.1")
  return(expression_stats)
}

##########
### add_reduction_seurat
##########

# TODO: also include clustering in this function

#' Adds the integration embedding and calculates umap. 
#' @param seurat_object
#' @param integration_files
#' @param global_seed seed
#' @return .

add_reduction_seurat = function(seurat_object,integration_name,new_name=NULL,integration_path,max_dim=50,global_seed=123,calc_umap=TRUE,k_param_umap=30,overwrite =FALSE,overwrite2=FALSE){
  
  require(Seurat)
  require(dplyr)
  message("Running add_reduction_seurat with calc_umap = ",calc_umap)
  # get current result
  current_file = list.files(integration_path,recursive = TRUE,pattern = integration_name)
  if(is.null(new_name)){new_name=integration_name}
  
  if(length(current_file)!=1){
    message("Did not find an unambigous file corresponding to the provided name. Files: ")
    message(paste0(current_file,collapse = " | "))
    stop("Stop")
  }else{
    message("Found file for ",integration_name)
  }
  if(! new_name %in% names(seurat_object@reductions) | overwrite){
    
    current_embedding = read_embedding(paste0(integration_path,current_file),seurat_object)
    # make dim red
    dimred <- Seurat::CreateDimReducObject(
      embeddings = as.matrix(current_embedding),
      stdev = as.numeric(apply(current_embedding, 2, stats::sd)),
      assay = "RNA",
      key = new_name
    )
    # add
    seurat_object@reductions[[new_name]] = dimred
    
  }else{
    message("Reduction with this name already found. Use overwrite = TRUE to overwrite")
  }
  if(calc_umap){
    message("Calculating umap for ",new_name)
    if((! paste0("umap_",new_name) %in% names(seurat_object@reductions)) | overwrite2){
      col_dimred = ncol(seurat_object@reductions[[new_name]])
      #print(seurat_object@reductions[[new_name]]@key)
      seurat_object = RunUMAP(seurat_object,reduction = new_name,seed.use=global_seed,dims=1:min(col_dimred,max_dim),n.neighbors=k_param_umap,
                              reduction.name=paste0("umap_",new_name),reduction.key = paste0("umap_",new_name),verbose=FALSE)
    }else{
      message("Umap with this name already found. Use overwrite2 = TRUE to overwrite")
    }
  }
  
  return(seurat_object)
  
}

##########
### read_embedding
##########

#' Load an emebedding with cells x lowDims from flatfile, ensuring consistency with a Seurat object (or metadata only for faster usage)
#' @param filename_withpath filepath
#' @param seurat_object seuratobject associated with current embedding. If specified metadata does not have to be set explicitly.
#' @param seurat_object_metadata metadata only of seuratobject associated with current embedding
#' @return 

read_embedding = function(filename_withpath,seurat_object=NULL,seurat_object_metadata=NULL){
  
  #get metadata
  if(!is.null(seurat_object_metadata)){
    metadata = seurat_object_metadata
    rownames(metadata) = metadata[,"Cell_ID"] # manually set rownames
  }else{
    if(!is.null(seurat_object)){
      metadata = seurat_object@meta.data
    }else{
      stop("Please provide either a dataframe with metadata or a seurat object with metadata that can be exctracted!")
    }
  }
  # load
  current_embedding = data.table::fread(filename_withpath,data.table = F)
  # use first col as rownames 
  if(is.character(current_embedding[,1])){
    # message("Using first column of loaded file as rownames for reduction")
    rnames = current_embedding[,1]
    current_embedding = current_embedding[,2:ncol(current_embedding)]
    rownames(current_embedding)=rnames
    # reorder to align with rest of object
    if(any(is.na(match(rownames(metadata),rownames(current_embedding))))){
      stop("Cell names from loaded reduction and new object are not matching exactly. Stopping import.")
    }
    current_embedding = current_embedding[match(rownames(metadata),rownames(current_embedding)),]
  }else{
    warning("First column of loaded file is not of type character, using rownames of metadata as rownames of added reduction. This can induce bugs if the order changed due to split/merge of the Seurat object!")
    rownames(current_embedding) = rownames(metadata)
  }
  return(current_embedding)
  
}

##########
### Custom DEG functions
##########

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
