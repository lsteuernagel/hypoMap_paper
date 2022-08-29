##########
### load_required_files
##########

#' This function is a general helper that loads required objects via the large data paths (has to be set to local system directory)
#' Currently loads the neuron and full hypoMap and the neuron version of the nuc-seq data.
#' @param large_data_path large_data_path
#' @param overwrite_existing which cells
#' @param filenames genes

# function 
load_required_files <- function(large_data_path, overwrite_existing = FALSE, filenames = c("hypoMap_v2_seurat" = "hypoMap_v2.rds")) {
  # go through all files
  for(i in 1:length(filenames)){
    objectname = names(filenames)[i]
    filename =  filenames[i]
    if(!exists(objectname,envir = .GlobalEnv) | overwrite_existing){
      message("Loading ", objectname," from ",large_data_path,filename)
      if(file.exists(paste0(large_data_path,filename))){
        assign(objectname, readRDS(paste0(large_data_path,filename)), envir = .GlobalEnv)
      }else{
        message("Warning: Cannot find file ",filename) 
      }
    }else{
      message(objectname, " already exists in .GlobalEnv. Skipping.")
    }
  }
}

##########
### load_plot_params
##########

#' Load color palettes and background colors into global environment.

# function 
load_plot_params <- function() {
  
  # general backgorund color (mostly for UMAPs etc)
  rasterize_px = 2048
  assign(x = "rasterize_px",rasterize_px, envir = .GlobalEnv)
  seurat_pt_size = 2
  assign(x = "seurat_pt_size",seurat_pt_size, envir = .GlobalEnv)
  rasterize_pixels = 2048
  assign(x = "rasterize_pixels",rasterize_pixels, envir = .GlobalEnv)
  rasterize_point_size = 2
  assign(x = "rasterize_point_size",rasterize_point_size, envir = .GlobalEnv)
  text_size = 20
  assign(x = "text_size",text_size, envir = .GlobalEnv)
}

##########
### load_colors
##########

#' Load color palettes and background colors into global environment.

# function 
load_colors <- function() {
  
  # general backgorund color (mostly for UMAPs etc)
  bg_col = "grey90"
  assign(x = "bg_col",bg_col, envir = .GlobalEnv)
  cols_for_feature_plot = c(bg_col,"#0b3ebd") # "#0b3ebd"
  assign(x = "cols_for_feature_plot",cols_for_feature_plot, envir = .GlobalEnv)
  
  fasting_color = "#28a7c9"
  assign(x = "fasting_color",fasting_color, envir = .GlobalEnv)
  adlib_color = "#772ac9"
  assign(x = "adlib_color",adlib_color, envir = .GlobalEnv)
  
  reference_sc_color = "#cc2118"
  assign(x = "reference_sc_color",reference_sc_color, envir = .GlobalEnv)
  query_sn_color = "#302ac9"
  assign(x = "query_sn_color",query_sn_color, envir = .GlobalEnv)
  # color-blind friendly color palettes:
  # see: https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
  # deleted greens!
  long_palette <- c(
    "dodgerblue2", "#E31A1C", # red
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1", "skyblue2", "#FB9A99", # lt pink
    "#CAB2D6", # lt purple
    "#FDBF6F", # lt orange
    "khaki2","maroon", "orchid1", "deeppink1", "blue1", "steelblue4","darkturquoise", "yellow4", "yellow3","darkorange4", "brown"
  )
  assign(x = "long_palette",long_palette, envir = .GlobalEnv)
  long_palette_strong <- c(
    "dodgerblue2", "#E31A1C", # red
    "#6A3D9A", # purple
    "#FF7F00", # orange
    "black", "gold1" ,"maroon", "deeppink1", "blue1", "steelblue4","darkturquoise", "yellow4", "yellow3","darkorange4", "brown"
  )
  assign(x = "long_palette_strong",long_palette_strong, envir = .GlobalEnv)
  
  alphabet_palette = as.character(palette.colors(palette = "Alphabet")[!names(palette.colors(palette = "Alphabet")) %in% c("forest","green","iron","jade","quagmire","ultraviolet","violet","wine","xanthin","yellow","zinnia") ])
  assign(x = "alphabet_palette",alphabet_palette, envir = .GlobalEnv)
  
  # see: https://stackoverflow.com/questions/57153428/r-plot-color-combinations-that-are-colorblind-accessible
  short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
  short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
  assign(x = "short_palette",short_palette, envir = .GlobalEnv)
}

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
### extend_sibling_markers
##########

## to add missing sibling markers

extend_sibling_markers = function(sibling_markers,edgelist,annotation){
  
  updated_sibling_markers = sibling_markers
  # for eah cluster find parent in edgelist
  for(i in 1:nrow(edgelist)){
    current_cluster = edgelist$to[i]
    current_cluster_name = annotation$cluster_name[annotation$cluster_id == current_cluster]
    # if parent has count==1 and sibling_markers does not contain id:
    if(!edgelist$from[i] %in% edgelist$to){next}
    if(edgelist$count[edgelist$to == edgelist$from[i]] == 1 & nrow(updated_sibling_markers[updated_sibling_markers$cluster_id==current_cluster,])==0){
      # extract parent sibling markers 
      parent_id = edgelist$from[i]
      new_sibling_markers = updated_sibling_markers[updated_sibling_markers$cluster_id == parent_id,]
      
      # oervwrite with current node id
      if(nrow(new_sibling_markers) > 0){
        new_sibling_markers$cluster_id = current_cluster
        new_sibling_markers$cluster_name = current_cluster_name
      }
      # add rows to all marker df
      updated_sibling_markers = dplyr::bind_rows(updated_sibling_markers,new_sibling_markers)
      
    }
    
    # order by edgelist
    updated_sibling_markers$cluster_id = factor(updated_sibling_markers$cluster_id,levels = edgelist$to)
    updated_sibling_markers = updated_sibling_markers %>% dplyr::arrange(cluster_id)
    updated_sibling_markers$cluster_id = as.character( updated_sibling_markers$cluster_id )
  }
  return(updated_sibling_markers)
}


##########
### RBO functions
##########

# I am using code from here:
# https://rdrr.io/bioc/gespeR/src/R/gespeR-concordance.R 
# https://www.bioconductor.org/packages/release/bioc/html/gespeR.html 
# The gesper package provides a nice implementation of rbo but is quite dependency heavy, so I extracted the function from there.

#' Rank biased overlap (Webber et al., 2010)
#' 
#' Evaluates the rank biased overlap (rbo) of two ranked lists based on formula based on (32) from 
#' "A Similarity Measure for Indefinite Rankings" (Webber et al.). Two ranked lists with high rbo are
#' very similar, wheras low rbo indicates dissimilar lists. rbo ranges between 0 and 1. In this method
#' the extrapolated version of rbo is implemented.
#' 
#' @author Fabian Schmich
#' @export
#' 
#' @param s List 1
#' @param t List 2
#' @param p Weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k Evaluation depth for extrapolation
#' @param side Evaluate similarity between the top or the bottom of the ranked lists
#' @param mid Set the mid point to for example only consider positive or negative scores
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @return rank biased overlap (rbo)
#' 
#' @seealso \code{\link{concordance}}
#' 
#' @examples
#' a <- rnorm(26)
#' b <- rnorm(26)
#' names(a) <- names(b) <- LETTERS
#' rbo(a, b, p = 0.95)
rbo2 <- function(s, t, p, k=floor(max(length(s), length(t))/2), side=c("top", "bottom"), mid=NULL, uneven.lengths = TRUE) {
  side <- match.arg(side)
  if (!is.numeric(s) | !is.numeric(t))
    stop("Input vectors are not numeric.")
  if (is.null(names(s)) | is.null(names(t)))
    stop("Input vectors are not named.")
  ids <- switch(side,
                "top"=list(s=.select.ids(s, "top", mid), t=.select.ids(t, "top", mid)),
                "bottom"=list(s=.select.ids(s, "bottom", mid), t=.select.ids(t, "bottom", mid))
  )
  min(1, .rbo.ext(ids$s, ids$t, p, k, uneven.lengths = uneven.lengths))
}

#' Select top or bottom names of ranked vector
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param x The ranked list
#' @param side The side to be evaluated ("top" or "bottom" of ranked list)
#' @param mid The mid point to split a list, e.g. to split between positive and negative values choose mid=0
#' @return A vector of selected identifiers
.select.ids <- function(x, side=c("top", "bottom"), mid=NULL) {
  side <- match.arg(side)
  if (side == "top")  {
    x <- sort(x, decreasing=TRUE)
    if (is.null(mid))
      return(names(x))
    else 
      return(names(x)[which(x > mid)])
  } else if (side == "bottom") {
    x <- sort(x, decreasing=FALSE)
    if (is.null(mid)) 
      return(names(x))
    else 
      return(names(x)[which(x < mid)])
  }
}

#' Rank biased overlap formula based on (32) from "A Similarity Measure for Indefinite Rankings" (Webber et al.)
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param x List 1
#' @param y List 2
#' @param p The weighting parameter in [0, 1]. High p implies strong emphasis on top ranked elements
#' @param k The evaluation depth
#' @param uneven.lengths Indicator if lists have uneven lengths
#' @return The rank biased overlap between x and y
.rbo.ext <- function(x, y, p, k, uneven.lengths = TRUE) {
  if (length(x) <= length(y)) {
    S <- x
    L <- y
  } else {
    S <- y
    L <- x
  }
  l <- min(k, length(L))
  s <- min(k, length(S))
  
  if (uneven.lengths) {
    Xd <- sapply(1:l, function(i) length(intersect(S[1:i], L[1:i])))
    ((1-p) / p) *
      ((sum(Xd[seq(1, l)] / seq(1, l) * p^seq(1, l))) +
         (sum(Xd[s] * (seq(s+1, l) - s) / (s * seq(s+1, l)) * p^seq(s+1, l)))) +
      ((Xd[l] - Xd[s]) / l + (Xd[s] / s)) * p^l  
  } else {
    #stopifnot(l == s)
    k <- min(s, k)
    Xd <- sapply(1:k, function(i) length(intersect(x[1:i], y[1:i])))
    Xk <- Xd[k]
    (Xk / k) * p^k + (((1-p)/p) * sum((Xd / seq(1,k)) * p^seq(1,k)))
  }
}

#' Jaccard index between two sets
#' 
#' @author Fabian Schmich
#' @noRd
#' 
#' @param x Set 1
#' @param y Set 2
#' @return The Jaccard index
.jaccard <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}
