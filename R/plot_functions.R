
##########
### rasterize_ggplot
##########

#' Rasterize ggplot
#'
#' Makes use of scattermore::geom_scattermore to substitute GeomPoint layer in ggplot objects with an equivalent rasterized version.
#' Allows to manually modify the resolution of rasterized plots from Seurat plots!
#'
#' @param plot ggplot2 or patchwork object. For example output of Seurat's DimPlot or FeaturePlot
#' @param pixel_raster integer: number of pixels passed to pixels argument (x and optionall y) from scattermore::geom_scattermore. Defaults to 1024
#' @param pixel_raster_y pixel_raster_y: to use a different y pixels than x. Defaults to NULL which will use the value from pixel_raster in x and y
#' @param interpolate whether to linearly interpolate the image. Defaults to FALSE
#' @param pointsize Radius of rasterized point from scattermore. defaults to 1
#'
#' @return ggplot object with GeomPoint layers changed to GeomScattermore with rastrized version
#'
#' @import ggplot2 scattermore

rasterize_ggplot = function(plot,pixel_raster = 1024,pixel_raster_y = NULL,interpolate=FALSE,pointsize = 1){
  # if NULL will use pixel_raster else use different for y
  if(is.null(pixel_raster_y)){pixel_raster_y = pixel_raster}
  # get the geom point mapping
  layer_idx_all = which(sapply(plot$layers, function(x) class(x$geom)[1])=="GeomPoint")
  # check if any GeomPoint mapping exists
  if(length(layer_idx_all)==0){
    stop("Error: Cannot find points. 'rasterize_ggplot' expects to find at least one valid 'GeomPoint' layer to build GeomScattermore from.")
  }
  rasterized_plot = plot
  # in case there are multiple GeomPoint mappings: do for each:
  for(layer_idx in layer_idx_all){
    geom_point_layer = plot$layers[[layer_idx]]
    # make a plot with a rasterized geom
    rasterized_plot = rasterized_plot + scattermore::geom_scattermore(
      mapping = aes_string(
        x = as_label(geom_point_layer$mapping$x),
        y = as_label(geom_point_layer$mapping$y),
        color = paste0("`", as_label(geom_point_layer$mapping$colour), "`"),
        shape = as_label(geom_point_layer$mapping$shape),
        alpha = as_label(geom_point_layer$mapping$alpha)
      ),
      interpolate = interpolate,
      pointsize = pointsize,
      pixels = c(pixel_raster,pixel_raster_y)
    )
    # move into geom_point_mapping level in object (to preserve order of layers)
    rasterized_plot$layers[[layer_idx]] = rasterized_plot$layers[[length(rasterized_plot$layers)]]
    rasterized_plot$layers[[length(rasterized_plot$layers)]] = NULL
    #check if the new geom has non-default data that should be added 
    if(length(geom_point_layer$data)>0){
      rasterized_plot$layers[[layer_idx]]$data = geom_point_layer$data
    }
  }
  return(rasterized_plot)
}



##########
### plot_cluster_tree
##########

#' Plot metadata entries in cluster tree using ggtree
#'
#' TODO: add description
#' TODO: add reference to ggtree vignette
#' https://github.com/YuLab-SMU/ggtree/issues/400#issuecomment-845781670
#'
#' @param edgelist TODO
#' @param heatmap_matrix a matrix with rownames corresponding to tip node names
#' @param heatmap_matrix a matrix with rownames corresponding to tip node names. if not NULL will set a second heatmap ring (for more pleas adjust manually)
#' @param leaf_level which level to use as leaves ?
#' @param anno_df provide a dataframe which contains mappings of all cluster ids and names in the cluster tree. Defaults to NULL which means automatic construction.
#' @param metadata metadata with cluster ids and names, if anno_df is NULL this has to be provided to construct anno_df
#' @param level_pattern regex for cluster level, if anno_df is NULL this has to be provided to construct anno_df. defaults to 'K[0-9]+'
#' @param cluster_id_pattern regex for cluster id column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_pruned'
#' @param cluster_name_pattern regex for cluster name column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_named'
#' @param label_size label size on tree. defaults to 2
#' @param show_genes TRUE
#' @param legend_title "Legend"
#' @param matrix_offset 0.2
#' @param matrix_width 0.2
#' @param heatmap_colnames TRUE
#' @param heatmap_colors colors to be passed to scale_fill_gradientn
#' @param returnData false
#'
#' @return ggtree object or plot
#'
#'

# load for tests:
# edgelist =  neuron_map_seurat@misc$pruned_edgelist
# leaf_level = 6
# metadata = neuron_map_seurat@meta.data
# level_pattern = "K[0-9]+"
# cluster_id_pattern = "_pruned"
# cluster_name_pattern = "_named"
# anno_df=NULL
# label_size = 2
# heatmap_data = metadata %>% dplyr::select(Cell_ID,K169_named) %>% dplyr::group_by(K169_named) %>%  #dplyr::filter(predicted_Campbell!="NA")
#   dplyr::add_count(name = "presence") %>% dplyr::distinct(K169_named,.keep_all=TRUE) %>%dplyr::ungroup() %>% dplyr::mutate(presence = presence / sum(presence)*100) %>% dplyr::ungroup()# %>% #%>%  dplyr::left_join(tree_data_tibble[,c("label","node")],by=c("K169"="label"))
# # tidyr::spread(key = 1,value=presence)
# heatmap_data = heatmap_data %>% dplyr::full_join(anno_df[,c("cluster_id","cluster_name")],by=c("K169_named"="cluster_name"))  %>% dplyr::filter(grepl("K169",cluster_id))
# heatmap_matrix = as.matrix(heatmap_data[,"presence"])
# rownames(heatmap_matrix) = heatmap_data$cluster_id
# heatmap_matrix[is.na(heatmap_matrix)] = 0

plot_cluster_tree = function(edgelist,heatmap_matrix=NULL,heatmap_matrix2 = NULL,leaf_level=NULL,anno_df=NULL,metadata=NULL,level_pattern = "K[0-9]+",cluster_id_pattern = "_pruned",
                             cluster_name_pattern = "_named",label_size = 2, show_genes = TRUE,heatmap_colors=c("white","darkred"),heatmap_colnames =TRUE, legend_title_1 = "Legend1",legend_title_2="Legend2",
                             matrix_offset = 0.2,matrix_width =0.2,matrix_width_2=0.2,colnames_angle=0,legend_text_size = 4,hjust_colnames=0.5, returnData =FALSE,manual_off_second = 5,heatmap_text_size=4,
                             annotate_reverse =TRUE){
  
  # I am loading packages here instead of using them as part of a pachage and playing around with namespaces because for some reason the fucks up the ggtree package
  # for example geom_nodelab behaves different when called via ggtree::geom_nodelab, indicating taht is uses some other version?
  # also ggtree has some bugs, like not working with dplyr > 1.06
  library(dplyr)
  library(igraph)
  library(ggplot2)
  library(scales)
  library(treeio)
  library(tidytree)
  library(ggtree)
  library(ggnewscale)
  
  # check that req columns exist
  if(length(setdiff(c("from","to","level"),colnames(edgelist)))>0){
    
    warning("Error: Wrong edgelist format. Requires columns: from, to, level")
    return(NULL)
  }
  # check that leaf_level is there
  if(!leaf_level %in% edgelist$level){
    warning("Error: leaf_level '",leaf_level,"' cannot be found in level column of edgelist")
    return(NULL)
  }
  
  # construct a dataframe with the required annotations
  if(is.null(anno_df)){
    if(is.null(metadata)){
      warning("Error: Please provide metadata with corresponding cluster ids and names that match provided edgelist")
      return(NULL)
    }
    if(!any(grepl(cluster_id_pattern,colnames(metadata)))){
      warning("Error: Cannot find columns with cluster_id_pattern '",cluster_id_pattern,"' in provided metadata")
      return(NULL)
    }
    if(!any(grepl(cluster_name_pattern,colnames(metadata)))){
      warning("Error: Cannot find columns with cluster_name_pattern '",cluster_name_pattern,"' in provided metadata")
      return(NULL)
    }
    # anno_df: cluster id, clean_names, clusterlevel, ncells, first_cluster_name
    pruned_ids =metadata[,c("Cell_ID",colnames(metadata)[grepl(cluster_id_pattern,colnames(metadata))])] %>%
      tidyr::gather(-Cell_ID,key="colname",value="id") %>% dplyr::mutate(clusterlevel = stringr::str_extract(colname,level_pattern)) %>% dplyr::distinct(id,clusterlevel,.keep_all=TRUE)
    named_ids =metadata[,c("Cell_ID",colnames(metadata)[grepl(cluster_name_pattern,colnames(metadata))])] %>%
      tidyr::gather(-Cell_ID,key="colname",value="id") %>% dplyr::mutate(clusterlevel = stringr::str_extract(colname,level_pattern))  %>% dplyr::distinct(id,clusterlevel,.keep_all=TRUE)
    both_map = dplyr::left_join(pruned_ids,named_ids,by=c("Cell_ID"="Cell_ID","clusterlevel"="clusterlevel")) %>% dplyr::select(cluster_id = id.x,cluster_name = id.y)
    
    anno_df = neuron_map_seurat@misc$pruned_edgelist %>% dplyr::select(cluster_id = to, clusterlevel = clusterlevel,ncells ) %>% dplyr::left_join(both_map,by="cluster_id")
    if(annotate_reverse){
      anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][1]})
    }else{
      anno_df$first_cluster_name = sapply(anno_df$cluster_name,function(x){strsplit(x,"\\.")[[1]][length(strsplit(x,"\\.")[[1]])]})
    }
  }else{
    # check that provided anno_df is valid:
    if(length(setdiff(c("cluster_id","clusterlevel","cluster_name","first_cluster_name"),colnames(anno_df)))>0){
      stop("Wrong anno_df format. Required columns: cluster_id, clusterlevel, cluster_name, first_cluster_name")
    }
  }
  
  # if a heatmap matrix is provided, this function tries to infer the leaflevel based on the matrix
  if(!is.null(heatmap_matrix) & is.null(leaf_level)){
    # TODO
  }
  
  # reduce edgelist to certain level and from and to cols
  edgelist$level = as.numeric(edgelist$level)
  edgelist = edgelist[edgelist$level<=as.numeric(leaf_level),1:2]
  edgelist = edgelist[edgelist$to %in% anno_df$cluster_id,] # remove edges/nodes that are not part of anno_df
  
  ## convert to treedata
  # only take
  tree_data_igraph = base::suppressWarnings(igraph::graph_from_edgelist(as.matrix(edgelist)))
  tree_data_phylo = base::suppressWarnings(treeio::as.phylo(tree_data_igraph))
  tree_data_tibble <- dplyr::as_tibble(tree_data_phylo)
  
  # add labels from annotation_df
  tree_data_tibble = dplyr::left_join(tree_data_tibble,anno_df,by=c("label"="cluster_id"))
  
  # update additional columns
  tree_data_tibble$first_cluster_name[is.na(tree_data_tibble$first_cluster_name)]=""
  tree_data_tibble$nodesize = 1 # default node size
  tree_data_tibble$n_children = sapply(tree_data_tibble$label,function(x,el){length(el$to[el$from==x])},el=edgelist) # count children number
  tree_data_tibble$n_siblings = sapply(tree_data_tibble$label,function(x,el){ # count siblings
    parent = el$from[el$to==x]
    return(length(el$to[el$from==parent])-1)
  },el=edgelist)
  tree_data_tibble$tip.label =NA
  tree_data_tibble$tip.label[tree_data_tibble$n_children==0] = tree_data_tibble$node[tree_data_tibble$n_children==0] # add tip labels if leaf
  tree_data_tibble$first_cluster_name[ tree_data_tibble$n_children<2 & is.na(tree_data_tibble$tip.label)] = "" # if only one child node: changeto ""
  
  # convert back to treedata
  tree_data = suppressWarnings(tidytree::as.treedata(tree_data_tibble))
  
  #plot circular tree
  circular_tree =ggtree(tree_data,layout = 'circular', branch.length='none')+
    geom_nodepoint(aes(subset = n_children > 1))#+geom_tippoint() + layout_circular()
  # add genes on edges if necessary
  if(show_genes){
    circular_tree = circular_tree +
      geom_nodelab(aes(x=branch, label=first_cluster_name), size=label_size,vjust=-.5, color="darkred")+
      geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name), size=label_size,vjust=-.5,color="darkred")
  }
  
  
  if(!is.null(heatmap_matrix)){
    ## add heatmap
    
    #circular_tree
    
    # optionally add custom cont scale
    if(is.numeric(heatmap_matrix[,1])){
      scale_limits = c(min(heatmap_matrix),max(heatmap_matrix))
      circular_tree_heat <- gheatmap(circular_tree,data=heatmap_matrix,offset = matrix_offset , width = matrix_width,colnames = heatmap_colnames,
                                     colnames_angle=colnames_angle,font.size=heatmap_text_size,hjust = hjust_colnames)+
        scale_fill_gradientn(colours = heatmap_colors,limits=scale_limits,oob=squish) +
        theme(legend.text=element_text(size=legend_text_size)) +
        guides(fill=ggplot2::guide_colourbar(title=legend_title_1)) # guide_colourbar for continous
    }else{
      circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix, offset=matrix_offset, width=matrix_width,colnames = heatmap_colnames,
                                     colnames_angle=colnames_angle,font.size=heatmap_text_size,hjust = hjust_colnames)+
        scale_fill_discrete() + 
        theme(legend.text=element_text(size=legend_text_size))+
        guides(fill=ggplot2::guide_legend(title=legend_title_1))
    }
    ## add a second heatmap
    if(!is.null(heatmap_matrix2)){
      # make new fill scale
      circular_tree_heat <- circular_tree_heat + ggnewscale::new_scale_fill()
      # optionally add custom cont scale
      if(is.numeric(heatmap_matrix2[,1])){
        scale_limits = c(min(heatmap_matrix2),max(heatmap_matrix2))
        circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix2, offset=manual_off_second, width=matrix_width_2,
                                       colnames = heatmap_colnames,colnames_angle=colnames_angle, colnames_offset_y = matrix_offset*2,
                                       font.size=heatmap_text_size,hjust = hjust_colnames)+
          scale_fill_gradientn(colours = heatmap_colors,limits=scale_limits,oob=squish) +
          guides(fill=ggplot2::guide_colourbar(title=legend_title_2)) + # guide_colourbar for continous  +
          theme(legend.text=element_text(size=legend_text_size))
      }else{
        circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix2, offset=manual_off_second,
                                       width=matrix_width_2,colnames = heatmap_colnames,colnames_angle=colnames_angle,
                                       colnames_offset_y = matrix_offset*2,font.size=heatmap_text_size,hjust = hjust_colnames)+
          scale_fill_discrete() +
          guides(fill=ggplot2::guide_legend(title=legend_title_2)) +
          theme(legend.text=element_text(size=legend_text_size))
      }
    }
  }else{
    circular_tree_heat = circular_tree
  }
  
  if(returnData){
    return(circular_tree_heat)
  }else{
    circular_tree_heat
  }
  
}

