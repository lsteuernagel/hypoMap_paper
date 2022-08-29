
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
    # if mapping is not in layer: substitute with global
    if(is.null(geom_point_layer$mapping)){
      geom_point_layer$mapping = plot$mapping
    }
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
#' @param leaf_level which level to use as leaves ?
#' @param anno_df provide a dataframe which contains mappings of all cluster ids and names in the cluster tree. Defaults to NULL which means automatic construction.
#' @param metadata metadata with cluster ids and names, if anno_df is NULL this has to be provided to construct anno_df
#' @param level_pattern regex for cluster level, if anno_df is NULL this has to be provided to construct anno_df. defaults to 'K[0-9]+'
#' @param cluster_id_pattern regex for cluster id column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_pruned'
#' @param cluster_name_pattern regex for cluster name column names in metadata, if anno_df is NULL this has to be provided to construct anno_df. defaults to '_named'
#' @param label_size label size on tree. defaults to 2
#' @param label_size_tip label size of lowest level on tree. label_size
#' @param show_genes TRUE
#' @param edge_color_manual defaults to NULL
#'
#' @return ggtree object or plot
#'

plot_cluster_tree = function(edgelist,leaf_level=NULL,anno_df=NULL,metadata=NULL,pruned_edgelist=NULL,level_pattern = "K[0-9]+",cluster_id_pattern = "_pruned",
                             cluster_name_pattern = "_named",label_size = 2,label_size_tip = label_size, show_genes = TRUE,vjust_label = -0.5,
                             annotate_reverse =TRUE,na_color = "grey90",edge_color="black",node_color="black"){
  
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
    
    # anno_df = neuron_map_seurat@misc$pruned_edgelist %>% dplyr::select(cluster_id = to, clusterlevel = clusterlevel,ncells ) %>% dplyr::left_join(both_map,by="cluster_id")
    anno_df =pruned_edgelist %>% dplyr::select(cluster_id = to, clusterlevel = clusterlevel,ncells ) %>% dplyr::left_join(both_map,by="cluster_id")
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
  # if(!is.null(heatmap_matrix) & is.null(leaf_level)){
  #   # TODO
  # }
  # 
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
  circular_tree =ggtree(tree_data,layout = 'circular', branch.length='none',color=edge_color)+
    geom_nodepoint(aes(subset = n_children > 1),color=node_color)#+geom_tippoint() + layout_circular()
  # add genes on edges if necessary
  if(show_genes){
    circular_tree = circular_tree +
      geom_nodelab(aes(x=branch, label=first_cluster_name), size=label_size,vjust=vjust_label, color="darkred")+
      geom_tiplab(ggplot2::aes(x=branch, label=first_cluster_name), size=label_size_tip,vjust=vjust_label,color="darkred")
  }
  
  return(circular_tree)
  
  # circular_tree_heat <- circular_tree
  # 
  # for(i in 1:length(heatmap_matrix_list)){
  #   
  #   ## add heatmap
  #   heatmap_matrix = heatmap_matrix_list[[i]]
  #   circular_tree_heat <- circular_tree_heat + ggnewscale::new_scale_fill()
  #   #circular_tree
  #   # make new fill scale
  #   circular_tree_heat <- circular_tree_heat + ggnewscale::new_scale_fill()
  #   # optionally add custom cont scale
  #   if(is.numeric(heatmap_matrix2[,1])){
  #     scale_limits = c(min(heatmap_matrix2),max(heatmap_matrix2))
  #     circular_tree_heat <- gheatmap(circular_tree, heatmap_matrix2, offset=manual_off_second, width=matrix_width_2,
  #                                    colnames = heatmap_colnames,colnames_angle=colnames_angle, colnames_offset_y = matrix_offset*2,
  #                                    font.size=heatmap_text_size,hjust = hjust_colnames)+
  #       scale_fill_gradientn(colours = heatmap_colors,limits=scale_limits,oob=squish) +
  #       guides(fill=ggplot2::guide_colourbar(title=legend_title_2)) + # guide_colourbar for continous  +
  #       theme(legend.text=element_text(size=legend_text_size))
  #   }else{
  #     circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix2, offset=manual_off_second,
  #                                    width=matrix_width_2,colnames = heatmap_colnames,colnames_angle=colnames_angle,
  #                                    colnames_offset_y = matrix_offset*2,font.size=heatmap_text_size,hjust = hjust_colnames)+
  #       scale_fill_discrete(na.value = na_color) +
  #       guides(fill=ggplot2::guide_legend(title=legend_title_2)) +
  #       theme(legend.text=element_text(size=legend_text_size))
  #   }
  # }
  # 
  # if(returnData){
  #   return(circular_tree_heat)
  # }else{
  #   circular_tree_heat
  # }
  
}

##########
### add_heatmap
##########

#' Plot metadata entries in cluster tree using ggtree
#'
#' @param circular_tree
#' @param heatmap_matrix a matrix with rownames corresponding to tip node names
#' @param heatmap_colnames olors to be passed to scale_fill_gradientn
#' @param legend_title "Legend"
#' @param matrix_offset 0.2
#' @param matrix_width 0.2
#' @param colnames_angle
#' @param legend_text_size
#' @param hjust_colnames colors to be passed to scale_fill_gradientn
#' @param returnData false
#'
#' @return ggtree object or plot
#'
#'

add_heatmap = function(circular_tree,heatmap_matrix,heatmap_colors=c("white","darkred"),scale_limits = NULL,heatmap_colnames =TRUE, legend_title = "Legend1",
                       matrix_offset = 0.2,matrix_width =0.2,colnames_angle=0,legend_text_size = 4,hjust_colnames=0.5, 
                       heatmap_text_size=4,na_color = "grey90"){
  
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
  library(scales)
  
  circular_tree_heat <- circular_tree
  # make new fill scale
  circular_tree_heat <- circular_tree_heat + ggnewscale::new_scale_fill()
  # optionally add custom cont scale
  if(is.numeric(heatmap_matrix[,1])){
    if(is.null(scale_limits)){
      scale_limits = c(min(heatmap_matrix),max(heatmap_matrix))
    }
    circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix, offset=matrix_offset, width=matrix_width,
                                   colnames = heatmap_colnames,colnames_angle=colnames_angle, colnames_offset_y = matrix_offset,
                                   font.size=heatmap_text_size,hjust = hjust_colnames)+
      scale_fill_gradientn(colours = heatmap_colors,limits=scale_limits,oob=squish,na.value = na_color) +
    #  guides(fill=ggplot2::guide_colourbar(title=legend_title)) + # guide_colourbar for continous  +
      theme(legend.text=element_text(size=legend_text_size))
  }else{
    circular_tree_heat <- gheatmap(circular_tree_heat, heatmap_matrix, offset=matrix_offset,
                                   width=matrix_width,colnames = heatmap_colnames,colnames_angle=colnames_angle,
                                   colnames_offset_y = matrix_offset,font.size=heatmap_text_size,hjust = hjust_colnames)+
      scale_fill_discrete(na.value = na_color) +
    #  guides(fill=ggplot2::guide_legend(title=legend_title)) +
      theme(legend.text=element_text(size=legend_text_size))
  }
  
  return(circular_tree_heat)
  
}


##########
### get_coordinates
##########

# get_coordinates for cluster/label centers 

get_coordinates = function(label_vector,label_column,seurat_object,reduction_name = "umap_scvi"){
  plotdata = bind_cols(seurat_object@reductions[[reduction_name]]@cell.embeddings,label=seurat_object@meta.data[,label_column]) %>%
    dplyr::filter(label %in% label_vector)
  label_centers = plotdata %>% dplyr::group_by(label) %>% dplyr::summarise_at(vars(matches("umap")), median)
  return(label_centers)
}

##########
### barplots_on_scatter
##########

# data_barplot requires columns: cluster and value

barplots_on_umap = function(scatter_plot,
                            data_barplot,
                            max_height = 1.5,
                            below_center_y = 0.2,
                            max_width = 0.75,
                            max_value_display = NULL,
                            color = "#228a03",
                            color_background = "grey50",
                            border_color = "black",
                            alpha_background = 0,
                            value_col = "rbo",
                            cluster_col = "cluster",
                            scatter_1_col = "umapscvi_1",
                            scatter_2_col= "umapscvi_2"){
  
  for(i in 1:nrow(data_barplot)){
    dataplot = data_barplot[i,]
    dataplot$dummy = "1"
    if(is.null(max_value_display)){
      max_value_display = max(data_barplot[,value_col]) 
    }
    add = data.frame(cluster_col = data_barplot[i,cluster_col],value_col=max_value_display)
    colnames(add) =c(cluster_col,value_col)
    dataplot = dplyr::bind_rows(dataplot,add)
    temp_barplot = ggplot(dataplot,aes_string(x=cluster_col,y=value_col,fill="dummy",alpha="dummy"))+geom_col(show.legend = FALSE,position = "identity",color = border_color)+
      scale_fill_manual(values = color,na.value = color_background)+
      scale_alpha_manual(values = 1,na.value = alpha_background)+
      theme_bw()  + theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank(),legend.position="none",
                          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
                          panel.grid.minor=element_blank(),plot.background=element_blank())+
      xlab("")+ylab("")
    # add
    scatter_plot = scatter_plot +
      annotation_custom(ggplotGrob(temp_barplot),
                        xmin = data_barplot[i,scatter_1_col]-(max_width/2),
                        xmax = data_barplot[i,scatter_1_col]+(max_width/2),
                        ymin = data_barplot[i,scatter_2_col] - below_center_y,
                        ymax = data_barplot[i,scatter_2_col] - below_center_y + max_height)
    
  }
  
  return(scatter_plot)
}

##########
### custom_co_expression
##########

## with opacity: does not work well!

coexpression_opacity <- function(seurat_object,genes,reduction_name,colors=c("#6b0801","#1d6b01","#01066b"),color_seven = FALSE,alpha = 0.5,text_size=15,bg_color = "grey90",pt_size= 0.5) {
  # require(scales)
  
  gene_data = Seurat::FetchData(seurat_object,vars = genes)
  plotdata =  dplyr::bind_cols(seurat_object@reductions[[reduction_name]]@cell.embeddings,gene_data)
  # get column names of plotdata
  col_labels = colnames(plotdata)
  
  # default colors to use:
  if(color_seven){colors=as.character(palette.colors(palette = "Okabe-Ito")[2:8])}
  if(length(colors) < ncol(gene_data)){
    stop("Error: More genes than colors. Recommended are a maximum of three genes!") 
  }
  
  # alpha_values 
  # init with alpha 1
  # plotdata$alpha = 1
  # # get gene occurences
  # 
  # plotdata$alpha = rowSums(gene_data_binary)
  # # claculate alpha as 1 / occurences
  # plotdata$alpha = 1 / plotdata$n_expressed_genes
  # # set background alpha ( all zeros get alpha 1, the others 0)
  # plotdata$alpha_bg = plotdata$alpha
  # plotdata$alpha_bg[plotdata$alpha_bg <= 1] = 0
  # plotdata$alpha_bg[!is.infinite(plotdata$alpha)] = 1
  # 
  # # infinite get alpha zero
  # plotdata$alpha[is.infinite(plotdata$alpha)] = 0
  # 
  # 
  ##
  gene_data_binary = as.matrix(gene_data)
  gene_data_binary[gene_data_binary > 0 ] = 1
  gene_data_alpha = as.data.frame(t(apply(gene_data_binary,1,function(x){
    y = x / sum(x)
    y[is.nan(y)]=0
    return(y)
  })))
  gene_data_alpha_rowsums = rowSums(gene_data_alpha)
  gene_data_alpha_rowsums = abs(1-gene_data_alpha_rowsums)
  
  gene_data_alpha = as.matrix(gene_data_alpha)
  gene_data_alpha[ gene_data_alpha < 1 & gene_data_alpha > 0] = 0.9
  
  # make basic plot
  g = ggplot(plotdata, aes_string(x = col_labels[1], y = col_labels[2])) +
    geom_point( size = pt_size,color=bg_color,alpha=gene_data_alpha_rowsums)+
    theme_bw() + theme(panel.background=element_blank(),
                       panel.grid.major=element_blank(),
                       panel.grid.minor=element_blank(),
                       plot.background=element_blank(),
                       text = element_text(size=text_size))
  global_max = max(gene_data)
  # add genes 
  for(i in 1:ncol(plotdata[,3:ncol(plotdata)])){
    g = g+ ggnewscale::new_scale_color()
    g = g + geom_point(aes_string(color=col_labels[i+2]), size = pt_size,alpha = gene_data_alpha[,i]) +
      scale_color_gradientn(colors = c("white",colors[i]),limits = c(0, global_max))
  }
  g
  return(g)
}



