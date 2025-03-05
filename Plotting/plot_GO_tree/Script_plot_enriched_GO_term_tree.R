# this script is to visualize enriched GO terms and their parent terms in GO tree structure

# Author: Qianhui Yu
# Date 2024-May-22

build_graph <- function(all_go_pair_path="/home/yuq22/ihb-intestine-evo/common_script/intestine_evolution_script_repository/plot_GO_tree/Dat_directed_GO_term_pairs.rds", 
                        igraph_obj_path="/home/yuq22/ihb-intestine-evo/common_script/intestine_evolution_script_repository/plot_GO_tree/Res_directed_GO_term_igraph_obj.rds",
                        root_terms="GO:0008150", # focus on BP by default
                        go_term_id=NULL,
                        GO_ID_and_name_path="/home/yuq22/ihb-intestine-evo/common_script/intestine_evolution_script_repository/plot_GO_tree/Dat_GO_ID_and_names.rds"){
    
    library(igraph)
    library(ggplot2)
    library(ggraph)
    library(tidygraph)
    library(dplyr)
    library(ggrepel)
    
    if(is.null(go_term_id)){
        print("Please provide the enriched GO term accession numbers")
        return(NULL)
    }

    if(root_terms=="GO:0008150"){
        print("Focus on biological process (BP)")
    }else if(root_terms=="GO:0003674"){
        print("Focus on molecular function (MF)")  
    }else if(root_terms=="GO:0005575"){ 
        print("Focus on cellular component (CC)")    
    }

    all_go_pairs <- readRDS(all_go_pair_path)
    igraph_obj <- readRDS(igraph_obj_path)
    graph_node_id <- V(igraph_obj)$name

    # step 0: build the basic graph
    print("Step 0/6: build the basic graph")
    # focus on the enriched GO terms and their parent terms which are in the annotation file
    selected_terms <- intersect(go_term_id, graph_node_id)
    print(paste(length(selected_terms), "(", round(length(selected_terms)/length(go_term_id)*100), "%) selected terms retained"))
    # get the nodes along the simple paths from each selected terms to the root terms
    node_list <- lapply(selected_terms, function(x){
      paths <- all_simple_paths(graph=igraph_obj, from=x, to=root_terms, mode="out")
      unique(names(unlist(paths))) 
    })
    kept_nodes <- unique(unlist(node_list))
    # get the edges between the selected nodes
    kept_pairs <- all_go_pairs[which(all_go_pairs$from %in% kept_nodes & all_go_pairs$to %in% kept_nodes),]
    print(paste("numbers of kept edge:", nrow(kept_pairs)))
    
    # build the graphs from the edges
    graph_obj <- graph_from_data_frame(kept_pairs)
    print("Get graph object before trimming")
    
    # layout in the tree structure
    layout <- layout.reingold.tilford(graph_obj, flip.y=T, root=which(V(graph_obj)$name %in% root_terms), mode="all")
    rownames(layout) <- V(graph_obj)$name
    print("Get graph object before trimming")

    output <- list("kept_pairs"=kept_pairs, "graph_obj"=graph_obj, "layout"=layout, "kept_nodes"=kept_nodes)

    # step 1: node clustering
    print("Step 1/6: node clustering")
    # for clustering purpose, remove the root terms
    cl_kept_nodes <- setdiff(kept_nodes, root_terms) 
    cl_kept_pairs <- all_go_pairs[which(all_go_pairs$from %in% cl_kept_nodes & all_go_pairs$to %in% cl_kept_nodes),]
    no_root_graph_obj <- graph_from_data_frame(cl_kept_pairs, directed=FALSE)
    adj_mat <- as_adjacency_matrix(no_root_graph_obj)
    cl_res <- Seurat::FindClusters(adj_mat, resolution=0.5, algorithm=2)
    length(unique(cl_res[,1]))
    cl_res_vec <- setNames(c(as.numeric(cl_res[,1]), length(unique(cl_res[,1]))+1), c(rownames(cl_res), root_terms))
    output[["cl_res_vec"]] <- cl_res_vec
    
    df <- readRDS(GO_ID_and_name_path)
    node_df <- data.frame(
      "name"=rownames(layout),
      "X"=layout[,1],
      "layer"=layout[,2],
      "Y"=layout[,2]+rnorm(nrow(layout), sd=0.1),
      "desc"=df$name[match(rownames(layout), df$term)],
      "group"=cl_res_vec[rownames(layout)],
      stringsAsFactors = F
    )
    output[["node_df_before_trimming"]] <- node_df

    # step 2: subclustering
    print("Step 2/6: node subclustering")
    # because usually there are more than 1 top terms on the same layer 
    # assign the nodes to the nearest top term 
    node_df$is_subgroup_top_id <- FALSE
    #id_vec <- c()
    for(i in sort(unique(node_df$group))){
        idx <- which(node_df$group==i)
        group_df <- node_df[idx,] 
        max_idx <- which(group_df$layer==max(group_df$layer))
        node_df[rownames(group_df[max_idx,]),"is_subgroup_top_id"] <- TRUE
    }

    # get the shortest distance from each subgroup inside node to the subgroup top node
    dist_table <- distances(
      graph=graph_obj, 
      v=node_df$name[which(!node_df$is_subgroup_top_id)],
      to=node_df$name[which(node_df$is_subgroup_top_id)], 
      mode="out")
    output[["dist_table"]] <- dist_table 

    # get the shortest distance from each node to a top term node
    min_dist_vec <- apply(dist_table, 1, min)

    # determine subgroup top term based on shortest distance
    node_df$subgroup_top_id <- NA
    node_df$subgroup_top_id[which(node_df$is_subgroup_top_id)] <- node_df$name[which(node_df$is_subgroup_top_id)]
    for(id in rownames(dist_table)){
      min_dist_id <- colnames(dist_table)[which(dist_table[id,]==min_dist_vec[id])]
      if(length(min_dist_id)==1){
        node_df[id,"subgroup_top_id"] <- min_dist_id
      }else{
        x_dist <- abs(node_df[id,"X"] - node_df[min_dist_id,"X"])
        node_df[id,"subgroup_top_id"] <- min_dist_id[which.min(x_dist)]
      }
    }
    # get subgroup top term description
    subgroup_desc <- setNames(node_df$desc[which(node_df$is_subgroup_top_id)], node_df$name[which(node_df$is_subgroup_top_id)])
    node_df$subgroup_desc <- subgroup_desc[node_df$subgroup_top_id]
    output[["node_df_after_subclustering"]] <- node_df
    
    n1 <- length(unique(node_df$subgroup_top_id))
    print(paste(n1, "subgroups after subclustering"))

    # step 3: merge the small subgroups to the nearest subgroup
    print("Step 3/6: merge the small subgroups to the nearest subgroup")
    # get the size of revised subgroups
    size <- table(node_df$subgroup_top_id)
    size_desc <- sort(table(node_df$subgroup_desc))
    # get the subgroup with size smaller than 50
    size_cutoff <- 50
    small_subgroup <- setdiff(names(size)[which(size<size_cutoff)], root_terms)
    node_df$merged_subgroup_top_id <- node_df$subgroup_top_id
    node_df$merged_subgroup_desc <- node_df$subgroup_desc
    while(length(small_subgroup)>0){
        # get the shortest distance from small subgroup to the root terms
        dist_table_small_group <- distances(
          graph=graph_obj, 
          v=small_subgroup,
          to=unique(c(root_terms,node_df$subgroup_top_id)), 
          mode="out")
        dist_to_root <- sort(dist_table_small_group[,root_terms], decreasing=T)
        ids <- names(dist_to_root)[which(dist_to_root>1)]
        if(length(ids)>0){
          print(paste(length(ids),"to merge"))
            for(x in ids){
                # get the nearest subgroup
                other_top_terms <- setdiff(colnames(dist_table_small_group), x)
                nearest_subgroup <- other_top_terms[which(dist_table_small_group[x,other_top_terms]==min(dist_table_small_group[x,other_top_terms]))]
                if(length(nearest_subgroup)>1){
                  x_dist <- abs(node_df[x,"X"] - node_df[nearest_subgroup,"X"])
                  nearest_subgroup <- nearest_subgroup[which.min(x_dist)]
                }
                # get the nodes and its child nodes
                nodes_to_merge <- unique(c(node_df$name[which(node_df$subgroup_top_id==x)],node_df$name[which(node_df$merged_subgroup_top_id==x)]))
                # assign the node and its child nodes to the nearest subgroup
                node_df[nodes_to_merge,"merged_subgroup_top_id"] <- nearest_subgroup
                node_df[nodes_to_merge,"merged_subgroup_desc"] <- subgroup_desc[nearest_subgroup]
            }
            # get the size of revised subgroups
            size <- table(node_df$merged_subgroup_top_id)
            # get the small groups after merging
            small_subgroup <- names(size)[which(size<size_cutoff)]
        }else{
            break
        }
    }
    # assign cluster index
    desc <- sort(unique(node_df$merged_subgroup_desc))
    non_root_desc <- setdiff(desc, node_df[root_terms, "desc"])
    idx <- c(setNames(paste0("G", seq(length(non_root_desc))), non_root_desc), setNames(paste0("G", length(desc)),node_df[root_terms, "desc"]))
    node_df$merged_subgroup_idx <- idx[node_df$merged_subgroup_desc]
    output[["node_df_after_merging"]] <- node_df
    
    # step 4: trim the graph based on the merged subgroups
    # only the top terms could keep edges with other term group, all other terms only keep the edges with the same group
    print("Step 4/6: trim the graph")
    kept_pairs$to_keep <- TRUE
    top_terms <- unique(node_df$merged_subgroup_top_id)
    for(i in sort(unique(node_df$merged_subgroup_top_id))){
      group_terms <- rownames(node_df)[which(node_df$merged_subgroup_top_id==i)]
      idx <- which(kept_pairs$from%in%setdiff(group_terms, top_terms) & !kept_pairs$to%in%group_terms)
      kept_pairs$to_keep[idx] <- FALSE
    }
    kept_pairs_after_trimming <- kept_pairs[which(kept_pairs$to_keep), 1:2]
    output[["kept_pairs_after_trimming"]] <-  kept_pairs_after_trimming

    # step 5: visualize the trimmed graph
    print("step 5/6: get ggraph object")
    graph_obj <- graph_from_data_frame(kept_pairs_after_trimming)
    saveRDS(graph_obj,file="Res_graph_obj_after_trimming.rds")
    layout <- layout.reingold.tilford(graph_obj, flip.y=T, root=which(V(graph_obj)$name %in% root_terms), mode="all")
    rownames(layout) <- V(graph_obj)$name
    saveRDS(layout, file="Res_layout_after_trimming.rds")
    output[["graph_obj_after_trmming"]]=graph_obj
    output[["layout_after_trimming"]]=layout

    # update the node_df with the updated graph coordinates 
    shared_terms <- intersect(rownames(layout), rownames(node_df))
    node_df <- node_df[shared_terms,]
    node_df$X <- layout[shared_terms,1]
    node_df$Y <- layout[shared_terms,2]+rnorm(length(shared_terms), sd=0.1)
    output[["final_node_df"]]=node_df

    ggraph_obj <- as_tbl_graph(graph_obj)
    ggraph_obj <- ggraph_obj %>% 
      tidygraph::activate(nodes) %>%
      left_join(node_df, by = c("name" = "name")) %>% 
      tidygraph::activate(edges) %>%
      mutate(from_name = (.N()$name[from])) %>%
      mutate(to_name = (.N()$name[to])) %>%
      mutate(pair_name = paste(from_name, to_name, sep=":")) %>%
      tidygraph::activate(nodes)
    output[["final_ggraph_obj"]]=ggraph_obj
    # get the description of the top terms
    aa <- unique(node_df[which(node_df$desc == node_df$merged_subgroup_desc),c("merged_subgroup_idx","merged_subgroup_desc")])
    cluster_num=nrow(aa)
    print(paste(cluster_num, "GO clusters"))
    output[["top_terms_desc"]]=aa

    
    # get colors
    ids_color_pool <- c("#F9E79F", "#4DAA99",
                "#A4D371","#138D75", "#F1C40F",
                "#08519C","#3C7AB6","#8E44AD",
                "#BB8FCE","#F0B27A",
                "#85C1E9",
                "#756BB1","#B5CAE5","#7BCAA4",
                "#758AB1","#CA6778", 
                "#322A84","#DCCB7C","#9E0142",
                "#B06754","#F781BF",
                "#41B6C4","#EDF8B1","#909090")

  prettyrainbow <- c("#DF4C42", "#AE2042", "#832D3A", "#8848A2","#663882","#832B8D",
                    "#A06E95","#AB99CF","#F588A1","#D31477","#EB91B6","#E70497",
                    "#EC8082","#F07C28","#F89D7E","#F9D108","#F0E515","#D1E069",
                    "#65BF46","#0FA243","#008B5C","#036C4D","#449F4A","#B0CC48",
                    "#6DC0AD","#9CD5D4","#0C88A7","#0B8BCB","#09BAC2")

  if(cluster_num<=length(ids_color_pool)){
    ids_colors <- ids_color_pool[1:cluster_num]
  }else{
    ids_colors <- colorRampPalette(prettyrainbow)(cluster_num)
  }
  names(ids_colors) <- paste0("G",seq(length(ids_colors)))
  output[["ids_colors"]]=ids_colors
  
  # get ggplot object
  print("step 6/6: get ggplot object")
  p1 <- ggraph(ggraph_obj, x=X, y=Y) +
  geom_edge_diagonal(color="#b0b0b0", width=0.05, alpha=0.5) +
  geom_node_point(aes(x=X, y=Y, fill=merged_subgroup_idx),
                  shape = 21, color="#303030", stroke=0.2, size=6) +
  scale_fill_manual(labels = sub("G", "", paste(sort(names(ids_colors)), aa$merged_subgroup_desc[match(sort(names(ids_colors)),aa$merged_subgroup_idx)], sep="-")), values = ids_colors) +
  geom_node_label(aes(label=merged_subgroup_idx, filter = node_df$name%in%rownames(aa), color=merged_subgroup_idx),
                 size=8, repel=T, max.overlaps = 999, show.legend=FALSE) +
  scale_color_manual(values = ids_colors) +
  theme_void()+
  theme(legend.position="bottom", legend.title=element_text(size=20), legend.text=element_text(size=20), legend.key.size = unit(1.5,"line"))+
  guides(fill=guide_legend(nrow=5,byrow=TRUE))
  output[["ggplot_obj"]]=p1
  return(output)

}
