# this script is to visualize enriched GO terms and their parent terms in GO tree structure

# Author: Qianhui Yu
# Date 2024-May-22

library(igraph)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggrepel)

# step 0: build the basic graph
all_go_pairs <- readRDS("Dat_directed_GO_term_pairs.rds")
igraph_obj <- readRDS("Res_directed_GO_term_igraph_obj.rds")
# focus on BP
root_terms <- "GO:0008150" 
# extract the node ID from an igraph object
graph_node_id <- V(igraph_obj)$name
# load enriched GO term accession numbers
go_term_id <- readRDS("Res_enriched_GO_term_id.rds")
selected_terms <- intersect(go_term_id, graph_node_id)
# get the nodes along the simple paths from each selected terms to the root terms
node_list <- lapply(selected_terms, function(x){
  paths <- all_simple_paths(graph=igraph_obj, from=x, to=root_terms, mode="out")
  unique(names(unlist(paths))) 
})
kept_nodes <- unique(unlist(node_list))
# get the edges between the selected nodes
kept_pairs <- all_go_pairs[which(all_go_pairs$from %in% kept_nodes & all_go_pairs$to %in% kept_nodes),]
dim(kept_pairs)
saveRDS(kept_pairs, file="Res_kept_pairs.rds")
# build the graphs from the edges
graph_obj <- graph_from_data_frame(kept_pairs)
saveRDS(graph_obj,file="Res_graph_obj_before_trimming.rds")
# layout in the tree structure
layout <- layout.reingold.tilford(graph_obj, flip.y=T, root=which(V(graph_obj)$name %in% root_terms), mode="all")
rownames(layout) <- V(graph_obj)$name
saveRDS(layout, file="Res_layout_before_trimming.rds")

# step 1: node clustering
# for clustering purpose, remove the root terms
cl_kept_nodes <- setdiff(kept_nodes, root_terms) 
cl_kept_pairs <- all_go_pairs[which(all_go_pairs$from %in% cl_kept_nodes & all_go_pairs$to %in% cl_kept_nodes),]
no_root_graph_obj <- graph_from_data_frame(cl_kept_pairs, directed=FALSE)
adj_mat <- as_adjacency_matrix(no_root_graph_obj)
cl_res <- Seurat::FindClusters(adj_mat, resolution=0.5, algorithm=2)
length(unique(cl_res[,1]))
cl_res_vec <- setNames(c(as.numeric(cl_res[,1]), length(unique(cl_res[,1]))+1), c(rownames(cl_res), root_terms))
saveRDS(cl_res_vec, file="Res_node_clustering_result.rds")

df <- readRDS("/home/yuq22/ihb-intestine-evo/Annotation/GO/Gene_ontology_relationship/Dat_GO_ID_and_names.rds")
node_df <- data.frame(
  "name"=rownames(layout),
  "X"=layout[,1],
  "layer"=layout[,2],
  "Y"=layout[,2]+rnorm(nrow(layout), sd=0.1),
  "desc"=df$name[match(rownames(layout), df$term)],
  "group"=cl_res_vec[rownames(layout)],
  stringsAsFactors = F
)
saveRDS(node_df, file="Res_enriched_GO_term_related_node_df_before_trimming.rds")

# step 2: subclustering
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
saveRDS(dist_table, file="Res_shortest_distance_between_nodes_inside_subgroup_and_subgroup_top_term.rds")

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
saveRDS(node_df, file="Res_node_df_after_subclustering.rds")
n1 <- length(unique(node_df$subgroup_top_id))
n1

# step 3: merge the small subgroups to the nearest subgroup
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
saveRDS(node_df, file="Res_node_df_after_merging_small_groups.rds")

# step 4: trim the graph based on the merged subgroups
# only the top terms could keep edges with other term group, all other terms only keep the edges with the same group
kept_pairs$to_keep <- TRUE
top_terms <- unique(node_df$merged_subgroup_top_id)
for(i in sort(unique(node_df$merged_subgroup_top_id))){
  group_terms <- rownames(node_df)[which(node_df$merged_subgroup_top_id==i)]
  idx <- which(kept_pairs$from%in%setdiff(group_terms, top_terms) & !kept_pairs$to%in%group_terms)
  kept_pairs$to_keep[idx] <- FALSE
}
saveRDS(kept_pairs, file="Res_kept_pairs.rds")
kept_pairs_after_trimming <- kept_pairs[which(kept_pairs$to_keep), 1:2]

# step 5: visualize the trimmed graph
graph_obj <- graph_from_data_frame(kept_pairs_after_trimming)
saveRDS(graph_obj,file="Res_graph_obj_after_trimming.rds")
layout <- layout.reingold.tilford(graph_obj, flip.y=T, root=which(V(graph_obj)$name %in% root_terms), mode="all")
rownames(layout) <- V(graph_obj)$name
saveRDS(layout, file="Res_layout_after_trimming.rds")

# update the node_df with the updated graph coordinates 
shared_terms <- intersect(rownames(layout), rownames(node_df))
node_df <- node_df[shared_terms,]
node_df$X <- layout[shared_terms,1]
node_df$Y <- layout[shared_terms,2]+rnorm(length(shared_terms), sd=0.1)
saveRDS(node_df, file="Res_node_df_trimmed.rds")

ggraph_obj <- as_tbl_graph(graph_obj)
ggraph_obj <- ggraph_obj %>% 
  tidygraph::activate(nodes) %>%
  left_join(node_df, by = c("name" = "name")) %>% 
  tidygraph::activate(edges) %>%
  mutate(from_name = (.N()$name[from])) %>%
  mutate(to_name = (.N()$name[to])) %>%
  mutate(pair_name = paste(from_name, to_name, sep=":")) %>%
  tidygraph::activate(nodes)
saveRDS(ggraph_obj, file="Res_ggraph_obj_after_trimming_subclustering.rds")

# get the description of the top terms
aa <- unique(node_df[which(node_df$desc == node_df$merged_subgroup_desc),c("merged_subgroup_idx","merged_subgroup_desc")])

# define colors for the subgroups
ids_colors <- c("#F9E79F", "#4DAA99",
                "#A4D371","#138D75", "#F1C40F",
                "#08519C","#3C7AB6","#8E44AD",
                "#BB8FCE","#F0B27A",
                "#85C1E9",
                "#756BB1","#B5CAE5","#7BCAA4",
                "#758AB1","#CA6778", 
                "#322A84","#DCCB7C","#9E0142",
                "#B06754","#F781BF",
                "#41B6C4","#EDF8B1","#909090")                
names(ids_colors) <- paste0("G",seq(length(ids_colors)))
saveRDS(ids_colors, file="Res_ids_colors.rds")

# plot the graph without highlighting specific nodes
p1 <- ggraph(ggraph_obj, x=X, y=Y) +
  geom_edge_diagonal(color="#b0b0b0", width=0.05, alpha=0.5) +
  geom_node_point(aes(x=X, y=Y, fill=merged_subgroup_idx),
                  shape = 21, color="#303030", stroke=0.2, size=6) +
  scale_fill_manual(labels = sub("G", "", paste(sort(names(ids_colors)), aa$merged_subgroup_desc[match(sort(names(ids_colors)),aa$merged_subgroup_idx)], sep="-")), values = ids_colors) +
  geom_node_label(aes(label=merged_subgroup_idx, filter = node_df$name%in%rownames(aa), color=merged_subgroup_idx),
                 size=8, repel=T, max.overlaps = 999, show.legend=FALSE) +
  scale_color_manual(values = ids_colors) +
  theme_void()+
  theme(legend.position="bottom")
p1

pdf("Plot_willow_plot_of_enriched_GO_terms.pdf", width=20, height=15)
p1
dev.off()
