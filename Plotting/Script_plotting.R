library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)


dat.list <- split(peak.assignments.scores.nodes,peak.assignments.scores.nodes$node)

plot.list <- list()
for(i in 1:13){

    dat1 <- dat.list[[i]]

    plot.list[[i]] <- ggplot(dat1, aes(x=default, y=tau)) +
    stat_density_2d(geom = "polygon", contour = TRUE,
                  aes(fill = after_stat(level)), colour = "black",
                  bins = 5)+
    geom_point(aes(alpha=ultraconserved, color=hit)) + 
    geom_hline(yintercept = 0.7, linetype="dashed", color = "black") +
    geom_hline(yintercept = 0.3, linetype="dashed", color = "black") +
    scale_fill_distiller(palette = "Greens", direction = 1) +
    theme_bw() +
    scale_color_manual(values=c("Broadly detected" = "blue", "Variable" = "red"), na.value = "black") +
    scale_alpha_discrete(range=c(0.1, 1)) + 
    xlim(0, 1) +
    ylim(0, 1) +
    facet_wrap(~node) +
    theme(legend.position="none") +
    labs(x="Cons. score", y="Tau score", color="")
                
}

p <- plot_grid(plotlist = plot.list, align = "hv",nrow = 2, axis = "l")


library(igraph)
igraph_obj <- graph_from_edgelist(el)
layout_nicely(igraph_obj) -> layout
rownames(layout) <- V(igraph_obj)$name
node_df <- data.frame(
  "name"=rownames(layout),
  "X"=layout[,1],
  "Y"=layout[,2],
  "size_pval"=pval_vec,
  "size_or"=or_vec,
  "group"=group_vec,
  stringsAsFactors = F
)


library(ggraph)
library(tidygraph)
ggraph_obj <- as_tbl_graph(igraph_obj)

ggraph_obj <- ggraph_obj %>% 
  tidygraph::activate(nodes) %>%
  left_join(node_df, by = c("name" = "name")) %>% 
  tidygraph::activate(edges) %>%
  mutate(from_name = (.N()$name[from])) %>%
  mutate(to_name = (.N()$name[to])) %>%
  mutate(pair_name = paste(from_name, to_name, sep=":")) %>%
  left_join(edge_df, by=c("pair_name"="pair_name")) %>%
  tidygraph::activate(nodes)



p1 <- ggraph(ggraph_obj, x=X, y=Y) +
  geom_edge_diagonal(aes(alpha=or, color=factor(direction)),
                     width=0.5, arrow = arrow(length = unit(1,"mm"))) +
  scale_edge_alpha_continuous(range=c(0.1,0.8), guide = "none") +
  scale_edge_color_manual(values = c('-1'='#7FB3D5', '1'='#EC7063')) +
  geom_node_point(aes(size = size_or, fill = group),
                  shape = 21, color='darkgrey') +
  scale_fill_manual(values = group_cols) +
  scale_size_continuous(range = c(1,5), trans = "sqrt") +
  geom_node_point(aes(size = size_or, filter = name%in%top_terms),
                  shape = 1, color='#303030') +
  geom_node_label(aes(label=name, filter = group=="Cell_type_label"),
                 size=3, repel=T, max.overlaps = 13) +
  geom_node_text(aes(label=name, filter = name %in% top_terms),
                  min.segment.length = unit(0, 'lines'), size=3, repel=T, max.overlaps = 99999) +
  theme_void()

pdf("Plot_ggraph_DAR_enriched_GO_terms.pdf")
p1
dev.off()