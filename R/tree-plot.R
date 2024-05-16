#' Plot ensemble tree
#' 
#' @export
#' @param trees list of tibbles of edges, each with columns edge, parent, child
plotEnsembleTree <- function(trees, palette=viridis::viridis) {
  am_chain <- lapply(trees, edgesToAmLong)
  post_am <- getPosteriorAmLong(am_chain)
  plotPosteriorAmLong(post_am, colorScheme(trees[[1]], palette))
}

plotPosteriorAmLong <- function(post_am, v_color, filter1 = TRUE, filter1.threshold = 0.1,
                                filter2 = TRUE, filter2.threshold = 0.1) {
  # filter1 filters columns (am wide format) for edges with posterior prob > (max(column) - filter1.threshold)
  admat <- prepPostAmForGraphing(post_am)
  
  # filter edges of low freq
  admat <- filterAdmat(admat, filter1 = filter1, filter1.threshold = filter1.threshold,
                       filter2 = filter2, filter2.threshold = filter2.threshold)
  
  ig <- igraph::graph_from_adjacency_matrix(admat, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE) 
  
  igraph::E(ig)$lty <- ifelse(igraph::E(ig)$weight < 0.25, 2, 1)
  
  # make edge black if only 1 edge to vertex
  e <- igraph::ends(ig, igraph::E(ig))
  numTo <- table(e[,2])
  edgeColors <- sapply(e[,2], function(x) ifelse(x %in% names(which(numTo==1)), "black", "darkgrey"))
  igraph::E(ig)$color <- edgeColors
  
  igraph::V(ig)$label.cex <- 1.5
  
  igraph::V(ig)$color <- as.list(v_color %>% arrange(match(v_sorted, names(V(ig)))) %>% select(colors))$colors
  
  par(mar=c(0,0,0,0)+.1)
  
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      vertex.label.family = "Helvetica", vertex.size=20,
                      edge.arrow.size = 0.5, edge.arrow.width = 2,
                      edge.width = igraph::E(ig)$weight*3)
}

#' Plot single tree 
#' 
#' @export
#' @param edges tibble of edges with columns edge, parent, child
plotTree <- function(edges, palette=viridis::viridis) {
  plotGraph(edgesToAmLong(edges), colorScheme(edges, palette))
}

#' generate colors for each vertice
#' @export
colorScheme <- function(edges, palette=viridis::viridis) {
  v_sorted = sort(unique(c(edges$parent, edges$child)))
  v_sorted = c(sort(as.integer(v_sorted[!v_sorted=='root'])), "root")
  # root_idx <- which(v_sorted=="root")
  colors <- c(palette(length(v_sorted)-1), "white")
  v_color <- tibble(v_sorted, colors)
  return(v_color)
}

plotGraph <- function(am.long, v_color){
  # make sure am.long is sorted by parent and child
  am.long <- mutate(am.long, child = as.numeric(am.long$child)) %>%
    arrange(parent, child)
  am.long <- mutate(am.long, child = as.character(am.long$child))
  
  # change to wide format and plot
  am <- toWide(am.long)
  rownames(am) <- c("root", colnames(am))
  am <- cbind(root=0, am) ## add column for root
  colnames(am) <- rownames(am)
  
  am[is.na(am)] <- 0
  
  ig <- igraph::graph_from_adjacency_matrix(am, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE) 
  V(ig)$color <- as.list(v_color %>% arrange(match(v_sorted, names(V(ig)))) %>% select(colors))$colors
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      vertex.size=24, vertex.frame.color = "#000000", vertex.label.cex = 1.5,
                      vertex.label.family = "Helvetica", vertex.label.color = "#000000",
                      edge.arrow.size = 0.5, edge.arrow.width = 2)
}