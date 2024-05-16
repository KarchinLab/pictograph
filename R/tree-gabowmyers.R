#' Enumerate all spanning trees using modified Gabow-Myers wrapper
#' 
#' @export
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @param lineage_precedence_thresh maximum allowed violation of lineage precedence (default = 0.1)
#' @param sum_filter_thresh thresh maximum allowed violation of Sum Condition (default = 0.2)
generateAllTrees <- function(mcf, purity, lineage_precedence_thresh=0.1, sum_filter_thresh=0.2) {
  mcf_mat <- estimateMCFs(mcf)
  mcf_mat <- assign("mcf_mat", mcf_mat, envir = .GlobalEnv)
  graph_G_pre <- prepareGraph(mcf_mat, lineage_precedence_thresh)
  graph_G <- filterEdgesBasedOnCCFs(graph_G_pre, mcf_mat, thresh = lineage_precedence_thresh)
  graph_G <- assign("graph_G", graph_G, envir = .GlobalEnv)
  enumerateSpanningTreesModified(graph_G, mcf_mat, purity, sum_filter_thresh = sum_filter_thresh)
}

#' Create tibble of possible edges from CCF values based on w_mat only
#' 
#' @export
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @return graph_G tibble of possible edges with columns edge, parent, child
prepareGraph <- function(mcf_mat, thresh) {
  graph_pre <- data.frame(edge = character(), parent = character(), child = character())
  for (i in seq_len(nrow(mcf_mat))) {
    graph_pre <- graph_pre %>% add_row(edge = paste("root->", i, sep = ""), parent = "root", child = as.character(i))
    for (j in seq_len(nrow(mcf_mat))) {
      if (i!=j) {
        i_row = mcf_mat[i, ]
        j_row = mcf_mat[j, ]
        if (all(j_row-i_row >= -thresh)) {
          graph_pre <- graph_pre %>% add_row(edge = paste(j, "->", i, sep = ""), parent = as.character(j), child = as.character(i))
        }
      }
    }
  }
  return(graph_pre)
}

#' Filter possible edges based on lineage precedence 
#' 
#' @export
#' @param graph_G tibble of possible edges with columns edge, parent, child
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @param thresh maximum allowed violation of lineage precedence (default = 0.1)
filterEdgesBasedOnCCFs <- function(graph_G, mcf, thresh = 0.1) {
  check_edges_logical <- apply(graph_G, 1, function(edge) checkEdge(edge, mcf, thresh))
  filtered_graph_G <- graph_G[check_edges_logical, ]
  return(filtered_graph_G)
}

checkEdge <- function(edge, mcf, thresh = 0.2) {
  # returns TRUE if satisfies lineage precedence with given threshold
  # returns FALSE if violates i.e. child_ccf - parent_ccf > thresh in any sample
  # edge is in the format c(edge_name, parent, child)
  
  # in case of factors
  p <- as.character(edge[2])
  c <- as.character(edge[3])
  
  if (p == "root") {
    parent_ccfs <- rep(1, ncol(mcf))
  } else {
    parent_ccfs <- mcf[as.numeric(p), ]
  }
  child_ccfs <- mcf[as.numeric(c), ]
  
  diff <- child_ccfs - parent_ccfs
  if (any(diff > thresh)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

#' Enumerate all spanning trees using modified Gabow-Myers
#' 
#' @export
#' @param graph_G tibble of possible edges with columns edge, parent, child
#' @param w matrix of CCF values (rows = clusters, columns = samples)
#' @param sum_filter_thresh thresh maximum allowed violation of Sum Condition (default = 0.2)
enumerateSpanningTreesModified <- function(graph_G, mcf, purity, sum_filter_thresh=0.2) {
  # all_spanning_trees must be set as an empty list, global variable, before function is called
  # graph_G must be set as global variable before function is called
  all_spanning_trees <- assign("all_spanning_trees", list(), envir = .GlobalEnv)
  #filtered_trees <- assign("filtered_trees", list(), envir = .GlobalEnv)
  F_tb <- assign("F_tb", filter(graph_G, parent == "root"), envir = .GlobalEnv)
  all_vertices <- verticesInGraph(graph_G)
  tree_T <- tibble(parent = character(), child = character())
  
  growModified(tree_T, all_vertices, mcf, purity, sum_filter_thresh)
}

verticesInGraph <- function(tb) {
  unique(c(tb$parent, tb$child))
}

growModified <- function(tree_T, all_vertices, w, purity, sum_thresh=0.2) {
  
  if (length(verticesInGraph(tree_T)) == length(all_vertices) & nrow(tree_T) == (length(all_vertices)-1)) {
    assign("all_spanning_trees", c(all_spanning_trees, list(tree_T)), envir = .GlobalEnv)
    
  } else {
    FF <- tibble(parent = character(), child = character())
    
    bridge <- FALSE
    while(!bridge) {
      # new tree edge
      if (nrow(F_tb) == 0) stop("F_tb is empty")
      edge_e <- pop(F_tb, "F_tb")
      v <- edge_e$child
      tree_T <- rbind(tree_T, edge_e)
      
      # check if adding this node does not violate the constraint
      if (satisfiesSumCondition(tree_T, w, purity, sum_thresh)) {
        # update F
        ## push each edge (v,w), w not in T onto F
        in_T <- verticesInGraph(tree_T)
        temp_add_to_F <- filter(graph_G, parent == v, !(child %in% in_T))
        # temp_add_to_F
        assign("F_tb", rbind(temp_add_to_F, F_tb), envir = .GlobalEnv)
        
        ## remove each edge (w,v), w in T from F
        w_in_T <- verticesInGraph(tree_T)
        removed_edges <- filter(F_tb, parent %in% w_in_T, child == v)
        assign("F_tb", filter(F_tb, !edge %in% removed_edges$edge), envir = .GlobalEnv)
        
        # recurse
        growModified(tree_T, all_vertices, w, purity, sum_thresh)
        
        # restore F
        # pop each edge (v,w), w not in T, from F
        not_in_T <- all_vertices[!all_vertices %in% verticesInGraph(tree_T)]
        if (length(not_in_T) > 0 & nrow(F_tb) > 0) {
          edges_to_remove_9 <- paste0(v, "->", not_in_T)
          assign("F_tb", filter(F_tb, !edge %in% edges_to_remove_9), envir = .GlobalEnv)
        }
        # restore each edge (w,v), w in T, in F
        assign("F_tb", rbind(removed_edges, F_tb), envir = .GlobalEnv)
        
      }
      # delete e from T and from G, add e to FF
      tree_T <- tree_T[tree_T$edge != edge_e$edge, ]
      assign("graph_G",graph_G[graph_G$edge != edge_e$edge, ], envir = .GlobalEnv)
      FF <- rbind(edge_e, FF)
      
      # bridge test
      bridge <- bridgeTestBFS(graph_G, edge_e)
    }
    
    # pop each edge e from FF, push e onto F,and add e to G
    if (nrow(FF) > 0) {
      
      # pop and push all edges at once (same order)
      # assign("F_tb", rbind(FF, F_tb), envir = .GlobalEnv)
      assign("graph_G", rbind(FF, graph_G), envir = .GlobalEnv)
      
      # pop and push edges one by one (rev order in F)
      while (nrow(FF) > 0) {
        assign("F_tb", rbind(FF[1, ], F_tb), envir = .GlobalEnv)
        FF <- FF[-1, ]
      }
    }
  }
}

satisfiesSumCondition <- function(edges, w, purity, thresh = 0.2) {
  # returns TRUE if sum condition is not violated with given threshold (default 0.2)
  
  edges$parent <- as.character(edges$parent)
  all_parents <- unique(edges$parent)
  
  for (p in all_parents) {
    # get parent CCF
    if (p == "root") {
      # parent_ccf <- rep(1, ncol(w))
      parent_ccf <- purity
    } else {
      parent_ccf <- w[as.numeric(p), ]
    }
    
    # get children CCF (sum if more than 1 child)
    children <- as.numeric(filter(edges, parent == p)$child)
    if (length(children) > 1) {
      children_ccf <- colSums(w[children, ,drop=FALSE])
    } else {
      children_ccf <- w[children, ]
    }
    
    diff <- children_ccf - parent_ccf
    if (any(diff > thresh)) return(FALSE)
  }
  
  # sum condition is never violated, return TRUE
  return(TRUE)
}

pop <- function(edges_tb, tb_name) {
  assign(tb_name, edges_tb[-1, ], envir = .GlobalEnv)
  return(edges_tb[1, ])
}

bridgeTestBFS <- function(graph_G, edge_e) {
  node_to_check <- edge_e$child
  
  nodes_connected_to_root <- bfsLong2(graph_G)
  !(node_to_check %in% nodes_connected_to_root)
}

bfsLong2 <- function(graph_G) {
  # returns vector of nodes in main tree (connected to root) including "root" 
  # starting at root
  # does not stop if there is a cycle in graph 
  graph_G$parent <- as.character(graph_G$parent)
  children <- graph_G[(graph_G$parent == "root"), ]$child
  nodes <- c("root", children)
  
  while(length(children) > 0) {
    c <- children[1]
    temp.children <- graph_G[(graph_G$parent == c), ]$child
    
    # remove children already seen
    if (any(temp.children %in% nodes)) {
      temp.children <- temp.children[! temp.children %in% nodes]
    }
    
    children <- c(children, temp.children)
    
    nodes <- c(nodes, temp.children)
    children <- children[-1]
  }
  return(nodes)
}