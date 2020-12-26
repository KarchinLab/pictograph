pop <- function(edges_tb, tb_name) {
  assign(tb_name, edges_tb[-1, ], envir = .GlobalEnv)
  return(edges_tb[1, ])
}

verticesInGraph <- function(tb) {
  unique(c(tb$parent, tb$child))
}

bridgeTest <- function(graph_G, edge_e) {
  node_to_check <- edge_e$child
  !(node_to_check %in% graph_G$child)
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

grow <- function(tree_T, all_vertices, w, thresh=0.1) {
  
  if (length(verticesInGraph(tree_T)) == length(all_vertices) & nrow(tree_T) == (length(all_vertices)-1)) {
    assign("all_spanning_trees", c(all_spanning_trees, list(tree_T)), envir = .GlobalEnv)
    
    if (satisfiesSumCondition(tree_T, w, thresh)) {
      assign("filtered_trees", c(filtered_trees, list(tree_T)), envir = .GlobalEnv)
    }
    
  } else {
    FF <- tibble(parent = character(), child = character())
      
    bridge <- FALSE
    while(!bridge) {
      # new tree edge
      if (nrow(F_tb) == 0) stop("F_tb is empty")
      edge_e <- pop(F_tb, "F_tb")
      v <- edge_e$child
      tree_T <- rbind(tree_T, edge_e)
      
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
      grow(tree_T, all_vertices, w, thresh)
      tree_L <- all_spanning_trees[[length(all_spanning_trees)]]
      
      # restore F
      # pop each edge (v,w), w not in T, from F
      not_in_T <- all_vertices[!all_vertices %in% verticesInGraph(tree_T)]
      if (length(not_in_T) > 0 & nrow(F_tb) > 0) {
        edges_to_remove_9 <- paste0(v, "->", not_in_T)
        assign("F_tb", filter(F_tb, !edge %in% edges_to_remove_9), envir = .GlobalEnv)
      }
      # restore each edge (w,v), w in T, in F
      assign("F_tb", rbind(removed_edges, F_tb), envir = .GlobalEnv)
     
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

prepareGraphForGabowMyers <- function(w, zero.thresh=0.01) {
  # input: matrix of cellular prevalences (rows = clusters, columns = samples)
  # output: tibble of possible edges with columns edge, parent, child
  graph_G <- constrainedEdges(w, zero.thresh=zero.thresh) %>%
    filter(possible_edge == TRUE) %>%
    mutate(parent = as.character(parent)) %>%
    select(edge, parent, child)
  return(graph_G)
}

enumerateSpanningTrees <- function(graph_G, w, sum_filter_thresh=0.1) {
  # all_spanning_trees must be set as an empty list, global variable, before function is called
  # graph_G must be set as global variable before function is called
  all_spanning_trees <- assign("all_spanning_trees", list(), envir = .GlobalEnv)
  filtered_trees <- assign("filtered_trees", list(), envir = .GlobalEnv)
  F_tb <- assign("F_tb", filter(graph_G, parent == "root"), envir = .GlobalEnv)
  all_vertices <- verticesInGraph(graph_G)
  tree_T <- tibble(parent = character(), child = character())
  
  
  
  grow(tree_T, all_vertices, w, sum_filter_thresh)
  
  #return(all_spanning_trees)
#   return(list(all_spanning_trees = all_spanning_trees,
#               filtered_trees = filtered_trees))
}

satisfiesSumCondition <- function(edges, w, thresh = 0.1) {
  # returns TRUE if sum condition is not violated with given threshold (default 0.1)
  
  edges$parent <- as.character(edges$parent)
  all_parents <- unique(edges$parent)
  
  for (p in all_parents) {
    # get parent CCF
    if (p == "root") {
      parent_ccf <- rep(1, ncol(w))
    } else {
      parent_ccf <- w[as.numeric(p), ]
    }
    
    # get children CCF (sum if more than 1 child)
    children <- as.numeric(filter(edges, parent == p)$child)
    if (length(children) > 1) {
      children_ccf <- colSums(w[children, ])
    } else {
      children_ccf <- w[children, ]
    }
    
    diff <- children_ccf - parent_ccf
    if (any(diff > thresh)) return(FALSE)
  }
  
  # sum condition is never violated, return TRUE
  return(TRUE)
}

filterEdgesBasedOnCCFs <- function(graph_G, w, thresh = 0.2) {
  check_edges_logical <- apply(graph_G, 1, function(edge) checkEdge(edge, w, thresh))
  filtered_graph_G <- graph_G[check_edges_logical, ]
  return(filtered_graph_G)
}

checkEdge <- function(edge, w, thresh = 0.2) {
  # returns TRUE if satisfies lineage precedence with given threshold
  # returns FALSE if violates i.e. child_ccf - parent_ccf > thresh in any sample
  # edge is in the format c(edge_name, parent, child)
  
  # in case of factors
  p <- as.character(edge[2])
  c <- as.character(edge[3])
  
  if (p == "root") {
    parent_ccfs <- rep(1, ncol(w))
  } else {
    parent_ccfs <- w[as.numeric(p), ]
  }
  child_ccfs <- w[as.numeric(c), ]
  
  diff <- child_ccfs - parent_ccfs
  if (any(diff > thresh)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

calcTreeSpaceUpperBound <- function(w, zero.thresh = 0.01, 
                                    lineage.precedence.filter = 0.1) {
  # calculates upper bound of tree space given CCF matrix (w)
  graph_G_pre <- prepareGraphForGabowMyers(w, zero.thresh = zero.thresh)
  graph_G <- filterEdgesBasedOnCCFs(graph_G_pre, w, thresh = lineage.precedence.filter)
  num_edges_to_each_child <- graph_G %>% 
    group_by(child) %>%
    summarize(n = n()) %>%
    pull(n)
  upper_bound <- prod(num_edges_to_each_child)
  return(upper_bound)
}

calcTreeSpaceUpperBound2 <- function(edges) {
  # calculates upper bound of tree space given tibble of edges (rows = edge, parent, child)
  num_edges_to_each_child <- edges %>% 
    group_by(child) %>%
    summarize(n = n()) %>%
    pull(n)
  upper_bound <- prod(num_edges_to_each_child)
  return(upper_bound)
}

splitGraphG <- function(graph_G, num_trees_per_run = 100000) {
  # returns list of graph_Gs
  # figure out which nodes for which to keep edges constant 
  # each split graph_G should have tree space upper bound < num_trees_per_run 
  num_edges_per_child <- graph_G %>%
    group_by(child) %>%
    summarize(n = n()) %>%
    arrange(desc(n))
  
  # split graph_G tree space is larger than num_trees_per_run 
  if (prod(num_edges_per_child$n) > num_trees_per_run) {
    nodes_to_hold <- c()
    for (i in seq_len(nrow(num_edges_per_child))) {
      upper_bound_per_run <- prod(num_edges_per_child$n[i:nrow(num_edges_per_child)])
      if (upper_bound_per_run <= num_trees_per_run) break
      nodes_to_hold <- c(nodes_to_hold, num_edges_per_child$child[i])
    }
    
    common_edges <- filter(graph_G, ! child %in% nodes_to_hold)
    edges_to_split <- filter(graph_G, child %in% nodes_to_hold) %>%
      group_by(child) %>%
      group_split()
    edges_list <- lapply(edges_to_split, function(x) x$edge) 
    all_combinations <- expand.grid(edges_list, stringsAsFactors = F)
    split_graph_Gs <- apply(all_combinations, 1, function(x) rbind(filter(graph_G, edge %in% x), common_edges))
    return(split_graph_Gs)
    
  } else {
    return(list(graph_G))
  }
}
