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

bfsLong2 <- function(am.long) {
  # returns vector of nodes in main tree (connected to root) including "root" 
  # starting at root
  # does not stop if there is a cycle in graph 
  am.long$parent <- as.character(am.long$parent)
  children <- am.long[(am.long$parent == "root" & am.long$connected == 1), ]$child
  nodes <- c("root", children)
  
  while(length(children) > 0) {
    c <- children[1]
    temp.children <- am.long[(am.long$parent == c & am.long$connected == 1), ]$child
    
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

grow <- function(tree_T, all_vertices) {
  
  if (length(verticesInGraph(tree_T)) == length(all_vertices) & nrow(tree_T) == (length(all_vertices)-1)) {
    assign("all_spanning_trees", c(all_spanning_trees, list(tree_T)), envir = .GlobalEnv)
    assign("i", i+1, envir = .GlobalEnv)
    return(tree_T)
    
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
      grow(tree_T, all_vertices)
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

