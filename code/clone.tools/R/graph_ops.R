constrainedEdgesMatrix <- function(wmat, zero.thresh=0.01) {
  ##
  ## Rules:
  ##  - cluster (node) cannot connect to itself
  ##  - a cluster with near-zero MCF cannot have children
  ##  - a cluster present in X multiple samples cannot connect to a cluster present in Y samples
  ##    if X < Y
  ##         - X < Y implies ...
  ##
  ##cluster.sample.presence <- apply(w, 1, function(x) which( x>= zero.thresh))
  cluster.sample.presence <- apply(wmat, 1, function(x) which(x>=zero.thresh))
  K <- nrow(wmat)
  S <- ncol(wmat)
  admat <- matrix(T, K, K)
  
  for(i in 1:K){
    for(j in 1:K){
      
      if (is.matrix(cluster.sample.presence)) {
        from.samples <- cluster.sample.presence[, i]
        to.samples <- cluster.sample.presence[, j]
      } else if (is.list(cluster.sample.presence)) {
        from.samples <- cluster.sample.presence[[i]]
        to.samples <- cluster.sample.presence[[j]]
      }
      
      if (setequal(from.samples, to.samples)) next()
      if(length(from.samples) < length(to.samples)) {
        admat[i, j] <- F
        next()
      }
      if (all(to.samples %in% from.samples)) next()
      admat[i, j] <- F 
    }            
  }
  diag(admat) <- F
  am2 <- rbind(T, admat)
  dimnames(am2) <- list(c("root", 1:K), 1:K)
  return(am2)
}


## refactored base.admat
constrainedEdges <- function(wmat, zero.thresh=0.01) {
    am2 <- constrainedEdgesMatrix(wmat, zero.thresh)
    am2.long <- as_tibble(am2) %>%
        mutate(parent=rownames(am2)) %>%
        pivot_longer(-parent,
                     names_to="child",
                     values_to="possible_edge") %>%
        filter(parent != child) %>%
        unite("edge", c("parent", "child"), sep="->",
              remove=FALSE) %>%
        mutate(parent=factor(parent, levels=unique(parent))) %>%
        mutate(connected=0)
    am2.long
}

calcConstrianedTreeSpace <- function(mcf_matrix, zero.thresh = 0.01) {
  # input: 
  #     - mcf_matrix = matrix of cell fraction values where rows are clusters, columns are samples
  #     - zero.thresh = minimum cell fraction to be considered "present" in sample (default = 0.01)
  # output: number of possible trees, given constraints
  ce <- constrainedEdgesMatrix(mcf_matrix, zero.thresh)
  possible_from_edges_per_node <- colSums(ce)
  tree_space <- prod(possible_from_edges_per_node)
  return(tree_space)
}

calcUnconstrainedTreeSpace <- function(mcf_matrix) {
  # input: 
  #     - mcf_matrix = matrix of cell fraction values where rows are clusters, columns are samples
  # output: number of possible trees, no constraints
  K <- nrow(mcf_matrix)
  tree_space <- K^K
  return(tree_space)
}

getAllNodes <- function(am.long) {
  # returns vector of all nodes in graph
  am.long$parent <- as.character(am.long$parent)
  unique(c(am.long$parent, am.long$child))
}

reversedConnection <- function(am) {
    connections <- setNames(am$connected, am$edge)
    reversed_connections <- connections[am$reverse_edge] %>%
        "["(!is.na(.))
    reversed <- setNames(rep(0, nrow(am)), am$reverse_edge)
    reversed[names(reversed_connections)] <- reversed_connections
    reversed
}

isBidirectional <- function(am) {
    am %>%
        mutate(bi_directional=(reverse_edge %in% edge) &
                   connected==1 &
                   reversed_connected == 1) %>%
        pull(bi_directional)
}

updateGraphElements <- function(am) {
    am %>%
        mutate(parent=factor(parent, levels=unique(parent))) %>%
        mutate(reversed_connected=reversedConnection(.)) %>%
        mutate(bi_directional=isBidirectional(.)) %>%
        mutate(root_connected=isRootConnected(.))
}


reversedEdges <- function(am) {
    am2 <- am %>%
        filter(!is.na(connected)) %>%
        unite("reverse_edge", c("child", "parent"), sep="->",
              remove=FALSE)
    am2
}

getEdgeName <- function(from, to) {
  paste0(from, "->", to)
}

numNodesConnectedToRoot <- function(am.long) {
  sum(am.long[am.long$parent == "root", ]$connected)
}

randAdmatUnchecked <- function(am.long, max.num.root.children) {
  # input: blank am.long (from either constrainedEdges or toLong(initEmptyAdmatFromK(K))
  #     - $connected = 0 
  #     - $possible_edge = T/F if from constrainedEdges
  # output: random graph (not necessarily valid)
  blank <- am.long # save copy of original am.long
  parent_levels <- levels(blank$parent)
  
  am.long$parent <- as.character(am.long$parent)
  
  all.nodes <- getAllNodes(am.long)
  node.pool <- all.nodes[all.nodes != "root"] # nodes left to connect in graph
  # possible edges may be limited by constraints
  if ("possible_edge" %in% colnames(am.long)) {
    possible.edges <- am.long %>%
      filter(possible_edge == TRUE)
    has_constraints <- TRUE
  } else {
    possible.edges <- am.long
    has_constraints <- FALSE
  }
  parent.pool <- unique(possible.edges$parent) # possible parent nodes 
  
  # choose node to connect to root
  # select "to" node from parent.pool to prevent getting stuck if max.num.root.children == 1
  temp.node <- sample(parent.pool[parent.pool != "root"], 1)
  am.long[am.long$edge == getEdgeName("root", temp.node), ]$connected <- 1
  node.pool <- node.pool[node.pool != temp.node]
  
  # connect nodes that are left 
  while(length(node.pool) > 0) {
  #for (n in node.pool) {
    if (length(node.pool) > 1) {
      n <- sample(node.pool, 1)
    } else {
      n <- node.pool
    }
    
    # all possible edges to node n
    # check if there are constraints present 
    if (has_constraints) {
      temp.possible.edges <- am.long %>%
        filter(possible_edge == T)
    } else {
      temp.possible.edges <- am.long
    }
    # can't connect to root if max.num.root.children quota is satisfied
    if (numNodesConnectedToRoot(am.long) >= max.num.root.children) {
      temp.possible.edges <- temp.possible.edges %>%
        filter(child == n, parent != "root")
    } else {
      temp.possible.edges <- temp.possible.edges %>%
        filter(child == n)
    }
    
    # choose edge to connect -- randomly sample if more than 1 possible edge
    if(nrow(temp.possible.edges) > 1) {
      temp.edge <- temp.possible.edges[sample(nrow(temp.possible.edges), 1), ]$edge
    } else if (nrow(temp.possible.edges) == 1) {
      temp.edge <- temp.possible.edges$edge
    } else {
      # if no possible edges, start over
      return(randAdmatUnchecked(blank, max.num.root.children))
    }
    
    # connect edge 
    am.long[am.long$edge == temp.edge, ]$connected <- 1
    
    node.pool <- node.pool[node.pool != n]
  }
  
  am.long <- am.long %>%
    mutate(parent = factor(am.long$parent, levels = parent_levels))
  
  am.long <- reversedEdges(am.long) %>%
    mutate(reversed_connected=reversedConnection(.),
           bi_directional=NA,
           root_connected=NA)    
  am.long <- updateGraphElements(am.long)
  
  return(am.long)
}

randAdmat <- function(am.long, max.num.root.children) {
  # input: blank am.long (from either constrainedEdges or toLong(initEmptyAdmatFromK(K))
  #     - $connected = 0 
  #     - $possible_edge = T/F if from constrainedEdges
  # output: random graph
  blank <- am.long
  has_constraints <- ifelse("possible_edge" %in% colnames(am.long), 
                            ifelse(any(!am.long$possible_edge), T, F),
                            F)
  
  am.long$parent <- as.character(am.long$parent)
  
  all.nodes <- getAllNodes(am.long)
  node.pool <- all.nodes[all.nodes != "root"] # nodes left to connect in graph
  possible.edges <- am.long %>% 
    filter(!is.na(connected))
  parent.pool <- unique(possible.edges$parent) # possible parent nodes 
  
  # choose node to connect to root
  # select "to" node from parent.pool to prevent getting stuck if max.num.root.children == 1
  temp.node <- sample(parent.pool[parent.pool != "root"], 1)
  am.long[am.long$edge == getEdgeName("root", temp.node), ]$connected <- 1
  node.pool <- node.pool[node.pool != temp.node]
  from.nodes <- c("root", temp.node)
  
  while(length(node.pool) > 0) {
    
    # remove "root" from possible parents if max.num.root.children quota satisfied
    if(numNodesConnectedToRoot(am.long) < max.num.root.children) {
      from.nodes.pool <- from.nodes
    } else {
      from.nodes.pool <- from.nodes[-1]
    }
    
    if (length(from.nodes.pool) == 1) {
      temp.from <- from.nodes.pool
    } else {
      temp.from <- sample(from.nodes.pool, 1)
    } 
    
    # remove possible "to" nodes based on NA constraints 
    temp.to.pool <- am.long %>% 
      filter(parent==temp.from) %>% 
      filter(!is.na(connected))
    # remove "to" nodes if not in node.pool
    temp.to.pool <- temp.to.pool$child
    temp.to.pool <- temp.to.pool[temp.to.pool %in% node.pool]
    # sample possible children nodes to be "to" node and connect edge 
    if (length(temp.to.pool) > 1) {
      temp.to <- sample(temp.to.pool, 1)
    } else if (length(temp.to.pool) == 1) {
      temp.to <- temp.to.pool
    } else {
      # remove temp.from from from.nodes
    }
    
    am.long[am.long$edge == getEdgeName(temp.from, temp.to), ]$connected <- 1
    
    # add temp.to to possible from.nodes if it is a possible parent
    if (temp.to %in% parent.pool) {
      from.nodes <- c(from.nodes, temp.to)
    }
    
    node.pool <- node.pool[node.pool != temp.to]
  }
  
  
  
  am.long <- reversedEdges(am.long) %>%
    mutate(reversed_connected=reversedConnection(.),
           bi_directional=NA,
           root_connected=NA)    
  am.long <- updateGraphElements(am.long)
  am.long
}

isParentConnected <- function(am) {
    am %>%
        mutate(parent=factor(parent, levels=unique(parent))) %>%
        group_by(parent) %>%
        summarize(n=sum(connected)) %>%
        pull(n) > 0
}

isRootConnected <- function(am) isParentConnected(am)[1]

isDirected <- function(am) !any(am$bi_directional)

isFullyConnected <- function(am.long) {
    # checks if graph (am.long format) is fully connected
    all_nodes <- getAllNodes(am.long)
    nodes_in_main_tree <- bfsLong(am.long)
    length(all_nodes) == length(nodes_in_main_tree)
}

containsCycle <- function(am.long) {
    # returns nodes in main tree (connected to root) including "root" 
    # starting at root
    am.long$parent <- as.character(am.long$parent)
    children <- am.long[(am.long$parent == "root" & am.long$connected == 1), ]$child
    nodes <- c("root", children)
    
    while(length(children) > 0) {
        c <- children[1]
        temp.children <- am.long[(am.long$parent == c & am.long$connected == 1), ]$child
        children <- c(children, temp.children)
        if (any(temp.children %in% nodes)) return(TRUE)
        nodes <- c(nodes, temp.children)
        
        children <- children[-1]
    }
    FALSE
}

validGraph <- function(am) {
    isDirected(am) &&
        isRootConnected(am) &&
            !containsCycle(am) &&
                isFullyConnected(am)
}

bfsLong <- function(am.long) {
  # returns vector of nodes in main tree (connected to root) including "root" 
  # starting at root
  am.long$parent <- as.character(am.long$parent)
  children <- am.long[(am.long$parent == "root" & am.long$connected == 1), ]$child
  nodes <- c("root", children)
  
  while(length(children) > 0) {
    c <- children[1]
    temp.children <- am.long[(am.long$parent == c & am.long$connected == 1), ]$child
    children <- c(children, temp.children)
    if (any(temp.children %in% nodes)) stop("graph has cycle")
    nodes <- c(nodes, temp.children)
    children <- children[-1]
  }
  return(nodes)
}

addEdge <- function(am, new_edge) {
  c <- new_edge$child
  # disconnect existing edge connecting to child 
  am[which(am$child == c & am$connected == 1), ]$connected <- 0
  # connect new edge
  am[which(am$edge == new_edge$edge), ]$connected <- 1
  return(am)
}

initializeGraph <- function(mcf, max.num.root.children=1, zero.thresh=0.01){
    # clusters <- seq_len(nrow(mcf))
    # nsamp <- ncol(mcf)
    # samples <- seq_len(nsamp)
    # nclust <- length(clusters)
    # mcf.long <- tibble(cluster_id=as.character(rep(clusters, nsamp)),
    #                    sample_id=as.character(rep(samples, each=nclust)),
    #                    mean=as.numeric(mcf))    
    am.long <- constrainedEdges(mcf, zero.thresh=zero.thresh)
    am.long2 <- randAdmatUnchecked(am.long, max.num.root.children)
    while (!validGraph(am.long2)) {
      am.long2 <- randAdmatUnchecked(am.long, max.num.root.children)
    }
    return(am.long2)
}


toWide <- function(am.long){
    am.long$child <- as.numeric(am.long$child)
    am.long %>% select(parent, child, connected) %>%
        tidyr::spread(child, connected) %>%
        select(-parent) %>%
        as.matrix()
}

toLong <- function(am) {
  am.long <- as_tibble(am) %>%
    mutate(parent=rownames(am)) %>%
    pivot_longer(-parent,
                 names_to="child",
                 values_to="connected") %>%
    filter(parent != child) %>%
    unite("edge", c("parent", "child"), sep="->",
          remove=FALSE) %>%
    mutate(parent=factor(parent, levels=unique(parent)))   
  return(am.long)
}

isMoveValid <- function(a, possible_move, max.num.root.children) {
  # a = am.long format of current graph
  # possible_move = a row in possible_moves tibble (am.long format)
  # max.num.root.children = maximum number of nodes allowed to be connected to root
  # returns TRUE or FALSE
  astar <- addEdge(a, possible_move)
  is_valid <- validGraph(astar) & (numNodesConnectedToRoot(astar) <= max.num.root.children)
  return(is_valid)
}

sampleNewEdge <- function(a, max.num.root.children, mc.cores=1){
    ## a move is connecting a new edge and disconnecting the pre-existing edge connected to the new edge's child
    possible_moves <- filter(a, connected==0, possible_edge==T)
    possible_moves_list <- possible_moves %>%
      group_by(edge) %>%
      group_split()
    is_valid <- unlist(parallel::mclapply(possible_moves_list, function(x) isMoveValid(a, x, max.num.root.children),
                                          mc.cores = mc.cores))
    move_set <- possible_moves[is_valid, ]
    ix <- sample(seq_len(nrow(move_set)), 1)
    astar <- addEdge(a, move_set[ix, ])
    return(astar)
}

initEmptyAdmatFromK <- function(K) {
  admat <- matrix(0, K, K)
  diag(admat) <- NA
  am2 <- rbind(0, admat)
  dimnames(am2) <- list(c("root", 1:K), 1:K)
  return(am2)
}

generateRandomGraphFromK <- function(K, max.num.root.children) {
  # input: number of mutation clusters, K
  # output: mutation tree; adjacency matrix
  am.long <- toLong(initEmptyAdmatFromK(K))
  rand.am.long <- randAdmat(am.long, max.num.root.children)
  if (!validGraph(rand.am.long)) warning("graph is not valid")
  return(rand.am.long)
}

plotGraph <- function(am.long){
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
  
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
              vertex.color = "white", vertex.label.family = "Helvetica",
              edge.arrow.size = 0.2, edge.arrow.width = 2)
}

getPosteriorAmLong <- function(am_chain) {
  # input: chain from tree MCMC of trees in am.long format
  # output: posterior am.long 
  num_trees <- length(am_chain)
  
  combined_am_chain <- am_chain %>%
    bind_rows
  post_am <- combined_am_chain %>%
    group_by(edge) %>%
    mutate(posterior_prob = sum(connected) / num_trees) %>%
    ungroup() %>%
    select(edge, parent, child, posterior_prob) %>%
    distinct()
  post_am
}

toWidePostAm <- function(post_am) {
  post_am %>% 
    mutate(child = as.numeric(post_am$child)) %>%
    select(parent, child, posterior_prob) %>%
    tidyr::spread(child, posterior_prob) %>%
    select(-parent) %>%
    as.matrix()
}

plotPosteriorAmLong <- function(post_am, filter1 = TRUE, filter1.threshold = 0.1) {
  # filter1 filters columns (am wide format) for edges with posterior prob > (max(column) - filter1.threshold)
  post_am_mat <- toWidePostAm(post_am)
  
  # add column for root
  post_am_mat <- cbind(0, post_am_mat) 
  colnames(post_am_mat)[1] <- "root"
  rownames(post_am_mat) <- colnames(post_am_mat)
  admat <- as.matrix(post_am_mat)
  admat[is.na(admat)] <- 0 

  # filter edges
  if (filter1) {
    #thresh <- apply(admat, 2, max) - filter1.threshold
    admat <- apply(admat, 2, function(x) ifelse(x > (max(x)-filter1.threshold), x, 0))
  } 
  
  ig <- igraph::graph_from_adjacency_matrix(admat, mode = "directed", weighted = TRUE,
                                    diag = FALSE, add.row = TRUE) 
  
  igraph::E(ig)$lty <- ifelse(igraph::E(ig)$weight < 0.25, 2, 1)
  
  # make edge black if only 1 edge to vertex
  e <- igraph::ends(ig, igraph::E(ig))
  numTo <- table(e[,2])
  edgeColors <- sapply(e[,2], function(x) ifelse(x %in% names(which(numTo==1)), "black", "darkgrey"))
  igraph::E(ig)$color <- edgeColors
  
  igraph::V(ig)$label.cex <- 0.5
  
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
              vertex.color = "white", vertex.label.family = "Helvetica",
              edge.arrow.size = 0.2, edge.arrow.width = 2,
              edge.width = igraph::E(ig)$weight*3)
}

labelMAPEdgesFromPostAM <- function(post_am) {
  post_am %>%
    group_by(child) %>%
    mutate(map_edge = posterior_prob == max(posterior_prob)) %>%
    ungroup()
}

getMAPGraphFromPostAM <- function(post_am) {
  map_am <- labelMAPEdgesFromPostAM(post_am) %>%
    mutate(connected = ifelse(map_edge, 1, 0))
  return(map_am)
}

edgeTibbleToAmLong <- function(edge_tb, root = 0) {
  K <- length(unique(c(edge_tb$From, edge_tb$To))) - 1
  am_long <- toLong(initEmptyAdmatFromK(K))
  edge_tb$From <- as.character(edge_tb$From)
  edge_tb[edge_tb == as.character(root)] <- "root"
  for (i in seq_len(nrow(edge_tb))) {
    temp_edge <- getEdgeName(edge_tb$From[i], edge_tb$To[i])
    am_long$connected[which(am_long$edge == temp_edge)] <- 1
  }
  return(am_long)
}
