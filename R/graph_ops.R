constrainedEdgesMatrix <- function(wmat, chains, input_data) {
  ##
  ## Rules:
  ##  - cluster (node) cannot connect to itself
  ##  - a cluster with near-zero MCF cannot have children
  ##  - a cluster present in X multiple samples cannot connect to a cluster present in Y samples
  ##    if X < Y
  ##         - X < Y implies ...
  ##
  ##cluster.sample.presence <- apply(w, 1, function(x) which( x>= zero.thresh))
  samp_pres <- matchSamplePresence(wmat, chains, input_data)
  if (input_data$S == 1) {
    cluster.sample.presence <- lapply(samp_pres, function(x) which(x==1))
  } else {
    cluster.sample.presence <- apply(samp_pres, 1, function(x) which(x==1))
  }
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
      } else if (is.vector(cluster.sample.presence)) {
        from.samples <- cluster.sample.presence[i]
        to.samples <- cluster.sample.presence[j]
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

matchSamplePresence <- function(w_mat, chains, input_data) {
  map_z <- estimateClusterAssignments(chains$z_chain) 
  
  # pull one variant for each cluster
  single_var_clust <- map_z[match(unique(map_z$value), map_z$value), ] %>%
    mutate(mut_ind = as.numeric(gsub("z\\[|]", "", Parameter))) %>%
    arrange(value)
  
  samp_pres <- ifelse(input_data$y[single_var_clust$mut_ind, ] >=1, 1, 0)
  return(samp_pres)
}


## refactored base.admat
constrainedEdges <- function(wmat, chains, input_data) {
    am2 <- constrainedEdgesMatrix(wmat, chains, input_data)
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
  temp.possible.root.children <- filter(possible.edges, parent == "root")$child
  if (length(temp.possible.root.children) > 1) {
    temp.node <- sample(temp.possible.root.children, 1)
  } else {
    temp.node <- temp.possible.root.children
  }
  
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
  # stops if there is a cycle present in graph
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
  # update graph elements
  am <- updateGraphElements(am)
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

initializeGraphFromPost <- function(post_am, max.num.root.children=1, thresh=0.1) {
  constrained_am <- post_am %>%
    group_by(child) %>%
    mutate(max_post_for_child = max(posterior_prob)) %>%
    ungroup() %>%
    mutate(possible_edge = (max_post_for_child-posterior_prob) <= thresh) %>%
    mutate(connected = 0) %>%
    select(edge, parent, child, possible_edge, connected)
  am <- randAdmatUnchecked(constrained_am, max.num.root.children)
  while (!validGraph(am)) {
    am <- randAdmatUnchecked(constrained_am, max.num.root.children)
  }
  return(am)
}

initializeGraphFromPost2 <- function(post_am, max.num.root.children=1) {
  K <- length(unique(post_am$child))
  thresh <- (1/K)/2
  constrained_am <- post_am %>%
    mutate(possible_edge = posterior_prob >= thresh) %>%
    mutate(connected = 0) %>%
    select(edge, parent, child, possible_edge, connected)
  am <- randAdmatUnchecked(constrained_am, max.num.root.children)
  while (!validGraph(am)) {
    am <- randAdmatUnchecked(constrained_am, max.num.root.children)
  }
  return(am)
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
    move_set <- possible_moves_list[is_valid]
    #ix <- sample(seq_len(length(move_set)), 1)
    ix <- tryCatch(sample(seq_len(length(move_set)), 1), error=function(e) NULL)
    if(is.null(ix)) {
      #print("no moves :(")
      #print(a, n=100)
      #saveRDS(list(am = a, move_set), "/mnt/disk005/data/projects/pictograph/scripts/method-comparison/pictograph/log-sample-error/test.rds")
      return(a)
    } else {
      astar <- addEdge(a, move_set[[ix]])
      return(astar)
    }
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
  
  par(mar=c(0,0,0,0)+.1)
  # vertex color generation
  gradient <- c("#ffffff", "#f7f8be", "#f3d9f7", "#f4c5ff", "#f9d496", "#f5b8b8", "#c4bdee",
                "#d2e9a4", "#acda89", "#87cb6f", "#7ad0c7", "#7db1d5", "#f09eec", "#fCf073", 
                "#21f8ff", "#4de4e3", "#39cdff", "#9daa57", "#54a071", "#3e643b", "#217572", 
                "#2088AF", "#1F9AEC", "#EC63E1", "#AB3893", "#9D1540", "#8a21e9", "#9822ff", 
                "#804dff", "#6878ff", "#51a2ff", "#1e5a9e", "#690d44", "#6f1fbd", "#541c90", 
                "#1d3a76", "#1c1a4f", "#eeeeee", "#c9c9c9", "#a3a3a3", "#7e7e7e", "#585858", 
                "#a6bdaa", "#a48f7f", "#83613f")
  # warning message if num_vertices > max_vertices
  max_vertices <- 45
  if (length(V(ig)) > max_vertices) {
    warning("Number of clusters exceeds maximum numbers of distinct colors.")
  }
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
              vertex.color = gradient[1:length(V(ig))], vertex.size=35, vertex.frame.color = "#000000", vertex.label.family = "Helvetica",
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
  post_am <- post_am %>% 
    mutate(child = as.numeric(post_am$child))
  if(!is.factor(post_am$parent)) {
    post_am <- post_am %>%
      mutate(parent = factor(parent, levels = c("root", 1:max(post_am$parent))))
  }
  post_am %>%
    select(parent, child, posterior_prob) %>% 
    tidyr::spread(child, posterior_prob) %>%
    select(-parent) %>%
    as.matrix()
}

filterAdmat <- function(admat, filter1 = TRUE, filter1.threshold = 0.1,
                        filter2 = TRUE, filter2.threshold = 0.1) {
  # filter1 filters columns (am wide format) for edges with posterior prob > (max(column) - filter1.threshold)
  # filter2 filters entire matrix for prob > filter2.threshold
  
  if (filter1) {
    admat <- apply(admat, 2, function(x) ifelse(x > (max(x)-filter1.threshold), x, 0))
  } 
  
  if (filter2) {
    admat[admat <= filter2.threshold] <- 0
  }
  
  return(admat)
}

prepPostAmForGraphing <- function(post_am) {
  post_am_mat <- toWidePostAm(post_am)
  
  # add column for root
  post_am_mat <- cbind(0, post_am_mat) 
  colnames(post_am_mat)[1] <- "root"
  rownames(post_am_mat) <- colnames(post_am_mat)
  admat <- round(as.matrix(post_am_mat), 2)
  admat[is.na(admat)] <- 0 
  
  return(admat)
}

plotPosteriorAmLong <- function(post_am, filter1 = TRUE, filter1.threshold = 0.1,
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
  
  igraph::V(ig)$label.cex <- 0.5
  
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
              vertex.color = "white", vertex.label.family = "Helvetica",
              edge.arrow.size = 0.2, edge.arrow.width = 2,
              edge.width = igraph::E(ig)$weight*3)
}

plotPosteriorAmLongSim <- function(post_am, filter1 = TRUE, filter1.threshold = 0.1,
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
  
  igraph::V(ig)$label.cex <- 1
  
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      vertex.color = "white", vertex.label.family = "Helvetica",
                      vertex.size = 30,
                      edge.arrow.size = 0.6, edge.arrow.width = 2,
                      edge.width = igraph::E(ig)$weight*5,
                      asp = 0)
}

plotPosteriorAmLong2 <- function(post_am, cluster_key_genes_tb,
                                filter1 = TRUE, filter1.threshold = 0.1,
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
  
  igraph::V(ig)$label.cex <- 0.5
  vertex_colors <- c("white", ifelse(cluster_key_genes_tb$contains_key_gene, "lightblue", "white"))
  igraph::V(ig)$color <- vertex_colors
  # highlight single sample vertices
  igraph::V(ig)$frame.color <- c("black", ifelse(cluster_key_genes_tb$single_sample, "#E69F00", "black"))
  
  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      #vertex.color = "white", 
                      vertex.label.family = "Helvetica",
                      edge.arrow.size = 0.2, edge.arrow.width = 2,
                      edge.width = igraph::E(ig)$weight*3)
}

plotPosteriorAmLong3 <- function(post_am, cluster_key_genes_tb,
                                 filter1 = TRUE, filter1.threshold = 0.1,
                                 filter2 = TRUE, filter2.threshold = 0.1) {
  # don't plot edges with low freq (<0.25)
  # filter1 filters columns (am wide format) for edges with posterior prob > (max(column) - filter1.threshold)
  admat <- prepPostAmForGraphing(post_am)
  
  # filter edges of low freq
  admat <- filterAdmat(admat, filter1 = filter1, filter1.threshold = filter1.threshold,
                       filter2 = filter2, filter2.threshold = filter2.threshold)
  
  # filter out edges < 0.2 prob
  admat[admat < 0.2] <- 0
  ig <- igraph::graph_from_adjacency_matrix(admat, mode = "directed", weighted = TRUE,
                                            diag = FALSE, add.row = TRUE) 
  # edges < 0.25 prob
  igraph::E(ig)$lty <- ifelse(igraph::E(ig)$weight < 0.25, 2, 1)
  
  # make edge black if only 1 edge to vertex
  e <- igraph::ends(ig, igraph::E(ig))
  numTo <- table(e[,2])
  edgeColors <- sapply(e[,2], function(x) ifelse(x %in% names(which(numTo==1)), "black", "darkgrey"))
  igraph::E(ig)$color <- edgeColors
  
  igraph::V(ig)$label.cex <- 0.5
  vertex_colors <- c("white", ifelse(cluster_key_genes_tb$contains_key_gene, "lightblue", "white"))
  igraph::V(ig)$color <- vertex_colors
  # highlight single sample vertices
  igraph::V(ig)$frame.color <- c("black", ifelse(cluster_key_genes_tb$single_sample, "#E69F00", "black"))

  par(mar=c(0,0,0,0)+.1)
  igraph::plot.igraph(ig, layout = igraph::layout_as_tree(ig),
                      #vertex.color = "white", 
                      vertex.label.family = "Helvetica",
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

# getParents <- function(am_long, node, prev_parents = c()) {
#   temp_parent <- filter(am_long, connected == 1, child == node)$parent
#   if (temp_parent == "root") {
#     return(c(prev_parents, "root"))
#   } else {
#     return(getParents(am_long, temp_parent, c(prev_parents, temp_parent)))
#   }
# }

getBTColumn <- function(am_long, node) {
  BT_column <- rep(0, length(unique(am_long$child)))
  path <- getPathFromRoot(node, am_long)
  BT_column[path] <- 1
  return(BT_column)
}

amToBT <- function(am_long) {
  # returns binary perfect phylogeny matrix
  # where columns are clones and rows are mutation clusters
  # clones have the same numeric ID as the last mutation cluster on their path from the root
  children <- unique(am_long$child)
  children_num <- sort(as.numeric(children))
  
  BT_cols <- sapply(children_num, function(node) getBTColumn(am_long, node))
  
  return(BT_cols)
}
