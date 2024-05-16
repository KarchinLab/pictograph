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
  ix <- tryCatch(sample(seq_len(length(move_set)), 1), error=function(e) NULL)
  if(is.null(ix)) {
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

summarizeWChain <- function(w.chain) {
  # output: mcf_stats 
  mcf_stats <- w.chain %>%
    group_by(Parameter) %>%
    summarize(sd=sd(value),
              mean=mean(value))
  return(mcf_stats)
}

create.cpov <- function(mcf_stats, alpha=0.05, zero.thresh=0.001, mcf_matrix = NULL, restriction.val = 1) {
  cpov <- NA
  MCF <- NA
  
  ## if mcf_matrix is supplied, use that to create cpov
  if (is.null(mcf_matrix)) {
    cpov <- initializeAdjacencyMatrix(mcf_stats = mcf_stats, zero.thresh = zero.thresh)
    cpov[is.na(cpov)] <- restriction.val
    MCF <- mcfMatrix(mcf_stats)
  } else {
    cpov <- initializeAdjacencyMatrix(mcf_matrix = mcf_matrix, zero.thresh = zero.thresh)
    cpov[is.na(cpov)] <- restriction.val
    MCF <- mcf_matrix
  }
  
  sds <- mcfMatrix(mcf_stats, parameter="sd")
  ##S <- ncol(mcmc_w) # number of samples
  S <- numberSamples(mcf_stats)
  ##cpov <- cpov[-1, ]
  ## root can go to anyone -- all 0's (default base admat value)
  for (r in 2:nrow(cpov)) {
    for (c in 1:ncol(cpov)) {
      if (cpov[r,c] == restriction.val) next # skip restricted position
      from <- r-1 # 'from' cluster node
      to <- c # 'to' cluster node
      statistic <- 0
      pval <- 0
      for(s in seq_len(S)) {
        ##d <- mcmc_w[from,s] - mcmc_w[to,s]
        d <- MCF[from, s] - MCF[to, s]
        d_sd <- sqrt(sds[from, s]^2 + sds[to, s]^2)
        ##d_sd <- sqrt((mcmc_w_sd[from,s])^2 + (mcmc_w_sd[to,s])^2)
        I <- sum(d < 0)
        ## cumulative sum of the
        ## number of standard deviations for the difference in
        ## MCFs between 2 samples
        if (d == 0 || is.nan(d / d_sd)) {
          next
        } else {
          statistic <- statistic + (d / d_sd)^2 * I
        }
        
        for (k in 0:S) {
          pval <- pval + ((1 - pchisq(statistic, k)) *
                            choose(S, k) / (2^S))
        }
      }
      
      pval <- ifelse(is.na(pval), 0, pval)
      pval <- ifelse(is.nan(pval), 0, pval)
      ##
      ## edge seems to be based on this ad-hoc statistic, not
      ## the probability of the tree
      ##
      cpov[r,c] <- decide.ht(pval, alpha)
    }
  }
  cpov
}

decide.ht <- function(pval, alpha=0.05) {
  # 1 signals rejection event for null of i -> j
  if (pval <= alpha) return(1)
  else return(0)
}

calcTopologyCost <- function(am, cpov, am_format = "long") {
  TC <- 0
  
  if (am_format == "long") {
    am <- toWide(am)
  } 
  
  edges <- which(am == 1, arr.ind=TRUE)
  N <- nrow(edges)
  for (i in seq_len(N)) {
    TC <- TC + cpov[edges[i,1], edges[i,2]]
  }
  
  TC
}

getEdges <- function(am.long) {
  am.long %>%
    filter(connected == 1) %>%
    mutate(parent = as.character(parent))
}

getChildren <- function(am.long, node) {
  # returns vector of children nodes
  edges <- am.long %>%
    mutate(parent = as.character(parent)) %>%
    filter(connected == 1) %>%
    filter(parent == node)
  return(edges$child)
}

calcMassCost <- function(am, mcf_matrix, purity, am_format="long") {
  num_samples <- ncol(mcf_matrix)
  
  if (am_format == "long") {
    edges <- getEdges(am)
    
    parent_nodes <- unique(edges$parent)
    mass_cost <- rep(0, length(parent_nodes)) # mass cost of each parent node
    
    for (i in seq_len(length(parent_nodes))) {
      parent_node <- parent_nodes[i]
      
      # root MCF is purity instead of 1
      if (parent_node == "root") {
        # parent_w <- rep(1, num_samples) # 1 replaced by purity
        parent_w <- purity
      } else {
        parent_w <- mcf_matrix[as.numeric(parent_node), ,drop=FALSE]
      }
      
      kids <- getChildren(am, parent_node)
      if (length(kids) > 1) {
        children_w <- colSums(mcf_matrix[as.numeric(kids), ,drop=FALSE])
      } else {
        children_w <- mcf_matrix[as.numeric(kids), ,drop=FALSE]
      }
      
      mc_s <- ifelse(parent_w >= children_w, 0, children_w - parent_w)
      #mass_cost[i] <- sqrt(sum(mc_s^2))
      mass_cost[i] <- max(mc_s) # take max across samples instead of euclidean distance
    }
    return(sum(mass_cost))
    
  } else if (am_format == "wide") {
    num_children <- rowSums(am, na.rm = T)
    nodes <- which(num_children > 0, arr.ind = T) # not leaves
    mc_node <- rep(0, length(nodes))
    
    for (i in 1:length(nodes)) {
      node <- nodes[i]
      
      # root node: MCF = 1
      parent_w <- rep(1, ncol(mcf_matrix))
      # not root node: look up MCF in mcf_matrix
      if (node != 1) { 
        parent_w <- mcf_matrix[node-1,]
      }
      
      kids <- which(am[node,] == 1, arr.ind = T)
      if (num_children[node] > 1) {
        children_w <- colSums(mcf_matrix[kids, ])
      } else {
        children_w <- mcf_matrix[kids, ]
      }
      
      mc_s <- ifelse(parent_w >= children_w, 0, children_w - parent_w)
      mc_node[i] <- sqrt(sum(mc_s^2))
    }
    return(sum(mc_node))
  }
}

edgesToAmLong <- function(edges) {
  am_wide <- initEmptyAdmatFromK(length(unique(edges$child)))
  edges[edges$parent == "root", "parent"] <- "0"
  edges <- edges %>%
    mutate(parent = as.numeric(parent) + 1,
           child = as.numeric(child)) %>%
    select(parent, child)
  edges <- as.matrix(edges)
  for (r in 1:nrow(edges)) {
    am_wide[edges[r,1], edges[r,2]] <- 1
  }
  admat <- toLong(am_wide)
  admat <- reversedEdges(admat) %>%
    mutate(reversed_connected=reversedConnection(.),
           bi_directional=NA,
           root_connected=NA)    
  admat <- updateGraphElements(admat)
  return(admat)
}

calcTreeFitness <- function(admat, cpov, mcf_matrix, purity, am_format = "long", weight_mass = 1, weight_topology = 1, scaling_coeff=5) {
  # if only edges are given, change into long format
  if (am_format == "edges") {
    admat <- edgesToAmLong(admat)
    am_format <- "long"
  }
  
  TC <- calcTopologyCost(admat, cpov, am_format)
  MC <- calcMassCost(admat, mcf_matrix, purity, am_format)
  Z <- weight_topology * TC + weight_mass * MC
  fitness <- exp(-scaling_coeff * Z)
  fitness
}

satisfiesCCFSumProperties <- function(am_long, mcf_matrix, threshold = 0.2) {
  # threshold = max value that children CCFs can be larger than parent;
  #   i.e. function returns false if (sum of children CCFs) - parent CCF > threshold
  
  edges <- getEdges(am_long)
  parent_nodes <- unique(edges$parent)
  
  for (i in seq_len(length(parent_nodes))) {
    parent_node <- parent_nodes[i]
    
    # root CCF is 1
    if (parent_node == "root") {
      parent_w <- rep(1, ncol(mcf_matrix))
    } else {
      parent_w <- mcf_matrix[as.numeric(parent_node), ]
    }
    
    kids <- getChildren(am_long, parent_node)
    if (length(kids) > 1) {
      children_w <- colSums(mcf_matrix[as.numeric(kids), ])
    } else {
      children_w <- mcf_matrix[as.numeric(kids), ]
    }
    
    # return false if violates sum properties
    if (any(children_w - parent_w > threshold)) return(FALSE)
  }
  
  return(TRUE)
}

#' Calculate SCHISM fitness scores for trees
#' 
#' @export
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @param trees list of tibbles, where each tibble contains edges of a tree with columns edge, parent, child
calcTreeScores <- function(mcf_chain, trees, purity, mc.cores = 1) {
  # calculate mean and sd of mcf for each cluster in each sample
  mcf_stats <- summarizeWChain(mcf_chain) 
  
  # create cpov matrix
  # first: using sample presence to create a binary matrix; 0 is i can be a ancestor of j; 1 if not
  # second: for i,j pair that pass the sample presence test, calc stats of the difference between mcf of all 
  # samples; return binary matrix
  # problem: is the stats calculation in create.cpov correct? only using cluster so is actually POV instead 
  # of CPOV
  cpov <- create.cpov(mcf_stats)
  mcf_mat <- estimateMCFs(mcf_chain)
  
  # first calculate topology cost: sum of the cpov matrix over all edges
  # potential problem: mcf[i]=0.4->mcf[j]=0.6 has same weight as mcf[i]=0.5->mcf[j]=0.6
  # second calculate mass cost: take the max mass violation among all samples; is there a better strategy?
  # fitness is exp(-5*(topology cost + mass cost))
  schism_scores <- unlist(parallel:::mclapply(trees, 
                                              function(x) calcTreeFitness(x, cpov, mcf_mat, purity, am_format = "edges"),
                                              mc.cores = mc.cores))
  return(schism_scores)
}

#' Calculate SCHISM fitness scores for trees
#' 
#' @export
calculateTreeScoreMutations <- function(mcf_chain, data, icnTable, cncfTable, multiplicityTable, clusterAssingmentTable, purity, trees, restriction.val = 1, mc.cores = 8) {
  
  mcfMutations1 <- (data$y * (icnTable$icn*cncfTable+2-2*cncfTable) - data$n * (multiplicityTable$Multiplicity-1)*cncfTable) / data$n
  mcfMutations1[mcfMutations1>1] <-1
  mcfMutations1[mcfMutations1<0] <-0
  mcfMutations2 <- (data$n - 2 * data$y) / (data$y * icnTable$icn - 2 * data$y - data$n*multiplicityTable$Multiplicity + data$n)
  mcfMutations2[mcfMutations2>1] <-1
  mcfMutations2[mcfMutations2<0] <-0
  
  mcfMutations <- matrix(0, nrow = nrow(mcfMutations1), ncol = ncol(mcfMutations1))
  mcfMutations[data$is_cn == 1, ] <- mcfMutations2[data$is_cn == 1, ]
  mcfMutations[data$is_cn == 0, ] <- mcfMutations1[data$is_cn == 0, ]
  mcfMutations[is.nan(mcfMutations)] <- 0
  
  colnames(mcfMutations) <- seq_len(ncol(mcfMutations))
  
  mutation_chain <- as_tibble(mcfMutations) %>% rownames_to_column("Row") %>%
    pivot_longer(
      cols = -Row,
      names_to = "Column",
      values_to = "value"
    ) %>% mutate(Row = as.integer(Row))

  mutation_chain <- clusterAssingmentTable %>% 
    mutate(row = row_number()) %>% 
    inner_join(mutation_chain, by = c("row" = "Row")) %>% 
    mutate(Parameter = paste("mcf[", row, ",", gsub("V", "", Column), "]", sep = "")) %>% 
    select(Parameter, value, Cluster)
  
  mutation_stats <- mutation_chain %>% 
    group_by(Cluster) %>% 
    summarize(sd=sd(value),mean=mean(value))
  
  mutation_stats <- mutation_chain %>% 
    left_join(mutation_stats, by="Cluster") %>%
    select(Parameter, sd, value) %>%
    mutate(mean = value) %>%
    select(Parameter, sd, mean)
  
  # create pov for mutaiton pairs
  pov <- create.cpov(mutation_stats)
  pov <- pov[2:nrow(pov),]
  
  cpov <- initializeAdjacencyMatrix(mcf_stats = summarizeWChain(mcf_chain))
  cpov[is.na(cpov)] <- restriction.val
  
  for (r in 2:nrow(cpov)) {
    for (c in 1:ncol(cpov)) {
      if (cpov[r,c] == restriction.val) next
      from <- r-1 # 'from' cluster node
      to <- c # 'to' cluster node
      fromMutations <- which(clusterAssingmentTable$Cluster==from)
      toMutations <- which(clusterAssingmentTable$Cluster==to)
      
      totPov <- 0
      for (fromIdx in fromMutations) {
        for (toIdx in toMutations) {
          totPov = totPov + pov[fromIdx, toIdx]
        }
      }
      totPov = totPov / (length(fromMutations)*length(toMutations))
      cpov[r,c] = totPov
    }
  }
  # scores <- calcTreeScores(chains$mcf_chain, all_spanning_trees, purity)
  mcf_mat <- estimateMCFs(mcf_chain)
  
  # first calculate topology cost: sum of the cpov matrix over all edges
  # potential problem: mcf[i]=0.4->mcf[j]=0.6 has same weight as mcf[i]=0.5->mcf[j]=0.6
  # second calculate mass cost: take the max mass violation among all samples; is there a better strategy?
  # fitness is exp(-5*(topology cost + mass cost))
  schism_scores <- unlist(parallel:::mclapply(trees, 
                                              function(x) calcTreeFitness(x, cpov, mcf_mat, purity, am_format = "edges"),
                                              mc.cores = mc.cores))
}

rand.admat <- function(admat) {
  for(col in 1:ncol(admat)) {
    ind.0 <- which(admat[,col] == 0) # possible positions (0's)
    
    # place 1 if only 1 position available
    if (length(ind.0) == 1) {
      admat[ind.0, col] <- 1
      
    } else { # else randomly pick which position to place 1
      rand.ind <- sample(ind.0, size=1)
      admat[rand.ind,col] <- 1
    }
  }
  
  if (sum(admat[1, ], na.rm=TRUE) == 0) {
    # pick random cluster to be connected to root
    rand.k <- sample(1:ncol(admat), size = 1)
    curr.1 <- which(admat[, rand.k] == 1)
    admat[curr.1, rand.k] <- 0
    admat[1, rand.k] <- 1
  }
  
  while (is.bidirectional(admat)) {
    admat <- fix.bidirectional(admat)
  }
  
  admat
}

restrictions.from.post.admat <- function(post.admat, threshold=0.1) {
  # input: posterior adjacency matrix
  # output: adjacency matrix with restrictions (0's in allowed positions, NAs in restricted positions)
  # For each column, posible positions are the position with the max posterior 
  #     probability and those within the threshold of the max (> (max - 0.1))
  
  thresh <- apply(post.admat, 2, max) - threshold
  apply(post.admat, 2, function(x) ifelse(x > (max(x)-threshold), 0, NA))
}

initializeAdjacencyMatrix <- function(mcf_stats=NULL, mcf_matrix=NULL, zero.thresh=0.01) {
  if (!is.null(mcf_stats)) {
    MCF <- mcfMatrix(mcf_stats)
  } else if (!is.null(mcf_matrix)) {
    MCF <- mcf_matrix
  } else stop("must supply either mcf_stats or mcf_matrix")
  
  ##
  ## for each row (cluster) of the MCF matrix,
  ## list indices of samples for which the cluster is present
  ##
  K <- nrow(MCF)
  S <- ncol(MCF)
  # cluster.sample.presence should be a list
  MCF_list <- split(MCF, seq(nrow(MCF)))
  cluster.sample.presence <- lapply(MCF_list, 
                                    function(x) which(x > zero.thresh))
  all.samples <- seq_len(S)
  ## initialize adjacency matrix
  ##admat <- matrix(data=0, nrow=(1+K), ncol=K)
  admat <- matrix(data=0, K, K)
  ## can't go to self
  diag(admat) <- NA  
  for(from in seq_len(K)){
    for(to in seq_len(K)){
      ## can't go to self 
      if (from == to) next()
      ## hierarchy restraints
      from.samples <- cluster.sample.presence[[from]]
      to.samples <- cluster.sample.presence[[to]]
      ## no restraints if same sample presence
      if (identical(from.samples, to.samples)) next()
      ## restraint if # from.samples < # to.samples
      if(length(from.samples) < length(to.samples)) {
        admat[from, to] <- NA 
        next()
      }
      ## no restraints if to.samples is subset of from.samples
      if (all(to.samples %in% from.samples)) {
        next
      } else {
        admat[from, to] <- NA 
      }
    }
  }
  ## Add root
  ## can go from root to anyone
  admat <- rbind(0, admat)
  dimnames(admat) <- list(c("root", paste0("cluster", seq_len(K))),
                          paste0("cluster", seq_len(K)))
  admat
}

init.admat <- function(w, zero.thresh) {
  base <- base.admat(w, zero.thresh)
  rand.admat(base)
}

#' @export
mutate.admat <- function(admat, ncol.to.mutate) {
  ## choose a column(s) to mutate
  K <- ncol(admat)
  rand.ks <- sample(seq_len(K), size=ncol.to.mutate)
  ## mutate columns
  new.admat <- admat 
  for(k in rand.ks){
    
    temp.admat <- mutate.column(new.admat, k)
    # make sure admat is fully connected
    while (!is.fully.connected(temp.admat)) {
      temp.admat <- mutate.column(new.admat, k)
    }
    new.admat <- temp.admat
  }
  
  # replace root edge if missing
  if (sum(new.admat[1, ]) == 0) {
    otherCol <- sample(seq_len(K)[-rand.ks], size = 1)
    ind.1 <- which(admat[, otherCol] == 1)
    new.admat[1, otherCol] <- 1
    new.admat[ind.1, otherCol] <- 0
  }
  
  # fix bidirectional edge if present
  while (is.bidirectional(new.admat)) {
    new.admat <- fix.bidirectional(new.admat)
  }
  
  new.admat
}

#' @export
mutate.admat.2 <- function(admat, ncol.to.mutate) {
  
  new.admat <- mutate.n.columns(admat, ncol.to.mutate)
  
  # make sure admat is fully connected
  while (!is.fully.connected(new.admat)) {
    new.admat <- mutate.n.columns(admat, ncol.to.mutate)
  }
  
  new.admat
}

#' @export
mutate.n.columns <- function(admat, ncol.to.mutate) {
  K <- ncol(admat)
  # columns with more than 1 possible position
  columns.to.choose.from <- seq_len(K)[apply(admat, 2, function(x) sum(is.na(x)) < K)]
  # sample columns from those with more than 1 possible position
  rand.ks <- sample(columns.to.choose.from, size=ncol.to.mutate)
  for(k in rand.ks){
    admat <- mutate.column(admat, k)
  }
  admat
}

#' @export
mutate.column <- function(admat, k) {
  ## possible positions (0's)
  possiblePos <- which(!is.na(admat[, k]) & admat[, k] != 1)
  ## current position with 1
  ind.1 <- which(admat[, k] == 1)
  ## select new position
  if (length(possiblePos) == 1) {
    new.1 <- possiblePos
  } else {
    new.1 <- sample(possiblePos, size=1)
  }
  
  admat[ind.1, k] <- 0
  admat[new.1, k] <- 1
  admat
}

#' @export
mutate.admat.3 <- function(admat, ncol.to.mutate, mcf_matrix) {
  
  mutate.prob.tb <- get.cluster.mutate.prob(mcf_matrix)
  
  new.admat <- mutate.n.columns.clusterprob(admat, ncol.to.mutate, mutate.prob.tb)
  
  # make sure admat is fully connected
  while (!is.fully.connected(new.admat)) {
    new.admat <- mutate.n.columns.clusterprob(admat, ncol.to.mutate, mutate.prob.tb)
  }
  
  new.admat
}

#' @export
mutate.n.columns.clusterprob <- function(admat, ncol.to.mutate, mutate.prob.tb) {
  K <- ncol(admat)
  rand.ks <- sample(seq_len(K), size=ncol.to.mutate, prob = mutate.prob.tb$cluster_prob)
  for(k in rand.ks){
    admat <- mutate.column(admat, k)
  }
  admat
}

get.cluster.mutate.prob <- function(mcf_matrix) {
  sample.presence <- get.sample.presence(mcf_matrix)
  tiers <- sort(unique(sample.presence$num_samples))
  tier.prob <- get.tier.probs(tiers)
  mutate.prob.tb <- sample.presence %>%
    mutate(tier_prob = tier.prob$tier_prob[match(sample.presence$num_samples, tier.prob$num_samples)])
  mutate.prob.tb <- mutate.prob.tb %>%
    group_by(num_samples) %>%
    mutate(cluster_prob = tier_prob / n())
  mutate.prob.tb
}

get.sample.presence <- function(mcf_matrix) {
  mcf <- round(mcf_matrix, 2)
  sample.presence <- tibble(cluster = paste0("cluster", 1:nrow(mcf_matrix)),
                            num_samples = rowSums(mcf > 0))
  sample.presence
}

get.tier.probs <- function(tiers) {
  numTiers <- length(tiers)
  tibble(num_samples = tiers, 
         tier_prob = round(rev(10^(numTiers))/sum(10^(numTiers)), 5))
  #tier_prob = rev(1:numTiers / sum(1:numTiers)))
}

is.fully.connected <- function(admat) {
  # checks if admat is fully connected
  numClusters <- nrow(admat)
  nodesInMainTree <- bfs(admat)
  numNodesInMainTree <- length(nodesInMainTree)
  numClusters == numNodesInMainTree
}

bfs <- function(admat) {
  # starting at root
  children <- names(which(admat[1, ] == 1))
  nodes <- c("root", children)
  
  while(length(children) > 0) {
    c <- children[1]
    children <- children[-1]
    temp.children <- names(which(admat[c, ] == 1))
    children <- c(children, temp.children)
    nodes <- c(nodes, temp.children)
  }
  nodes
}

is.bidirectional <- function(admat) {
  numClusters <- ncol(admat)
  for (i in seq_len(numClusters)) {
    for (j in seq_len(numClusters)) {
      s <- sum(admat[i+1, j], admat[j+1, i], na.rm = TRUE)
      if (s == 2) return(TRUE)
    }
  }
  
  FALSE
}

fix.bidirectional <- function(admat) {
  numClusters <- ncol(admat)
  for (i in seq_len(numClusters)) {
    for (j in seq_len(numClusters)) {
      s <- sum(admat[i+1, j], admat[j+1, i], na.rm = TRUE)
      if (s == 2) {
        # pick one edge to change
        to <- sample(c(i,j), size = 1)
        
        new.admat <- admat
        ## possible positions (0's)
        possiblePos <- which(!is.na(admat[, to]) & admat[, to] != 1)
        ## current position with 1
        ind.1 <- which(admat[, to] == 1)
        ## select new position
        if (length(possiblePos) == 1) {
          new.1 <- possiblePos
        } else {
          new.1 <- sample(possiblePos, size=1)
        }
        
        new.admat[ind.1, to] <- 0
        new.admat[new.1, to] <- 1
        return(new.admat)
      }
    }
  }
}


plotDAG <- function(admat){
  admat.untouched <- admat
  admat <- cbind(0, admat) ## add column for root
  dimnames(admat)[[2]][1] <- "root"
  dimnames(admat) <- lapply(dimnames(admat), function(x) gsub("cluster", "", x))
  
  admat[is.na(admat)] <- 0
  
  net <- network::network(admat, directed=TRUE)
  GGally::ggnet2(net, label=TRUE, arrow.size=12,
                 arrow.gap=0.025, mode = get.DAG.coords.2(admat.untouched))
}

get.DAG.coords <- function(admat) {
  dat <- data.frame(label = rownames(admat), 
                    x = 0,
                    y = 0)
  
  # fix root position
  dat[dat$label=="root", ]$x <- 0.5
  dat[dat$label=="root", ]$y <- 1
  
  
  lvls <- getNodesInLevels(admat)
  #lvls <- lapply(lvls, function(x) gsub("cluster", "", x))
  
  yvals <- seq(1, 0, by = -1/length(lvls))[-1]
  for (i in 1:length(lvls)) {
    nodes <- lvls[[i]]
    dat[match(nodes, dat$label), ]$y <- yvals[i]
    
    xvals <- seq(0, 1, by = 1/(length(nodes) + 1))[-c(1, length(nodes) + 2)]
    dat[match(nodes, dat$label), ]$x <- xvals
  }
  cbind(dat$x, dat$y)
}


get.DAG.coords.2 <- function(admat) {
  nodeInfo <- getNodeInfo(admat)
  nodeInfo$x <- 0
  nodeInfo$y <- 0
  
  # fix root position
  nodeInfo[nodeInfo$node=="root", ]$x <- 0.5
  nodeInfo[nodeInfo$node=="root", ]$y <- 1
  
  yvals <- seq(1, 0, by = -1/max(nodeInfo$level))[-1]
  
  
  for (i in 1:max(nodeInfo$level)-1) {
    
    parents <- nodeInfo[nodeInfo$level == i, ]$node 
    for (parent in parents) {
      if(nodeInfo[nodeInfo$node == parent, ]$numKids == 0) next
      kids <- nodeInfo[which(nodeInfo$parent == parent), ]$node
      
      # set y vals
      nodeInfo[match(kids, nodeInfo$node), ]$y <- yvals[i+1]
      
      # set x vals
      if (parent == "root") {
        xvals <- seq(0, 1, by = 1/(length(kids) + 1))[-c(1, length(kids) + 2)]
        nodeInfo[match(kids, nodeInfo$node), ]$x <- xvals
      } else {
        p.x <- nodeInfo[which(nodeInfo$node == parent), ]$x
        if (length(kids) == 1) {
          nodeInfo[match(kids, nodeInfo$node), ]$x <- p.x
        } else {
          r <- 0.025*(max(nodeInfo$level)-i)* length(kids)
          xvals <- seq(p.x - (r/2), p.x + (r/2), length.out = length(kids))
          nodeInfo[match(kids, nodeInfo$node), ]$x <- xvals
        }
        
      }
    }
  }
  cbind(nodeInfo$x, nodeInfo$y)
}

getNodeInfo <- function(admat) {
  nodeInfo <- getLevels(admat)
  nodeInfo$parent <- NA
  for (r in 2:nrow(nodeInfo)) {
    nodeInfo$parent[r] <- names(which(admat[, r-1] == 1))
    
  }
  nodeInfo$numKids <- rowSums(admat, na.rm = T)
  nodeInfo
}

getLevels <- function(admat) {
  nodeNames <- rownames(admat)
  numNodes <- nrow(admat)
  
  lvl <- data.frame(node = nodeNames,
                    level = 0, 
                    stringsAsFactors = F)
  
  currParents <- "root"
  currKids <- unname(unlist(sapply(currParents, function(x) names(which(admat[x, ] == 1)))))
  currlvl <- 1
  
  while (length(currKids) > 0) {
    lvl[match(currKids, lvl$node), ]$level <- currlvl
    currParents <- currKids
    currKids <- unname(unlist(sapply(currParents, function(x) names(which(admat[x, ] == 1)))))
    
    currlvl <- currlvl + 1
  }
  lvl
}

getNodesInLevels <- function(admat) {
  nodeNames <- rownames(admat)
  numNodes <- nrow(admat)
  
  lvls <- list()
  
  currParents <- "root"
  currKids <- unname(unlist(lapply(currParents, function(x) names(which(admat[x, ] == 1)))))
  currlvl <- 1
  
  while (length(currKids) > 0) {
    lvls[[currlvl]] <- currKids
    currParents <- currKids
    currKids <- unname(unlist(lapply(currParents, function(x) names(which(admat[x, ] == 1)))))
    currlvl <- currlvl + 1
  }
  lvls
}



numericRepresentation <- function(x){
  x[is.na(x)] <- 0
  x <- as.numeric(x)
  paste(x, collapse="")
}

plainAdmat <- function(admat) {
  plain <- cbind(root = rep(0, nrow(admat)), admat)
  plain[is.na(plain)] <- 0
  plain
}

calc.tree.prop.true <- function(admat, true.admat) {
  # inputs:
  #   - admat = an adjacency matrix
  #   - true.admat = the true adjacency matrix 
  # output: proportion of edges that are correct [0,1]
  true.edges <- which(true.admat == 1)
  admat.edges <- which(admat == 1)
  sum(admat.edges %in% true.edges) / length(true.edges)
}

plotEnsembleDAG <- function(post.admat, filter1 = TRUE, filter1.threshold = 0.1) {
  # filter1 filters columns for edges with posterior prob > (max(column) - filter1.threshold)
  
  admat <- cbind(0, post.admat) ## add column for root
  dimnames(admat)[[2]][1] <- "root"
  dimnames(admat) <- lapply(dimnames(admat), function(x) gsub("cluster", "", x))
  admat <- as.matrix(admat)
  
  ad <- admat
  # filter edges
  if (filter1) {
    #thresh <- apply(admat, 2, max) - filter1.threshold
    ad <- apply(admat, 2, function(x) ifelse(x > (max(x)-filter1.threshold), x, 0))
  }
  
  ig <- graph_from_adjacency_matrix(ad, mode = "directed", weighted = TRUE,
                                    diag = FALSE, add.row = TRUE) 
  
  E(ig)$lty <- ifelse(E(ig)$weight < 0.25, 2, 1)
  
  # make edge black if only 1 edge to vertex
  e <- ends(ig, E(ig))
  numTo <- table(e[,2])
  edgeColors <- sapply(e[,2], function(x) ifelse(x %in% names(which(numTo==1)), "black", "darkgrey"))
  E(ig)$color <- edgeColors
  
  V(ig)$label.cex <- 0.5
  
  plot.igraph(ig, layout = layout_as_tree(ig),
              vertex.color = "white", vertex.label.family = "Helvetica",
              edge.arrow.size = 0.2, edge.arrow.width = 2,
              edge.width = E(ig)$weight*3)
}

mcfMatrix <- function(mcf_stats, parameter="mean"){
  K <- numberClusters(mcf_stats)
  S <- numberSamples(mcf_stats)
  if(parameter=="mean")
    MCF <- matrix(mcf_stats$mean, K, S, byrow=TRUE)
  if(parameter=="sd")
    MCF <- matrix(mcf_stats$sd, K, S, byrow=TRUE)
  MCF
}
