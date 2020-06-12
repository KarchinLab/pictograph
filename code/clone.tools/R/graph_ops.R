## refactored base.admat
constrainedEdges <- function(w, zero.thresh=0.01) {
    ##
    ## Rules:
    ##  - cluster (node) cannot connect to itself
    ##  - a cluster with near-zero MCF cannot have children
    ##  - a cluster present in X multiple samples cannot connect to a cluster present in Y samples
    ##    if X < Y
    ##         - X < Y implies ...
    ##
    ##cluster.sample.presence <- apply(w, 1, function(x) which( x> zero.thresh))
    wmat <- spread(w, sample_id, mean) %>%
        ungroup() %>%
        select(-cluster_id) %>%
        as.matrix()
    cluster.sample.presence <- apply(wmat, 1, function(x) which(x>zero.thresh))
    K <- length(unique(w$cluster_id))
    S <- length(unique(w$sample_id))
    admat <- matrix(0, K, K)
    for(i in 1:K){
        for(j in 1:K){
            from.samples <- cluster.sample.presence[[i]]
            to.samples <- cluster.sample.presence[[j]]
            if (setequal(from.samples, to.samples)) next()
            if(length(from.samples) < length(to.samples)) {
                admat[i, j] <- NA
                next()
            }
            if (all(to.samples %in% from.samples)) next()
            admat[i, j] <- NA 
      }            
    }
    diag(admat) <- NA
    am2 <- rbind(0, admat)
    dimnames(am2) <- list(c("root", 1:K), 1:K)
    am2.long <- as_tibble(am2) %>%
        mutate(parent=rownames(am2)) %>%
        pivot_longer(-parent,
                     names_to="child",
                     values_to="connected") %>%
        filter(parent != child) %>%
        unite("edge", c("parent", "child"), sep="->",
              remove=FALSE) %>%
        mutate(parent=factor(parent, levels=unique(parent)))    
    am2.long
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

randAdmat <- function(am.long, max.num.root.children) {
  # input: blank am.long
  # output: random graph
  
  am.long$parent <- as.character(am.long$parent)
  
  all.nodes <- getAllNodes(am.long)
  not.root.nodes <- all.nodes[all.nodes != "root"]
  
  # choose node to connect to root
  temp.node <- sample(not.root.nodes, 1)
  am.long[am.long$edge == getEdgeName("root", temp.node), ]$connected <- 1
  not.root.nodes <- not.root.nodes[not.root.nodes != temp.node]
  from.nodes <- c("root", temp.node)
  
  while(length(not.root.nodes) > 0) {
    
    # remove "root" from possible parents if max.num.root.children quota satisfied
    if(numNodesConnectedToRoot(am.long) > max.num.root.children) {
      from.nodes.pool <- from.nodes
    } else {
      from.nodes.pool <- from.nodes[-1]
    }
      
    temp.to <- sample(not.root.nodes, 1)
    temp.from <- sample(from.nodes.pool, 1)
    am.long[am.long$edge == getEdgeName(temp.from, temp.to), ]$connected <- 1
    from.nodes <- c(from.nodes, temp.to)
    not.root.nodes <- not.root.nodes[not.root.nodes != temp.to]
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
  nodes
}


addEdge <- function(am, new_edge){
    ## Add the edge
    new_edge$connected <- 1
    disconnect_edge <- filter(am, child == new_edge$child)  %>%
        filter(connected==1) %>%
        mutate(connected=0)
    updated_edges <- bind_rows(new_edge, disconnect_edge)
    am2 <- filter(am, ! edge %in% updated_edges$edge ) %>%
        bind_rows(updated_edges) %>%
        arrange(parent) %>%
        updateGraphElements()
    am2
}

initializeGraph <- function(mcf, zero.thresh=0.01){
    clusters <- seq_len(nrow(mcf))
    nsamp <- ncol(mcf)
    samples <- seq_len(nsamp)
    nclust <- length(clusters)
    mcf.long <- tibble(cluster_id=as.character(rep(clusters, nsamp)),
                       sample_id=as.character(rep(samples, each=nclust)),
                       mean=as.numeric(mcf))    
    am.long <- constrainedEdges(mcf.long, zero.thresh=zero.thresh)
    am.long2 <- randAdmat(am.long)
    am.long2
}


toWide <- function(am.long){
    am.long$child <- as.numeric(am.long$child)
    am.long %>% select(parent, child, connected) %>%
        spread(child, connected) %>%
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
  am.long
}

sampleNewEdge <- function(a){
    condition1 <- a %>% group_by(child) %>%
        summarize(n=sum(connected==0)) %>%
        mutate(possible=which(n >= 1))
    ## condition 2:
    ##   for each child and for each zero of that child
    ##     determine whether changing the 1 to a zero and the zero to a 1 would be valid
    possible_moves <- filter(a, connected==0, child %in% condition1$child)
    is_valid <- rep(NA, nrow(possible_moves))
    for(i in seq_len(nrow(possible_moves))){
        astar <- addEdge(a, possible_moves[i, ])
        is_valid[i] <- validGraph(astar)
    }
    move_set <- possible_moves[is_valid, ]
    ix <- sample(seq_len(nrow(move_set)), 1)
    astar <- addEdge(a, move_set[ix, ])
    astar
}

initEmptyAdmatFromK <- function(K) {
  admat <- matrix(0, K, K)
  diag(admat) <- NA
  am2 <- rbind(0, admat)
  dimnames(am2) <- list(c("root", 1:K), 1:K)
  am2
}

generateRandomGraphFromK <- function(K, max.num.root.children) {
  # input: number of mutation clusters, K
  # output: mutation tree; adjacency matrix
  am.long <- toLong(initEmptyAdmatFromK(K))
  rand.am.long <- randAdmat(am.long, max.num.root.children)
  if (!validGraph(rand.am.long)) warning("graph is not valid")
  rand.am.long
}

plotGraph <- function(am.long){
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
