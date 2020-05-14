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
    dimnames(am2) <- list(c("root", 1:4), 1:4)
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

reversedConnection <- function(am){
    connections <- setNames(am$connected, am$edge)
    reversed_connections <- connections[am$reverse_edge] %>%
        "["(!is.na(.))
    reversed <- setNames(rep(0, nrow(am)), am$reverse_edge)
    reversed[names(reversed_connections)] <- reversed_connections
    reversed
}

isBidirectional <- function(am){
    am %>%
        mutate(bi_directional=(reverse_edge %in% edge) &
                   connected==1 &
                   reversed_connected == 1) %>%
        pull(bi_directional)
}

updateGraphElements <- function(am){
    am %>%
        mutate(reversed_connected=reversedConnection(.)) %>%
        mutate(bi_directional=isBidirectional(.)) %>%
        mutate(root_connected=isRootConnected(.))
}


reversedEdges <- function(am){
    am2 <- am %>%
        filter(!is.na(connected)) %>%
        unite("reverse_edge", c("child", "parent"), sep="->",
              remove=FALSE)
    am2
}

randAdmat <- function(am.long){
    new_edges <- am.long %>%
        filter(connected==0) %>%
        group_by(child) %>%
        sample_n(1) %>%
        mutate(connected=1) %>%
        ungroup()
    am3.long <- filter(am.long, !edge %in% new_edges$edge) %>%
        bind_rows(new_edges) %>%
        arrange(parent, child) 
    am <- reversedEdges(am3.long) %>%
        mutate(reversed_connection=reversedConnection(.),
               bi_directional=NA,
               root_connected=NA)    
    am  <- updateGraphElements(am)
    am
}

isParentConnected <- function(am){
    am %>%
        group_by(parent) %>%
        summarize(n=sum(connected)) %>%
        pull(n) > 0
}

isRootConnected <- function(am) isParentConnected(am)[1]

isDirected <- function(am) !any(am$bi_directional)

validGraph <- function(am){
    isDirected(am) &&
        isRootConnected(am)
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
    am.long %>% select(parent, child, connected) %>%
        spread(child, connected) %>%
        select(-parent) %>%
        as.matrix()
}


sampleNewEdge <- function(a){
    condition1 <- a %>% group_by(child) %>%
        summarize(n=sum(connected==0)) %>%
        mutate(possible=which(n >= 1))
    ## condition 2:
    ##   for each child and for each zero of that child
    ##     determine whether changing the 1 to a zero and the zero to a 1 would be valid
    possible_moves <- filter(a, connected==0, child %in% condition1$child)
    is_valid <- rep(NA, nrow(atmp))
    for(i in seq_len(nrow(atmp))){
        astar <- addEdge(a, atmp[i, ])
        is_valid[i] <- validGraph(astar)
    }
    move_set <- possible_moves[is_valid, ]
    ix <- sample(seq_len(nrow(move_set)), 1)
    astar <- addEdge(a, move_set[ix, ])
    astar
}
