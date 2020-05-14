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
    admat <- rbind(0, admat)
    admat
}

isParentConnected <- function(am){
    am %>%
        group_by(parent) %>%
        summarize(n=sum(connected)) %>%
        pull(n) > 0
}

isRootConnected <- function(am) isParentConnected(am)[1]

addEdge <- function(am, new_edge){
    ## Add the edge
    new_edge$connected <- 1
    disconnect_edge <- filter(am, child == new_edge$child)  %>%
        filter(connected==1) %>%
        mutate(connected=0)
    updated_edges <- bind_rows(new_edge, disconnect_edge)
    am2 <- filter(am, ! edge %in% updated_edges$edge ) %>%
        bind_rows(updated_edges) %>%
        arrange(parent)
    am2
}
