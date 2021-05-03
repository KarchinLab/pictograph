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

mutate.admat.2 <- function(admat, ncol.to.mutate) {
  
  new.admat <- mutate.n.columns(admat, ncol.to.mutate)
  
  # make sure admat is fully connected
  while (!is.fully.connected(new.admat)) {
    new.admat <- mutate.n.columns(admat, ncol.to.mutate)
  }
  
  new.admat
}

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

mutate.admat.3 <- function(admat, ncol.to.mutate, mcf_matrix) {

  mutate.prob.tb <- get.cluster.mutate.prob(mcf_matrix)
  
  new.admat <- mutate.n.columns.clusterprob(admat, ncol.to.mutate, mutate.prob.tb)
  
  # make sure admat is fully connected
  while (!is.fully.connected(new.admat)) {
    new.admat <- mutate.n.columns.clusterprob(admat, ncol.to.mutate, mutate.prob.tb)
  }
  
  new.admat
}

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