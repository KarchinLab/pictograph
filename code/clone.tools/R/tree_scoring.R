summarizeWChain <- function(w.chain) {
  # output: mcf_stats 
  mcf_stats <- w.chain %>%
    group_by(Parameter) %>%
    summarize(sd=sd(value),
              mean=mean(value))
  return(mcf_stats)
}

create.cpov <- function(mcf_stats, alpha=0.05, zero.thresh=0.01, mcf_matrix = NULL, restriction.val = 1) {
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

calcMassCost <- function(am, mcf_matrix, am_format="long") {
  num_samples <- ncol(mcf_matrix)
  
  if (am_format == "long") {
    edges <- getEdges(am)
    
    parent_nodes <- unique(edges$parent)
    mass_cost <- rep(0, length(parent_nodes)) # mass cost of each parent node
    
    for (i in seq_len(length(parent_nodes))) {
      parent_node <- parent_nodes[i]
      
      # root CCF is 1
      if (parent_node == "root") {
        parent_w <- rep(1, num_samples)
      } else {
        parent_w <- mcf_matrix[as.numeric(parent_node), ]
      }
      
      kids <- getChildren(am, parent_node)
      if (length(kids) > 1) {
        children_w <- colSums(mcf_matrix[as.numeric(kids), ])
      } else {
        children_w <- mcf_matrix[as.numeric(kids), ]
      }
      
      mc_s <- ifelse(parent_w >= children_w, 0, children_w - parent_w)
      mass_cost[i] <- sqrt(sum(mc_s^2))
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
  am_wide <- initEmptyAdmatFromK(length(edges$child))
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

calcTreeFitness <- function(admat, cpov, mcf_matrix, am_format = "long", weight_mass = 1, weight_topology = 1, scaling_coeff=5) {
  # if only edges are given, change into long format
  if (am_format == "edges") {
    admat <- edgesToAmLong(admat)
    am_format <- "long"
  }
  
  TC <- calcTopologyCost(admat, cpov, am_format)
  MC <- calcMassCost(admat, mcf_matrix, am_format)
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