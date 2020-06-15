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
        statistic <- statistic + (d / d_sd)^2 * I
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

calc.topology.cost <- function(admat, cpov) {
  TC <- 0
  edges <- which(admat == 1, arr.ind=TRUE)
  N <- nrow(edges)
  for (i in seq_len(N)) {
    TC <- TC + cpov[edges[i,1], edges[i,2]]
  }
  TC
}

calc.mass.cost <- function(admat, mcf_matrix) {
  numChildren <- rowSums(admat, na.rm = T)
  nodes <- which(numChildren > 0, arr.ind = T) # not leaves
  mc.node <- rep(0, length(nodes))
  
  for (i in 1:length(nodes)) {
    node <- nodes[i]
    
    # root node: MCF = 1
    parent.w <- rep(1, ncol(mcf_matrix))
    # not root node: look up MCF in mcf_matrix
    if (node != 1) { 
      parent.w <- mcf_matrix[node-1,]
    }
    
    kids <- which(admat[node,] == 1, arr.ind = T)
    if (numChildren[node] > 1) {
      children.w <- colSums(mcf_matrix[kids,])
    } else {
      children.w <- mcf_matrix[kids,]
    }
    
    mc.s <- ifelse(parent.w >= children.w, 0, children.w - parent.w)
    mc.node[i] <- sqrt(sum(mc.s^2))
  }
  sum(mc.node)
}

calc.tree.fitness <- function(admat, cpov, mcf_matrix, weight.mass = 1, weight.topology = 1, scaling.coeff=5) {
  TC <- calc.topology.cost(admat, cpov)
  MC <- calc.mass.cost(admat, mcf_matrix)
  Z <- weight.topology * TC + weight.mass * MC
  fitness <- exp(-scaling.coeff * Z)
  fitness
}