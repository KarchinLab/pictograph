simulateData <- function(I, K, S){
    pi <- rep(0.1, 10)
    ##z <- sample(1:K, size = I, replace = T, prob = pi)
    ## True cancer cell fraction
    z <- rep(1:10, each=10)
    w <- matrix(c(0.98, 0.99, 0.97, 
                  0.98, 0.90, 0.82,
                  0.55, 0.00, 0.80, 
                  0.20, 0.00, 0.50,
                  0.30, 0.00, 0.30, 
                  0.43, 0.90, 0.00,
                  0.30, 0.70, 0.00,
                  0.20, 0.00, 0.00,
                  0.00, 0.00, 0.30,
                  0.00, 0.50, 0.00),
                byrow=T,
                nrow=K, ncol=S)

    colnames(w) <-  paste0("sample", 1:S)
    w

    tcn <- matrix(2, nrow=I, ncol=S)
    m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
                nrow=I, ncol=S)
    W <- w[z, ]
    calcTheta <- function(m, tcn, w) {
      (m * w) / (tcn * w + 2*(1-w))
    }
    theta <- calcTheta(m, tcn, W)
    ##
    ## Simulate altered reads
    ##
    n <- replicate(S, rpois(I, 100))
    y <- matrix(NA, nrow=I, ncol=S)
    for (i in 1:I) {
      for (s in 1:S) {
        y[i, s] <- rbinom(1, n[i, s], theta[i,s])
      }
    }
    
    p <- y/n
    colnames(w) <- colnames(p) <- paste0("sample", 1:S)
    colnames(m) <- colnames(p)
    test.data <- list("I" = I, "S" = S, "K" = K, 
                      "y" = y, "n" = n,
                      "m" = m, "tcn" = tcn,
                      z=z,
                      w=w,
                      theta=theta,
                      p=p)
    test.data
}

runMCMC <- function(data, K, jags.file, inits, params,
                    n.iter, thin, n.chains=1,
                    n.adapt=1000) {
    data$K <- K
    data <- data[c("I", "S", "y", "n",
                   "m", "tcn", "K")]
    jags.m <- jags.model(jags.file,
                         data,
                         n.chains = n.chains,
                         inits = inits,
                         n.adapt = n.adapt)
    samps <- coda.samples(jags.m, params, n.iter=n.iter, thin=thin)
    samps
}

getParamChain <- function(samps, param) {
  chains <- do.call(rbind, samps)
  chain <- chains[, grep(param, colnames(chains))]
}

reshapeW <- function(w, S, K) {
  w.mat <- matrix(w, nrow = K)
  colnames(w.mat) <- paste0("sample", 1:S)
  w.mat
}

calcLogLik <- function(z.iter, w.iter, data) {
  W <- w.iter[z.iter, ]
  theta <- calcTheta(data$m, data$tcn, W)
  sum(dbinom(data$y, data$n, theta, log=T))
}

calcChainLogLik <- function(samps, data, K) {
  z.chain <- getParamChain(samps, "z")
  w.chain <- getParamChain(samps, "w")
  lik <- c()
  for(iter in 1:nrow(z.chain)) {
    z.iter <- z.chain[iter,]
    w.iter <- reshapeW(w.chain[iter,], data$S, K)
    lik <- c(lik, calcLogLik(z.iter, w.iter, data))
  }
  mean(lik)
}

calcBIC <- function(n, k, ll) log(n)*k - 2*ll

relabelZ <- function(z.chain){
  mcmc_z <- z.chain %>%
      group_by(Parameter, value) %>%
      summarize(n=n(),
                maxiter=max(Iteration)) %>%
      mutate(probability=n/maxiter)
  map_z  <- mcmc_z %>%
      group_by(Parameter) %>%
      summarize(map=value[probability==max(probability)]) %>%
      mutate(original_label=rep(1:10, each=10)) %>%
      ungroup() %>%
      group_by(original_label) %>%
      summarize(value=names(sort(table(map), decreasing=TRUE))[1],
                value=as.numeric(value)) %>%
      ungroup()
  sim_cluster <- mcmc_z %>%
      group_by(Parameter) %>%
      summarize(n=n()) %>%
      mutate(sim_cluster=rep(1:10, each=10)) %>%
      select(-n)
  mcmc_z <- z.chain %>%
      left_join(map_z, by="value") %>%
      left_join(sim_cluster, by="Parameter") %>%
      mutate(chain_value=value,
             value=original_label) %>%
      select(-original_label) %>%
      group_by(Parameter, value) %>%
      summarize(n=n(),
                maxiter=max(Iteration),
                sim_cluster=unique(sim_cluster)) %>%
      mutate(probability=n/maxiter)
      ##filter(value==sim_cluster)
}

orderW <- function(w, map_z){
    mcmc_cluster_numbering  <- matrix(map_z$value, nrow = 10)
    get_mode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
    }
    mcmc_cluster_numbering <- apply(mcmc_cluster_numbering, 2, get_mode)
    true_to_mcmc_w_ordering <- match(1:K, mcmc_cluster_numbering)
    w_ordered <- w[true_to_mcmc_w_ordering, ]
    mcmc_vals <- w.chain %>%
        group_by(Parameter) %>%
        summarize(mean=mean(value),
                  q1=quantile(value, 0.025),
                  q3=quantile(value, 0.975)) %>%
        mutate(truth=as.numeric(t(w_ordered)))    
    mcmc_vals
}

rand.admat <- function(admat) {
  for(col in 1:ncol(admat)) {
    ind.0 <- which(admat[,col] == 0) # possible positions (0's)
    rand.ind <- sample(ind.0, size=1)
    admat[rand.ind,col] <- 1
  }
  
  while (sum(admat[1, ]) == 0) {
    admat <- mutate.admat(admat)
  }
  
  admat
}

numberClones <- function(mcf_stats){
    params <- as.character(mcf_stats$Parameter)
    K <- strsplit(params, ",") %>%
        sapply("[", 1) %>%
        str_replace("w\\[", "") %>%
        as.numeric() %>%
        max()
}

numberSamples <- function(mcf_stats){
    params <- as.character(mcf_stats$Parameter)    
    nSamples <- strsplit(params, ",") %>%
        sapply("[", 2) %>%
        str_replace("\\]", "") %>%
        as.numeric() %>%
        max()
}

mcfMatrix <- function(mcf_stats, parameter="mean"){
    K <- numberClones(mcf_stats)
    S <- numberSamples(mcf_stats)
    if(parameter=="mean")
        MCF <- matrix(mcf_stats$mean, K, S, byrow=TRUE)
    if(parameter=="sd")
        MCF <- matrix(mcf_stats$sd, K, S, byrow=TRUE)
    MCF
}

##base.admat <- function(MCF, zero.thresh=0.01) {
initializeAdjacencyMatrix <- function(mcf_stats, zero.thresh=0.01) {
    ##
    ## for each row (clone) of the MCF matrix,
    ## list indices of samples for which the clone is present
    ##
    MCF <- mcfMatrix(mcf_stats)
    K <- nrow(MCF)
    S <- ncol(MCF)
    cluster.sample.presence <- apply(MCF, 1, function(x)
        which(x > zero.thresh))
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
    dimnames(admat) <- list(c("root", paste0("clone", seq_len(K))),
                            paste0("clone", seq_len(K)))
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
    ##for(rand.k in rand.ks) {
    for(k in rand.ks){
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
        new.admat[ind.1, k] <- 0
        new.admat[new.1, k] <- 1
    }
    ## Recursion -- why?
    ## - requires that the root be connected
    ##    while (sum(new.admat[1, ]) == 0) {
    ##        new.admat <- mutate.admat(admat)
    ##    }
    new.admat
}

decide.ht <- function(pval, alpha=0.05) {
  # 1 signals rejection event for null of i -> j
  if (pval <= alpha) return(1)
  else return(0)
}

##create.cpov <- function(mcmc_w, mcmc_w_sd, alpha=0.05, zero.thresh=0.01) {
##
## What is the purpose of this function?
##
create.cpov <- function(mcf_stats, alpha=0.05, zero.thresh=0.01) {
    ##cpov <- initializeAdjacencyMatrix(mcmc_w, zero.thresh = 0.01)
    cpov <- initializeAdjacencyMatrix(mcf_stats, zero.thresh = zero.thresh)
    MCF <- mcfMatrix(mcf_stats)
    sds <- mcfMatrix(mcf_stats, parameter="sd")
    ##S <- ncol(mcmc_w) # number of samples
    S <- numberSamples(mcf_stats)
    ##cpov <- cpov[-1, ]
    ## root can go to anyone -- all 0's (default base admat value)
    for (r in 2:nrow(cpov)) {
        for (c in 1:ncol(cpov)) {
            if (is.na(cpov[r,c])) next # skip restricted position
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

calc.topology.cost <- function(admat, cpov) {
    TC <- 0
    edges <- which(admat == 1, arr.ind=TRUE)
    N <- nrow(edges)
    for (i in seq_len(N)) {
        TC <- TC + cpov[edges[i,1], edges[i,2]]
    }
  TC
}

calc.mass.cost <- function(admat, mcmc_w) {
    numChildren <- rowSums(admat, na.rm = T)
    nodes <- which(numChildren > 0, arr.ind = T) # not leaves
    mc.node <- rep(0, length(nodes))

    for (i in 1:length(nodes)) {
      node <- nodes[i]

      # root node: MCF = 1
      parent.w <- rep(1, ncol(mcmc_w))
      # not root node: look up MCF in mcmc_w
      if (node != 1) { 
        parent.w <- mcmc_w[node-1,]
      }

      kids <- which(admat[node,] == 1, arr.ind = T)
      if (numChildren[node] > 1) {
        children.w <- colSums(mcmc_w[kids,])
      } else {
        children.w <- mcmc_w[kids,]
      }

      mc.s <- ifelse(parent.w >= children.w, 0, children.w - parent.w)
      mc.node[i] <- sqrt(sum(mc.s^2))
    }
    sum(mc.node)
}

calc.tree.fitness <- function(admat, cpov, mcf_stats, scaling.coeff=5) {
    TC <- calc.topology.cost(admat, cpov)
    MC <- calc.mass.cost(admat, mcfMatrix(mcf_stats))
    Z <- TC + MC
    fitness <- exp(-scaling.coeff * Z)
    fitness
}


plotDAG <- function(admat){
    dimnames(admat) <- NULL
    admat <- cbind(0, admat) ## add column for root
    admat[is.na(admat)] <- 0
    ## directed graph is not possible without removing root
    ##admat2 <- admat2[-1, ]
    net <- network(admat, directed=TRUE)
    ggnet2(net, label=TRUE, arrow.size=12,
           arrow.gap=0.025)
}
