clusterSep <- function(input_data,
                       n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                       inits = list(".RNG.name" = "base::Wichmann-Hill",
                                    ".RNG.seed" = 123),
                       max_K = 5,
                       model_type = "spike_and_slab",
                       beta.prior = FALSE) {
  
  # 2a. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
  params = c("z", "w", "ystar")
  
  # Max number of clusters cannot be more than number of mutations
  temp_max_K <- min(max_K, input_data$I)
  
  all_set_results <- runMutSetMCMC(input_data, 
                                   n.iter = n.iter, n.burn = n.burn, thin = thin, 
                                   mc.cores = mc.cores,
                                   inits = inits,
                                   temp_max_K = temp_max_K,
                                   model_type = model_type,
                                   params = params,
                                   beta.prior = beta.prior)
  
  return(all_set_results)
}


runMCMCForABox <- function(box, 
                           n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                           inits = list(".RNG.name" = "base::Wichmann-Hill",
                                        ".RNG.seed" = 123),
                           model = model_type,
                           params = c("z", "w", "ystar"),
                           max_K = 5, beta.prior = FALSE) {
  # returns samps_list 
  box_input_data <- getBoxInputData(box)
  
  #extdir <- system.file("extdata", package="pictograph")
  extdir <- "C:/Users/Leste/OneDrive - Johns Hopkins/Documents/GitHub/pictograph/inst/extdata"
  jags.file.K1 <- file.path(extdir, "single_sample_K1.jags")
  jags.file <- file.path(extdir, "single_sample.jags")
  
  # choose sample in which mutations are present
  
  samps_K1 <- runMCMC(box_input_data, 1, jags.file.K1, 
                      inits, params, n.iter=n.iter, thin=thin, n.burn=n.burn)
  
  
  # Max number of clusters cannot be more than number of mutations
  max_K <- min(max_K, length(box$mutation_indices)) 
  if (max_K > 1) {
    samps_2 <- parallel::mclapply(2:max_K,
                                  function(k) runMCMC(box_input_data, k,
                                                      jags.file, inits, params,
                                                      n.iter=n.iter, thin=thin,
                                                      n.burn=n.burn,
                                                      beta.prior=beta.prior),
                                  mc.cores=mc.cores)
    
    samps_list <- c(list(samps_K1), samps_2)
    names(samps_list) <- paste0("K", 1:max_K)
    return(samps_list)
    
  } else {
    names(samps_K1) <- "K1"
    return(samps_K1)
  }
}


reshapeW <- function(w.chain.iter) {
  w.mat <- w.chain.iter %>%
    mutate(cluster=gsub("w\\[(.*)\\]","\\1",Parameter),
           sample = 1)%>%
    select(cluster,sample, value) 
  S <- 1
  w.mat <- w.mat %>%
    pivot_wider(names_from = sample, 
                values_from = value)
  w.mat$cluster <- NULL
  w.mat <- as.matrix(w.mat)
  colnames(w.mat) <- paste0("sample", 1:S)
  return(w.mat)
}


calcChainBIC <- function(samp_list,k,input.data) {
  n <- input.data$I * input.data$S
  ll <- calcChainLogLik(samp_list[[k]], input.data, k)
  BIC <- calcBIC(n, k, ll)
  return(BIC)
}


runMutSetMCMC <- function(temp_box, 
                          n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                          inits = list(".RNG.name" = "base::Wichmann-Hill",
                                       ".RNG.seed" = 123),
                          temp_max_K = 5,
                          model_type = "spike_and_slab",
                          params = c("z", "w", "ystar"),
                          beta.prior = FALSE) {
  
  # Run MCMC
  if (temp_max_K == 1) {
    # only 1 possible cluster
    temp_samps_list <- runMCMCForABox(temp_box,
                                      n.iter = n.iter, n.burn = n.burn, 
                                      thin = thin, mc.cores = mc.cores,
                                      inits = inits,
                                      params = params,
                                      max_K = temp_max_K)
  } else {
    # assess range of K: [1, temp_max_K]
    temp_samps_list <- runMCMCForABox(temp_box,
                                      n.iter = n.iter, n.burn = n.burn, 
                                      thin = thin, mc.cores = mc.cores,
                                      inits = inits,
                                      params = params,
                                      max_K = temp_max_K,
                                      model = model_type,
                                      beta.prior = beta.prior)
  }
  
  temp <- temp_samps_list 
  # Format chains
  format_chain <- function(temp_samps_list, k){
    
    w <- temp_samps_list[[k]][[1]]%>%
      as.data.frame()%>%
      select(contains("w"))%>%
      mutate(Iteration = seq.int(1,nrow(.),1),
             Chain = rep(1, nrow(.)))%>%
      gather(key = Parameter, value = value, 1:(ncol(.)-2))%>%
      mutate(Parameter= ifelse(Parameter == "w","w[1]",Parameter))
    y <- temp_samps_list[[k]][[1]]%>%
      as.data.frame()%>%
      select(contains("y"))%>%
      mutate(Iteration = seq.int(1,nrow(.),1),
             Chain = rep(1, nrow(.)))%>%
      gather(key = Parameter, value = value, 1:(ncol(.)-2))
    z <- temp_samps_list[[k]][[1]]%>%
      as.data.frame()%>%
      select(contains("z"))%>%
      mutate(Iteration = seq.int(1,nrow(.),1),
             Chain = rep(1, nrow(.)))%>%
      gather(key = Parameter, value = value, 1:(ncol(.)-2))
    
    temp[[k]] <- list(w_chain = w,
                      z_chain = z,
                      ystar_chain = y)
    
  }
  
  samp_list <- lapply(names(temp_samps_list), function(k) format_chain(temp_samps_list,k))
  names(samp_list) <- names(temp)
  
  # Calculate BIC
  K_tested <- seq_len(temp_max_K)
  if (temp_max_K > 1) {
    box_indata <- getBoxInputData(temp_box)
    bic_vec <- unname(unlist(parallel::mclapply(1:length(names(samp_list)), 
                                                function(k) calcChainBIC(samp_list, k, box_indata),
                                                mc.cores = mc.cores)))
    bic_tb <- tibble(K_tested = K_tested,
                     BIC = bic_vec)
    best_chains <- samp_list[[which.min(bic_vec)]]
    res_list <- list(all_chains = samp_list,
                     BIC = bic_tb,
                     best_chains = best_chains,
                     best_K = which.min(bic_vec))
  } else {
    # only 1 variant, so must be 1 cluster and don't need to check BIC
    res_list <- list(all_chains = samps_list,
                     BIC = NA,
                     best_chains = samps_list[[1]],
                     best_K = 1)
  }
  
  return(res_list)
}

estimateCCFs <- function(w_chain) {
  K <- w_chain%>%
    mutate(Cluster = as.numeric(gsub("w\\[(.*)\\]","\\1",Parameter)))%>%
    pull(Cluster)%>%
    max()
  # density plot
  w.dens <- ggplot(w_chain, aes(x = value)) +
    geom_density() +
    facet_wrap(~Parameter, ncol = 1, scales = "free_y") +
    theme_light()
  # find peak for MAP w
  w.dens.p <- ggplot_build(w.dens)$data[[1]]
  w.map <- w.dens.p %>%
    as_tibble() %>%
    group_by(PANEL) %>%
    arrange(desc(y))%>%
    slice(1)%>%
    ## need ungroup()
    ungroup()%>%
    mutate(Parameter = unique(w_chain$Parameter),
           value_rounded = round(x, 2))%>%
    select(PANEL, x,Parameter, value_rounded)
  ########there is no max when everythign is zero############
  #   summarize(value = x[max(y) == y])
  # w.map <- w.map %>%
  #   mutate(Parameter = unique(w_chain$Parameter),
  #          value_rounded = round(value, 2))
  # return w matrix
  w.map.matrix <- matrix(w.map$value_rounded, K, 1, byrow=TRUE)
  return(w.map.matrix)
}

mcfMatrix <- function(mcf_stats, parameter="mean"){
  # changed k and s assignment
  K <- nrow(mcf_stats)
  S <- 1
  if(parameter=="mean")
    MCF <- matrix(mcf_stats$mean, K, S, byrow=TRUE)
  if(parameter=="sd")
    MCF <- matrix(mcf_stats$sd, K, S, byrow=TRUE)
  MCF
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
  #cluster.sample.presence <- lapply(MCF_list, 
  #                                  function(x) which(x > zero.thresh))
  #all.samples <- seq_len(S)
  ## initialize adjacency matrix
  ##admat <- matrix(data=0, nrow=(1+K), ncol=K)
  admat <- matrix(data=0, K, K)
  ## can't go to self
  diag(admat) <- NA
  # removed nested for loop to check box hierarchy
  ## Add root
  ## can go from root to anyone
  admat <- rbind(0, admat)
  dimnames(admat) <- list(c("root", paste0("cluster", seq_len(K))),
                          paste0("cluster", seq_len(K)))
  admat
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
  S <- 1
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


calcTreeScores <- function(w_chain, trees, mc.cores = 1) {
  mcf_stats <- summarizeWChain(w_chain)
  cpov <- create.cpov(mcf_stats)
  w_mat <- estimateCCFs(w_chain)%>%
    cbind(matrix(rep(0, nrow(.)),nrow = nrow(.), ncol = 1))
  schism_scores <- unlist(parallel:::mclapply(trees, 
                                              function(x) calcTreeFitness(x, cpov, w_mat, am_format = "edges"),
                                              mc.cores = mc.cores))
  return(schism_scores)
}