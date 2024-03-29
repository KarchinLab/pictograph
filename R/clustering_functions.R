#' @import rjags
#' @importFrom epiR epi.betabuster
runMCMC <- function(data, K, jags.file, inits, params,
                    n.iter=10000, thin=10, n.chains=1,
                    n.adapt=1000, n.burn=1000,
                    beta.prior=FALSE) {
    if (K > 1) data$K <- K
    jags.m <- jags.model(jags.file,
                         data,
                         n.chains = n.chains,
                         inits = inits,
                         n.adapt = n.adapt)
    if (n.burn > 0) update(jags.m, n.burn)
    samps <- coda.samples(jags.m, params, n.iter=n.iter, thin=thin)
    
    if (beta.prior & K > 1) {
      # use initial MCMC to estimate beta priors for identified clusters
      initial_chains <- formatChains(samps)
      est_ccfs <- initial_chains$w_chain %>%
        estimateCCFs %>%
        as.data.frame %>%
        magrittr::set_colnames(1:ncol(.)) %>%
        tibble %>%
        mutate(cluster = 1:nrow(.)) %>%
        pivot_longer(cols = colnames(.)[colnames(.) != "cluster"], 
                     names_to = "sample", 
                     values_to = "ccf") %>%
        mutate(sample = as.numeric(sample))
      
      
      beta_shapes <- lapply(est_ccfs$ccf, estimateBetaPriors)
      cluster_beta_params <- est_ccfs %>%
        mutate(shape1 = sapply(beta_shapes, function(x) x[1]),
               shape2 = sapply(beta_shapes, function(x) x[2]))
      
      cluster_shape1 <- cluster_beta_params %>%
        select(cluster, shape1, sample) %>%
        pivot_wider(names_from = sample,
                    values_from = shape1) %>%
        select(-c(cluster)) %>%
        as.matrix
      cluster_shape2 <- cluster_beta_params %>%
        select(cluster, shape2, sample) %>%
        pivot_wider(names_from = sample,
                    values_from = shape2) %>%
        select(-c(cluster)) %>%
        as.matrix
      
      # run second MCMC with specified beta priors
      data$cluster_shape1 <- cluster_shape1
      data$cluster_shape2 <- cluster_shape2

      extdir <- system.file("extdata", package="pictograph")
      jags.file.beta <- file.path(extdir, "model-simple-set-beta.jags")
      
      jags.m.beta <- rjags::jags.model(jags.file.beta,
                                  data,
                                  n.chains = n.chains,
                                  inits = inits,
                                  n.adapt = n.adapt)
      if (n.burn > 0) update(jags.m.beta, n.burn)
      samps.beta <- rjags::coda.samples(jags.m.beta, params, n.iter=n.iter, thin=thin)
      return(samps.beta)
    }

    return(samps)
}

#' @importFrom epiR epi.betabuster
estimateBetaPriors <- function(ccf) {
  if (ccf == 0) ccf <- 0.001
  suppressWarnings(
    param_list <- epiR::epi.betabuster(ccf, 
                         conf = 0.5, 
                         greaterthan = TRUE,
                         ccf, 
                         conf.level = 0.8,
                         max.shape1 = 10)
  )
  
  if (length(param_list$shape1) > 1) {
    shape1 <- mean(param_list$shape1)
  } else {
    shape1 <- param_list$shape1
  }
  if (length(param_list$shape2) > 1) {
    shape2 <- mean(param_list$shape2)
  } else {
    shape2 <- param_list$shape2
  }
  return(c(shape1, shape2))
}


getParamChain <- function(samps, param) {
  chains <- do.call(rbind, samps)
  chain <- chains[, grep(param, colnames(chains))]
  chain
}

relabelZ <- function(z.chain, I, K){
  mcmc_z <- z.chain %>%
      group_by(Parameter, value) %>%
      summarize(n=n(),
                maxiter=max(Iteration)) %>%
      mutate(probability=n/maxiter)
  map_z  <- mcmc_z %>%
      group_by(Parameter) %>%
      summarize(map=value[probability==max(probability)]) %>%
      mutate(original_label=rep(1:K, each=I/K)) %>%
      ungroup() %>%
      group_by(original_label) %>%
      summarize(value=names(sort(table(map), decreasing=TRUE))[1],
                value=as.numeric(value)) %>%
      ungroup()
  sim_cluster <- mcmc_z %>%
      group_by(Parameter) %>%
      summarize(n=n()) %>%
      mutate(sim_cluster=rep(1:K, each=I/K)) %>%
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

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

orderW <- function(w, map_z){
    mcmc_cluster_numbering  <- matrix(map_z$value, nrow = 10)
    
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


#' @importFrom stringr str_replace
numberClusters <- function(mcf_stats){
    params <- as.character(mcf_stats$Parameter)
    K <- strsplit(params, ",") %>%
        sapply("[", 1) %>%
        str_replace("w\\[", "") %>%
        as.numeric() %>%
        max()
    K
}

#' @importFrom stringr str_replace
numberSamples <- function(mcf_stats){
    params <- as.character(mcf_stats$Parameter)    
    nSamples <- strsplit(params, ",") %>%
        sapply("[", 2) %>%
        stringr::str_replace("\\]", "") %>%
        as.numeric() %>%
        max()
    nSamples
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

##base.admat <- function(MCF, zero.thresh=0.01) {
chec.mcf_stats.format <- function(mcf_stats) {
  if(all.equal(colnames(mcf_stats), c("Parameter", "sd", "mean"))) TRUE
  else FALSE
}


sample.w <- function(w.chain, K) {
  numIter <- max(w.chain$Iteration)
  randIter <- sample(numIter, size = 1)
  w.sample <- w.chain[w.chain$Iteration == randIter, ]
  matrix(w.sample$value, nrow = K, byrow = T)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability
#' 
#' @export
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
estimateClusterAssignments <- function(z_chain) {
  it <- max(z_chain$Iteration)
  mcmc_z <- z_chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              maxiter=it) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(value=value[probability==max(probability)])
  
  # choose first cluster if equal probability
  map_z_count <- map_z %>% 
    group_by(Parameter) %>%
    summarize(map_count = n()) %>%
    ungroup()
  if (any(map_z_count$map_count > 1)) {
    mut_ind <- which(map_z_count$map_count > 1)
    for (i in mut_ind) {
      dup_var <- as.numeric(gsub("z\\[|]", "", map_z_count$Parameter[i]))
      map_z_dups <- which(gsub("z\\[|]", "", map_z$Parameter) == dup_var)
      dup_ind <- map_z_dups[-1]
      map_z <- map_z[-dup_ind, ]
    }
  }
  return(map_z)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability. 
#' 
#' @export
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
#' @param Mut_ID Vector of mutation IDs, same order as provided as input data (e.g. indata$Mut_ID)
#' @return A tibble listing mutation IDs and their cluster assignments
writeClusterAssignmentsTable <- function(z_chain, Mut_ID = NULL) {
  map_z <- estimateClusterAssignments(z_chain) 
  if (is.null(Mut_ID)) {
    Mut_ID <- paste0("Mut", 1:nrow(map_z))
  }
  map_z <- map_z %>%
    mutate(Mut_ID = Mut_ID, Cluster = value) %>%
    select(Mut_ID, Cluster) %>%
    arrange(Cluster)
  # map_z <- map_z %>%
  #   add_column(Mut_ID = Mut_ID) %>%
  #   mutate(Cluster = value) %>%
  #   select(Mut_ID, Cluster) %>%
  #   arrange(Cluster)
  return(map_z)
}

#' Determine the most probable cluster CCF values by taking the mode of the posterior distributions
#' 
#' @export
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @return matrix of estimated cluster CCFs
estimateCCFs <- function(w_chain) {
  S <- numberSamples(w_chain)
  K <- numberClusters(w_chain)
  # density plot 
  w.dens <- ggplot(w_chain, aes(x = value)) +
    geom_density() +
    facet_wrap(~Parameter, ncol = S, scales = "free_y") +
    theme_light()
  # find peak for MAP w
  w.dens.p <- ggplot_build(w.dens)$data[[1]]
  w.map <- w.dens.p %>%
    as_tibble() %>%
    group_by(PANEL) %>%
    summarize(value = x[max(y) == y])
  w.map <- w.map %>%
    mutate(Parameter = unique(w_chain$Parameter),
           value_rounded = round(value, 2))
  # return w matrix
  w.map.matrix <- matrix(w.map$value_rounded, K, S, byrow=TRUE)
  return(w.map.matrix)
}

#' Determine the most probable cluster CCF values by taking the mode of the posterior distributions
#' 
#' @export
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @param Sample_ID Vector of sample IDs, same order as provided as input data (e.g. indata$Sample_ID)
#' @return A tibble of estimated cluster CCFs in each sample 
writeClusterCCFsTable <- function(w_chain, Sample_ID = NULL) {
  map_w <- as.data.frame(estimateCCFs(w_chain))
  
  if (is.null(Sample_ID)) {
    Sample_ID <- paste0("Sample ", 1:ncol(map_w))
  }
  colnames(map_w) <- Sample_ID
  map_w <- map_w %>%
    as_tibble() %>%
    bind_cols(tibble(Cluster = 1:nrow(map_w)), .)
  return(map_w)
}

# function to check if clustering respects sample presence 
test.pres <- function(samps, pres.tb) {
  # input: samps = one item in samps.list
  #        pres.tb = tibble with columns mut_id and sample_presence (label of sample-presence set)
  # output: TRUE or FALSE
  chains <- ggs(samps)
  z.chain <- get.parameter.chain("z", chains)
  map.z <- get.map.z(z.chain)
  
  # check that clusters respect sample presence
  pres.tb$cluster <- map.z$Cluster[match(pres.tb$mut_id, map.z$Mutation)]
  clusts <- unique(pres.tb$cluster)
  pres.equal <- all(sapply(clusts, function(x) 
    length(unique(pres.tb$sample_presence[pres.tb$cluster == x])) == 1))
}

plot_ppd <- function(samps, test.data, K) {
  ##
  ## 50 mutations x 10 samples
  chains <- do.call(rbind, samps)
  ystar <- chains[, grep("ystar", colnames(chains))]
  ## each row of MCMC is in column-major order
  orig.order <- tibble(statistic=colnames(ystar))
  ppd.summaries <- ystar %>%
      as_tibble() %>%
      gather("statistic", "value") %>%
      group_by(statistic) %>%
      summarize(mean=mean(value),
                q1=quantile(value, 0.025),
                q3=quantile(value, 0.975))
  ppd.summaries2 <- left_join(orig.order,
                              ppd.summaries, by="statistic") %>%
      mutate(observed=as.numeric(test.data$y)) %>%
      mutate(sample=paste0("sample", rep(1:test.data$S, each=test.data$I)),
             variant=rep(1:test.data$I, test.data$S))
  ggplot(ppd.summaries2, aes(x=statistic, y=mean,
                             ymin=q1,
                             ymax=q3)) +
      geom_errorbar() +
      geom_point(aes(x=statistic, y=observed),
                 size=1, color="steelblue") +
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.background=element_rect(fill="white",
                                          color="black")) +
      ylab("Middle 95% of posterior\npredictive distribution") +
      xlab("observation index (column-major order)") +
      facet_wrap(~sample) +
      ggtitle(paste0("K = ", K))
}


simulateVAF <- function(mcf, nvarClust, avg_depth=100, sd_depth=20){
    nClust <- nrow(mcf)
    nSamp <- ncol(mcf)
    nVariants <- sum(nvarClust)
    nObs <- nVariants * nSamp
    stopifnot(length(nvarClust) == nClust)
    ##pi <- rep(1/3, 3)
    ## 10 mutations for each of the 3 clusters
    z <- rep(1:nClust, nvarClust)

    MCF <- mcf[z, ]
    dimnames(MCF) <- list(paste0("variant", seq_len(nVariants)),
                          paste0("sample", seq_len(nSamp)))
    mult <- tcn <- MCF;
    tcn[,] <- 2
    mult[,] <- sample(1:2, nVariants*nSamp, replace=TRUE)
    vaf <- (mult * MCF)/(tcn * MCF + 2*(1-MCF))
    ## add more variation to sequencing depth
    mus <- rnorm(nVariants, avg_depth, sd_depth)
    mus <- ifelse(mus < 0, 0, mus)
    depth <- matrix(rpois(nObs, mus), nVariants, nSamp)
    y <- rbinom(nObs, as.numeric(depth), as.numeric(vaf))
    y <- ifelse(y < 0, 0, y)
    tibble(variant=rep(seq_len(nVariants), nSamp),
           sample=rep(seq_len(nSamp), each=nVariants),
           cluster=rep(z, nSamp),
           y=y,
           n=as.numeric(depth),
           multiplicity=as.numeric(mult),
           copy_number=as.numeric(tcn),
           mcf=as.numeric(MCF))
}


listJagInputs <- function(dat){
    w <- group_by(dat, cluster, sample) %>%
        summarize(mcf=unique(mcf)) %>%
        ungroup() %>%
        spread(sample, mcf) %>%
        select(-cluster) %>%
        as.matrix()
    K <- length(unique(dat$cluster))
    tcn <- dat %>%
        select(variant, sample, copy_number) %>%
        spread(sample, copy_number) %>%
        select(-variant) %>%
        as.matrix()
    S <- length(unique(dat$sample))
    m <- dat %>%
        select(variant, sample, multiplicity) %>%
        spread(sample, multiplicity) %>%
        select(-variant) %>%
        as.matrix()
    y <- dat %>%
        select(variant, sample, y) %>%
        spread(sample, y) %>%
        select(-variant) %>%
        as.matrix()
    n <- dat %>%
        select(variant, sample, n) %>%
        spread(sample, n) %>%
        select(-variant) %>%
        as.matrix()                
    I <- length(unique(dat$variant))

    list(I=I, S=S, K=K,
         y=y, n=n, m=m, tcn=tcn)
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

get.z.mapping <- function(true.z, z.chain) {
  mcmc_z <- z.chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              maxiter=max(Iteration)) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(value=value[probability==max(probability)]) %>%
    mutate(Variant = 1:length(true.z), true_z = true.z)
  z_mappings <- map_z %>%
    select(Variant, true_z, value) %>%
    rename("mcmc_z" = "value")
  z_mappings <- z_mappings %>%
    group_by(true_z) %>%
    summarize(mcmc_z = getmode(mcmc_z))
  z_mappings
}

relabel.w.z.chains <- function(true.z, chains) {
  w.chain <- get.parameter.chain("w", chains)
  z.chain <- get.parameter.chain("z", chains)
  
  z_mappings <- get.z.mapping(true.z, z.chain)
  
  # relabel w.chain
  w.chain.relabeled <- w.chain %>%
    mutate(k = as.numeric(gsub("w\\[", "", sapply(w.chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = gsub("\\]", "", sapply(w.chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][2])))
  w.chain.relabeled <- w.chain.relabeled %>%
    mutate(new_k = match(w.chain.relabeled$k, z_mappings$mcmc_z))
  w.chain.relabeled <- w.chain.relabeled %>%
    mutate(new_Parameter = paste0("w[", new_k, ",", s, "]")) %>%
    arrange(new_k, s)
  w.chain.relabeled <- w.chain.relabeled %>%
    select(Iteration, Chain, new_Parameter, value) %>%
    rename("Parameter" = "new_Parameter")
  w.chain.relabeled <- w.chain.relabeled %>%
    mutate(Parameter = factor(w.chain.relabeled$Parameter, levels = unique(w.chain.relabeled$Parameter)))
  
  # relabel z.chain
  z.chain.relabeled <- z.chain %>%
    mutate(new_value = match(value, z_mappings$mcmc_z))
  z.chain.relabeled <- z.chain.relabeled %>%
    select(Iteration, Chain, Parameter, new_value) %>%
    rename("value" = "new_value")
  
  list(w.chain=w.chain.relabeled,
       z.chain=z.chain.relabeled)
}

