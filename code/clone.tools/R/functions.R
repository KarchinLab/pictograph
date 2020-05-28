calcTheta <- function(m, tcn, w) {
    (m * w) / (tcn * w + 2*(1-w))
}

calcTheta2 <- function(m, tcn, w, p) {
  (m * w * p) / (tcn * p + 2*(1-p))
}

runMCMC <- function(data, K, jags.file, inits, params,
                    n.iter=20000, thin=10, n.chains=1,
                    n.adapt=1000, n.burn=10000) {
    data$K <- K
    jags.m <- jags.model(jags.file,
                         data,
                         n.chains = n.chains,
                         inits = inits,
                         n.adapt = n.adapt)
    update(jags.m, n.burn)
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

calcLogLik <- function(z.iter, w.iter, input.data) {
  W <- w.iter[z.iter, ]
  
  if (is.null(input.data$purity)) {
    theta <- calcTheta(input.data$m, input.data$tcn, W)
  } else {
    purity <- input.data$purity
    P <- matrix(rep(purity, each = I), nrow = nrow(W), ncol = ncol(W))
    theta <- calcTheta2(input.data$m, input.data$tcn, W, P)
  }
  
  sum(dbinom(as.matrix(input.data$y), as.matrix(input.data$n), as.matrix(theta), log=T))
}

calcChainLogLik <- function(samps, input.data, K) {
  z.chain <- getParamChain(samps, "z\\[")
  w.chain <- getParamChain(samps, "w\\[")
  lik <- c()
  for(iter in 1:nrow(z.chain)) {
    z.iter <- z.chain[iter, ]
    w.iter <- reshapeW(w.chain[iter, ], input.data$S, K)
    lik <- c(lik, calcLogLik(z.iter, w.iter, input.data))
  }
  mean(lik)
}

calcBIC <- function(n, k, ll) log(n)*k - 2*ll

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



numberClusters <- function(mcf_stats){
    params <- as.character(mcf_stats$Parameter)
    K <- strsplit(params, ",") %>%
        sapply("[", 1) %>%
        str_replace("w\\[", "") %>%
        as.numeric() %>%
        max()
    K
}

numberSamples <- function(mcf_stats){
    params <- as.character(mcf_stats$Parameter)    
    nSamples <- strsplit(params, ",") %>%
        sapply("[", 2) %>%
        str_replace("\\]", "") %>%
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

get.map.z <- function(z.chain) {
  numIter <- max(z.chain$Iteration)
  mcmc_z <- z.chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              numIter=numIter) %>%
    mutate(probability=n/numIter) %>%
    ungroup()
  
  map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    select(Parameter, value, probability)
  map_z <- filter(map_z, probability==max(probability)) %>%
    ungroup()
  map_z <- map_z %>%
    mutate(variant_ind = as.numeric(
      gsub("z\\[", "", gsub("\\]", "", as.character(map_z$Parameter)))))
  
  muts <- y$mut_id
  muts.z <- map_z %>%
    mutate(Mutation = muts[map_z$variant_ind]) %>%
    arrange(value) %>%
    rename(Cluster = value, Posterior_probability = probability)
  muts.z <- muts.z %>%
    select(Cluster, Posterior_probability, Mutation)
  muts.z
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

#' @export
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


#' @export
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



