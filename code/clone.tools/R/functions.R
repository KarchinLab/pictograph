calcTheta <- function(m, tcn, w) {
  (m * w) / (tcn * w + 2*(1-w))
}

runMCMC <- function(data, K, jags.file, inits, params,
                    n.iter, thin, n.chains=1,
                    n.adapt=1000) {
    data$K <- K

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

calcLogLik <- function(z.iter, w.iter, input.data) {
  W <- w.iter[z.iter, ]
  theta <- calcTheta(input.data$m, input.data$tcn, W)
  sum(dbinom(as.matrix(input.data$y), as.matrix(input.data$n), as.matrix(theta), log=T))
}

calcChainLogLik <- function(samps, input.data, K) {
  z.chain <- getParamChain(samps, "z\\[")
  w.chain <- getParamChain(samps, "w")
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

