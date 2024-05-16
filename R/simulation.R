simulation_type3_jags <- function(eta=0.8, S=3, K=2, I=50, num_cn=10, seed=12345, depth=50) {
  set.seed(12345)
  eta <- 0.85
  S <- 2
  K <- 3
  I <- 70
  num_cn=10
  depth = 50
  mcf <- ifelse(runif(K*S, 0, 1) < eta, 1, 0) *
    rbeta(K*S, 1, 1) %>%
    matrix(K, S)
  pi <- rdirichlet(1, rep(1, K))
  z <- replicate(I, sample(seq_len(K), size=1,
                           replace=TRUE, prob=pi))
  ## indicator for whether mutation is in copy-altered or LOH region
  h <- rep(c(0, 1), c(I-num_cn, num_cn))
  
  ## integer total copy number in the tumor
  icn <- rep(NA, I)
  icn[ h==0 ] <- 2
  icn[ h==1 ] <- rpois(sum(h==1), 2)
  
  ## m is the number of copies of the chromosome for an allele/variant
  m <- rep(NA, I)
  m[ h == 0 ] <- 1
  m[ h == 1 ] <- pmin(icn[h==1], rpois(sum(h==1), 1))
  
  indices <- sample(1:2, sum(h==1), replace = TRUE)
  m[ h==1 ] <- ifelse(indices == 1, m[h==1], icn[h==1]-m[h==1])
  
  ## tcn is the fraction copy number of mixture of normal and tumor cells
  epsilon <- 0.1
  tcn <- matrix(2, I, S)
  for (i in 1:I) {
    for (s in 1:S) {
      if (h[i] == 1) {
        tcn[i,s] <- rnorm(1, icn[i] * mcf[z[i], s] + 2 * (1-mcf[z[i], s]), epsilon)
      }
    }
  }
  
  ## indicator for whether a SSM is diploid
  ## cncf if the fraction of copy number alteration that overlaps a SSM
  q <- rep(0, I)
  for (i in 1:I) {
    if (i <= (I-num_cn)) {
      prob <- runif(1)
      if (prob > 0.9) {
        q[i] <- sample((I-num_cn+1):I,1)
        tcn[i,] <- tcn[q[i],]
        icn[i] <- icn[q[i]]
        m[i] <- sample(c(m[q[i]], icn[q[i]]-m[q[i]]),1)
      } else {
        q[i] <- i
      }
    } else {
      q[i] <- i
    }
  }
  
  vaf <- matrix(NA, I, S)

  for(i in seq_len(I)){
    for(s in seq_len(S)){
      if(h[i] == 0){
        vaf[i,s] <- max(0, (mcf[z[i], s] + (m[i] - 1) * mcf[z[q[i]],s]) / tcn[i,s])
        vaf[i,s] <- min(1, vaf[i,s])
      } else{
        vaf[i,s] <- max(0, (mcf[z[i], s]*m[i] + 1 - mcf[z[i], s]) / tcn[i,s])
        vaf[i,s] <- min(1, vaf[i,s])
      }
    }
  }
  
  n <- matrix(rpois(n = I*S, lambda = depth), I, S)
  y <- matrix(rbinom(n=I*S, size=as.numeric(n), prob=as.numeric(vaf)), I, S)
  
  model <- jags.model(file = "~/JHU/scripts/R/pictograph2/inst/extdata/type3.jags",
                      data = list(K=K, S=S, I=I, y=y, is_cn=h, n=n, tcn=tcn, q=q),
                      inits = list(".RNG.name" = "base::Wichmann-Hill",
                                   ".RNG.seed" = 123))
  
  update(model, n.iter = 1000)
  
  Nrep = 10000
  
  posterior_sample <- coda.samples(model,
                                   variable.names = c("mcf", "z", "icn", "m"),
                                   n.iter = Nrep)
  # plot(posterior_sample[,1])
  # summary(posterior_sample)
  
  # get.parameter.chain("z", ggmcmc::ggs(posterior_sample))
  
  chains <- formatChains(posterior_sample)
  
  mcf
  mcf_est <- writeClusterMCFsTable(chains$mcf_chain)
  mcf_est
  plotChainsCCF(chains$mcf_chain)
  
  # which(z==1)
  # which(z==2)
  # z
  z_est <- writeClusterAssignmentsTable(chains$z_chain)
  # z_est %>% filter(Cluster==2)
  z_est <- z_est %>% mutate(idx = as.numeric(gsub("Mut(\\d+)","\\1", Mut_ID))) %>% arrange(idx)
  # z_est$Cluster
  table(z, z_est$Cluster)
  
  m
  m_est <- writeMultiplicityTable(chains$m_chain)
  m_est$Multiplicity
  table(m, m_est$Multiplicity)
  # m_tmp <- get.parameter.chain("m", ggmcmc::ggs(posterior_sample))
  # m_tmp <- m_tmp %>% filter(Parameter=="m[46]")
  
  icn
  icn_est <- writeIcnTable(chains$icn_chain)
  icn_est$icn
  table(icn, icn_est$icn)
  # icn_tmp <- get.parameter.chain("icn", ggmcmc::ggs(posterior_sample))
  # icn_tmp <- icn_tmp %>% filter(Parameter=="icn[50]")
  # icn_tmp <- chains$icn_chain %>% group_by(Parameter) %>% reframe(val=round(mean(value))) %>% mutate(idx = as.numeric(gsub("icn\\[(\\d+)\\]","\\1", Parameter))) %>% arrange(idx)
}

simulation_type1_jags <- function(ets=0.8, S=3, K=2, I=50, num_cn=10, seed=12345) {
  set.seed(123)
  eta <- 0.8
  S <- 4
  K <- 3
  I <- 70
  num_cn=10
  mcf <- ifelse(runif(K*S, 0, 1) < eta, 1, 0) *
    rbeta(K*S, 1, 1) %>%
    matrix(K, S)
  pi <- rdirichlet(1, rep(1, K))
  z <- replicate(I, sample(seq_len(K), size=1,
                           replace=TRUE, prob=pi))
  ## indicator for whether mutation is a copy number alteration
  h <- rep(c(0, 1), c(I-num_cn, num_cn))
  ## integer total copy number in the tumor
  icn <- rep(NA, I)
  icn[ h==0 ] <- 2
  # multinom <- function(p){
  #   sample(0:(p-1), sum(h==1), replace=TRUE,
  #          prob=rep(1/p, p))
  # }
  # icn[ h==1 ] <- multinom(8)
  icn[ h==1 ] <- rpois(sum(h==1), 2)
  
  ## m is the number of copies of the chromosome for an allele/variant
  m <- rep(NA, I)
  m[ h == 0 ] <- 1
  # m[ h == 1 ] <- pmin(icn[h==1], multinom(4))
  m[ h == 1 ] <- pmin(icn[h==1], rpois(sum(h==1), 1))
    
  indices <- sample(1:2, sum(h==1), replace = TRUE)
  m[ h==1 ] <- ifelse(indices == 1, m[h==1], icn[h==1]-m[h==1])
  
  vaf <- matrix(NA, I, S)
  tcn <- matrix(2, I, S)
  epsilon <- 0.1
  for(i in seq_len(I)){
    for(s in seq_len(S)){
      if(h[i] == 0){
        vaf[i, s] <- mcf[z[i], s]/tcn[i,s]
      } else{
        tcn[i,s] <- max(0.01, rnorm(1, icn[i] * mcf[z[i], s] + 2 * (1-mcf[z[i], s]), epsilon))
        vaf[i,s] <- (mcf[z[i], s]*m[i] + 1 - mcf[z[i], s])/tcn[i,s]
      }
    }
  }
  
  n <- matrix(rpois(n = I*S, lambda = 50), I, S)
  y <- matrix(rbinom(n=I*S, size=as.numeric(n), prob=as.numeric(vaf)), I, S)

  model <- jags.model(file = "~/JHU/scripts/R/pictograph2/inst/extdata/type1.jags",
                      data = list(K=K, S=S, I=I, y=y, is_cn=h, n=n, tcn=tcn),
                      inits = list(".RNG.name" = "base::Wichmann-Hill",
                                   ".RNG.seed" = 123))
  
  update(model, n.iter = 1000)
  
  Nrep = 10000
  
  posterior_sample <- coda.samples(model,
                                   variable.names = c("mcf", "z", "icn", "m"),
                                   n.iter = Nrep)
  # plot(posterior_sample[,1])
  # summary(posterior_sample)
  
  # get.parameter.chain("z", ggmcmc::ggs(posterior_sample))
  
  chains <- formatChains(posterior_sample)
  
  mcf
  mcf_est <- writeClusterMCFsTable(chains$mcf_chain)
  mcf_est
  plotChainsMCF(chains$mcf_chain)
  
  # which(z==1)
  # which(z==2)
  # z
  z_est <- writeClusterAssignmentsTable(chains$z_chain)
  # z_est %>% filter(Cluster==2)
  z_est <- z_est %>% mutate(idx = as.numeric(gsub("Mut(\\d+)","\\1", Mut_ID))) %>% arrange(idx)
  # z_est$Cluster
  table(z, z_est$Cluster)

  m
  m_est <- writeMultiplicityTable(chains$m_chain)
  m_est$Multiplicity
  table(m, m_est$Multiplicity)
  # m_tmp <- get.parameter.chain("m", ggmcmc::ggs(posterior_sample))
  # m_tmp <- m_tmp %>% filter(Parameter=="m[46]")
  
  icn
  icn_est <- writeIcnTable(chains$icn_chain)
  icn_est$icn
  table(icn, icn_est$icn)
  # icn_tmp <- get.parameter.chain("icn", ggmcmc::ggs(posterior_sample))
  # icn_tmp <- icn_tmp %>% filter(Parameter=="icn[50]")
  # icn_tmp <- chains$icn_chain %>% group_by(Parameter) %>% reframe(val=round(mean(value))) %>% mutate(idx = as.numeric(gsub("icn\\[(\\d+)\\]","\\1", Parameter))) %>% arrange(idx)
}

simulation_type2_jags <- function(ets=0.8, S=3, K=2, I=50, num_cn=10, seed=12345) {
  set.seed(1234)
  eta <- 0.8
  S <- 3
  K <- 2
  I <- 50
  num_cn <- 10
  mcf <- ifelse(runif(K*S, 0, 1) < eta, 1, 0) *
    rbeta(K*S, 1, 1) %>%
    matrix(K, S)
  pi <- rdirichlet(1, rep(1, K))
  z <- replicate(I, sample(seq_len(K), size=1,
                           replace=TRUE, prob=pi))
  ## indicator for whether mutation is in copy-altered or LOH region
  h <- rep(c(0, 1), c(I-num_cn, num_cn))
  
  ## integer total copy number in the tumor
  icn <- rep(NA, I)
  icn[ h==0 ] <- 2
  icn[ h==1 ] <- rpois(sum(h==1), 2)
  
  ## m is the number of copies of the chromosome for an allele/variant
  m <- rep(NA, I)
  m[ h == 0 ] <- 1
  m[ h == 1 ] <- pmin(icn[h==1], rpois(sum(h==1), 1))
  
  indices <- sample(1:2, sum(h==1), replace = TRUE)
  m[ h==1 ] <- ifelse(indices == 1, m[h==1], icn[h==1]-m[h==1])
  
  
  ## tcn is the fraction copy number of mixture of normal and tumor cells
  epsilon <- 0.1
  tcn <- matrix(2, I, S)
  for (i in 1:I) {
    for (s in 1:S) {
      if (h[i] == 1) {
        tcn[i,s] <- rnorm(1, icn[i] * mcf[z[i], s] + 2 * (1-mcf[z[i], s]), epsilon)
      }
    }
  }
  
  ## indicator for whether a SSM is diploid
  ## cncf if the fraction of copy number alteration that overlaps a SSM
  q <- rep(0, I)
  cncf <- matrix(0, I, S)
  for (i in 1:I) {
    if (i <= (I-num_cn)) {
      prob <- runif(1)
      if (prob > 0.9) {
        q[i] <- sample((I-num_cn+1):I,1)
        if (icn[q[i]]!=0) {
          cncf[i,] <- mcf[z[q[i]],]
          icn[i] <- icn[q[i]]
          m[i] <- sample(c(m[q[i]], icn[q[i]]-m[q[i]]),1)
          tcn[i,] <- tcn[q[i],]
        } else {
          q[i] <- 0
        }
      } 
    } else {
      cncf[i,] <- mcf[z[i],]
    }
  }
  
  vaf <- matrix(NA, I, S)
  for(i in seq_len(I)){
    for(s in seq_len(S)){
      if(h[i] == 0){
        vaf[i, s] <- max(0, (mcf[z[i], s] + (m[i] - 1) * cncf[i,s]) / tcn[i,s])
      } else{
        vaf[i,s] <- max(0, (mcf[z[i], s]*m[i] + 1 - mcf[z[i], s]) / tcn[i,s])
      }
    }
  }
  
  n <- matrix(rpois(n = I*S, lambda = 50), I, S)
  y <- matrix(rbinom(n=I*S, size=as.numeric(n), prob=as.numeric(vaf)), I, S)
  
  model <- jags.model(file = "~/JHU/scripts/R/pictograph2/inst/extdata/type2.jags",
                      data = list(K=K, S=S, I=I, y=y, is_cn=h, n=n, tcn=tcn, mtp=m, cncf=cncf, icn=icn),
                      inits = list(".RNG.name" = "base::Wichmann-Hill",
                                   ".RNG.seed" = 123))
  
  update(model, n.iter = 1000)
  
  Nrep = 10000
  
  posterior_sample <- coda.samples(model,
                                   variable.names = c("mcf", "z", "icn", "m"),
                                   n.iter = Nrep)
  # plot(posterior_sample[,1])
  # summary(posterior_sample)
  
  # get.parameter.chain("z", ggmcmc::ggs(posterior_sample))
  
  chains <- formatChains(posterior_sample)
  
  mcf
  mcf_est <- writeClusterMCFsTable(chains$mcf_chain)
  mcf_est
  plotChainsCCF(chains$mcf_chain)
  
  # which(z==1)
  # which(z==2)
  # z
  z_est <- writeClusterAssignmentsTable(chains$z_chain)
  # z_est %>% filter(Cluster==2)
  z_est <- z_est %>% mutate(idx = as.numeric(gsub("Mut(\\d+)","\\1", Mut_ID))) %>% arrange(idx)
  # z_est$Cluster
  table(z, z_est$Cluster)
  
  m
  m_est <- writeMultiplicityTable(chains$m_chain)
  m_est$Multiplicity
  table(m, m_est$Multiplicity)
  # m_tmp <- get.parameter.chain("m", ggmcmc::ggs(posterior_sample))
  # m_tmp <- m_tmp %>% filter(Parameter=="m[46]")
  
  icn
  icn_est <- writeIcnTable(chains$icn_chain)
  icn_est$icn
  table(icn, icn_est$icn)
  # icn_tmp <- get.parameter.chain("icn", ggmcmc::ggs(posterior_sample))
  # icn_tmp <- icn_tmp %>% filter(Parameter=="icn[50]")
  # icn_tmp <- chains$icn_chain %>% group_by(Parameter) %>% reframe(val=round(mean(value))) %>% mutate(idx = as.numeric(gsub("icn\\[(\\d+)\\]","\\1", Parameter))) %>% arrange(idx)
}


simulation <- function(mcf=0, icn=2, minor_cn=1, depth=30, num_SNV=30, seed=NULL, normal_prop=0) {
  # mcf=0.5
  # icn=4
  # minor_cn=1
  # depth=100
  # num_SNV=100
  # normal_prop=0 # proportion of SNV is actually from CN-neutral region
  # seed=123

  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  tcn = 2 * (1-mcf) + icn * mcf
  
  # if a proportion of SNVs on CNA is actually on CN-neutral region
  num_neutral = round(num_SNV * normal_prop)
  SNV_depth_neutral = rpois(n = num_neutral, lambda = depth)
  SNV_alt_neutral = numeric(num_neutral)
  for (i in seq_len(num_neutral)) {
    SNV_alt_neutral[i] <- rbinom(n=1, size=SNV_depth_neutral[i], prob=0.5)
  }
  
  vaf_neutral = SNV_alt_neutral / SNV_depth_neutral
  
  # Actual tumor proportion
  num_SNV = num_SNV - num_neutral
  
  # generate depth for each SNV
  depth_total = rpois(n = num_SNV, lambda = depth * tcn / 2)
  
  # assign each SNV to one copy
  SNV_assignment = sample(c(1, 2), size = num_SNV, replace = TRUE) # which segment
  
  # generate depth for germline
  SNV_depth_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_depth_germline[i] <- rbinom(n=1, size=depth_total[i], prob=2 * (1-mcf)/tcn)
  }
  
  SNV_alt_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_alt_germline[i] <- rbinom(n=1, size=SNV_depth_germline[i], prob=0.5)
  }
  
  vaf_germline <- SNV_alt_germline/SNV_depth_germline
  
  # generate depth for tumor
  SNV_depth_tumor = depth_total - SNV_depth_germline 
  SNV_alt_tumor = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    if (icn == 0) {
      SNV_alt_tumor[i] <- 0
    } else {
      tumor_vaf = minor_cn / icn
      if (SNV_assignment[i] == 1) {
        tumor_vaf = 1 - (minor_cn / icn)
      } 
      SNV_alt_tumor[i] <- rbinom(n=1, size=SNV_depth_tumor[i], prob=tumor_vaf)
    }
  }
  
  # vaf_tumor <- ifelse(SNV_depth_tumor==0, 0, SNV_alt_tumor/SNV_depth_tumor)
  vaf_tumor = c((SNV_alt_germline + SNV_alt_tumor) / depth_total, vaf_neutral)
  # return(vaf_neutral)
  title = paste("mcf: ", mcf, ", num_SNV: ", num_SNV, ", depth: ", depth, ", \nICN: ", icn, ", minor CN: ", minor_cn, 
                ", true_prop: ", 1-normal_prop, ", unimodal: ", is.unimodal(vaf_tumor), sep="")
  plot(density(vaf_tumor), xlim=c(0,1), main = title, cex.main=0.9)

  # plot(density(vaf_tumor), xlim=c(0,1), main = "VAF_tumor")
  # qqnorm(vaf_tumor, main = "tumor")
  # qqline(vaf_tumor, col="grey")
  
  is.unimodal(vaf_tumor)
  # return(um)
}

simulation_normal <- function(depth=30, num_SNV=30) {
  
  # generate depth for each SNV
  SNV_depth = rpois(n = num_SNV, lambda = depth)

  # generate alt read counts for each SNV
  SNV_alt = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_alt[i] <- rbinom(n=1, size=SNV_depth[i], prob=0.5)
  }
  
  # calculate VAF for all SNVs
  vaf <- SNV_alt/SNV_depth
  
  # plot the SNVs
  title = paste("num_SNV: ", num_SNV, ", depth: ", depth, ", unimodal: ", is.unimodal(vaf), sep="")
  plot(density(vaf), xlim=c(0,1), main = title)
}

simulation_CNA_germline_SNV <- function(mcf=0.9, icn=2, minor_cn=1, depth=30, num_SNV=30, seed=NULL, normal_prop=1) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  tcn = 2 * (1-mcf) + icn * mcf
  
  # if a proportion of SNVs on CNA is actually on CN-neutral region
  num_neutral = round(num_SNV * normal_prop)
  SNV_depth_neutral = rpois(n = num_neutral, lambda = depth)
  SNV_alt_neutral = numeric(num_neutral)
  for (i in seq_len(num_neutral)) {
    SNV_alt_neutral[i] <- rbinom(n=1, size=SNV_depth_neutral[i], prob=0.5)
  }
  
  vaf_neutral = SNV_alt_neutral / SNV_depth_neutral
  
  # Actual tumor proportion
  num_SNV = num_SNV - num_neutral
  
  # generate depth for each SNV
  depth_total = rpois(n = num_SNV, lambda = depth * tcn / 2)
  
  # assign each SNV to one copy
  SNV_assignment = sample(c(1, 2), size = num_SNV, replace = TRUE) # which segment
  
  # generate depth for germline
  SNV_depth_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_depth_germline[i] <- rbinom(n=1, size=depth_total[i], prob=2 * (1-mcf)/tcn)
  }
  
  SNV_alt_germline = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    SNV_alt_germline[i] <- rbinom(n=1, size=SNV_depth_germline[i], prob=0.5)
  }
  
  vaf_germline <- SNV_alt_germline/SNV_depth_germline
  
  # generate depth for tumor
  SNV_depth_tumor = depth_total - SNV_depth_germline 
  SNV_alt_tumor = numeric(num_SNV)
  for (i in seq_len(num_SNV)) {
    if (icn == 0) {
      SNV_alt_tumor[i] <- 0
    } else {
      tumor_vaf = minor_cn / icn
      if (SNV_assignment[i] == 1) {
        tumor_vaf = 1 - (minor_cn / icn)
      } 
      SNV_alt_tumor[i] <- rbinom(n=1, size=SNV_depth_tumor[i], prob=tumor_vaf)
    }
  }
  
  vaf_tumor = c((SNV_alt_germline + SNV_alt_tumor) / depth_total, vaf_neutral)
  
  is.unimodal(vaf_tumor)
}

