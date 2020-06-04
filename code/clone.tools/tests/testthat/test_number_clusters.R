context("Number of clusters")

test_that("Best K as min BIC", {
  library(rjags)
  library(ggmcmc)
  library(stringr)
  library(parallel)
  
  set.seed(123)
  extdir <- system.file("extdata", package="clone.tools")
  jags.file <- file.path(extdir, "spike_and_slab_purity_2.jags")
  #sim.data <- simTestCase2(5)
  #sim.data <- simulateDataPurity() # this sim should pass all tests
  sim.data <- simTestCase1(50, avg.cov=100) # 10 (avg.cov=100) fails, avg.cov=200 passes all
  input.data <- sim.data[c("y", "n", "purity", "tcn", "m", "I", "S")]
  
  inits <- list(".RNG.name" = "base::Wichmann-Hill",
                ".RNG.seed" = 123)
  params <- c("z", "w")
  kToTest <- 5:15
  kToTest <- 10
  samps.list <- mclapply(kToTest,
                         function(k) runMCMC(input.data, k, jags.file, inits, params,
                                           n.iter=1000, thin=1,
                                           n.burn=100),
                         mc.cores=8)
  
  # check if best K (min BIC) is same as truth (K=10)
  BIC <- mapply(function(samps, k)
    calcBIC(input.data$I*input.data$S, k, calcChainLogLik(samps, input.data, k)),
    samps = samps.list, k = kToTest)
  min.BIC.k <- kToTest[which(BIC == min(BIC))]
  expect_equivalent(10, min.BIC.k)


  ## try fitting mixture to pattern where columns 2 and 3 zero 
  devtools::load_all()
  sim.data2 <- simTestCaseZeros(20, avg.cov=100)
  input.data2 <- sim.data2[c("y", "n", "purity", "tcn", "m", "I", "S")]
  kToTest <- 10
  samps.list <- mclapply(kToTest,
                         function(k) runMCMC(input.data2, k, jags.file, inits, params,
                                             n.iter=1000, thin=1,
                                             n.burn=100),
                         mc.cores=8)
  ## this fails
  ## try K 1 - smallnumber (4)  on just the mutations with 2 columns (S1, S2) of zeros (clusters 4 and 9)
  input.data2.clust49 <- input.data2
  select.muts <- rowSums(input.data2.clust49$y[,1:2]) == 0
  input.data2.clust49$I <- sum(select.muts)
  input.data2.clust49$y <- input.data2.clust49$y[select.muts, ]
  input.data2.clust49$n <- input.data2.clust49$n[select.muts, ]
  input.data2.clust49$tcn <- input.data2.clust49$tcn[select.muts, ]
  input.data2.clust49$m <- input.data2.clust49$m[select.muts, ]
  maxK <- 4
  kToTest <- 2:maxK
  samps.list <- mclapply(kToTest,
                         function(k) runMCMC(input.data2.clust49, k, jags.file, inits, params,
                                             n.iter=1000, thin=1,
                                             n.burn=100),
                         mc.cores=8)
  # need different jags model for K=1
  jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1.jags")
  samps.K1 <- runMCMC(input.data2.clust49, 1, jags.file.K1, inits, params, n.iter=1000, thin=1, n.burn=100)
  samps.list <- c(list(samps.K1), samps.list)
  # check BIC
  BIC <- mapply(function(samps, k)
    calcBIC(input.data2.clust49$I*input.data2.clust49$S, k, calcChainLogLik(samps, input.data2.clust49, k)),
    samps = samps.list, k = 1:maxK)
  min.BIC.k <- which(BIC == min(BIC))
  min.BIC.k
  expect_equivalent(2, min.BIC.k)
  
  
  # BIC_tb <- tibble(k = kToTest,
  #                  BIC = BIC)
  # ggplot(BIC_tb, aes(x = k, y = BIC)) +
  #   geom_line() +
  #   geom_hline(yintercept=min(BIC), lty="dashed", colour="darkgrey") +
  #   theme_light() +
  #   scale_x_continuous(breaks = kToTest)
  
  # ---------------------------------------------------------------------------
  # chains for K=10
  chains <- ggs(samps.list[[which(kToTest == 10)]])
  
  # are there really 10 MAP clusters?
  map_z <- get.map.z(get.parameter.chain("z", chains))
  expect_equivalent(length(unique(map_z$value)), 10)

  # compare cluster assignment (z) to truth
  # map cluster numbers back to truth (for K=10)
  w.z.relabeled.chains <- relabel.w.z.chains(sim.data2$z[select.muts], chains)
  w.chain <- w.z.relabeled.chains[["w.chain"]]
  z.chain <- w.z.relabeled.chains[["z.chain"]]
  map_z_remapped <- get.map.z(z.chain)
  
  expect_equivalent(map_z_remapped$value, sim.data$z)

  zs <- map_z_remapped %>%
      mutate(truez=sim.data$z)
  table(zs$value, zs$truez)

  
  
  # compare cancer cell fraction (w) values to truth
  # MAP w -------------------
  # density plot 
  w.dens <- ggplot(w.chain, aes(x = value)) +
    geom_density() +
    facet_wrap(~Parameter, ncol = sim.data$S, scales = "free_y") +
    theme_light()
  # find peak for MAP w
  w.dens.p <- ggplot_build(w.dens)$data[[1]]
  w.map <- w.dens.p %>%
    as_tibble() %>%
    group_by(PANEL) %>%
    summarize(value = x[max(y) == y])
  w.map <- w.map %>%
    mutate(Parameter = unique(w.chain$Parameter),
           value_rounded = round(value, 2))
  # is MAP w within certain range of truth?
  w.map.matrix <- matrix(w.map$value_rounded, sim.data$K, sim.data$S, byrow=TRUE)
  thresh <- 0.025
  w.diff <- abs(sim.data$w - w.map.matrix)
  for (k in 1:sim.data$K) {
    for (s in 1:sim.data$S) {
      expect_lte(w.diff[k,s], thresh)
    }
  }
})
