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

  # ---------------------------------------------------------------------------
  ## try fitting mixture to pattern where columns 2 and 3 zero 
  devtools::load_all()
  set.seed(1)
  sim.data2 <- simTestCaseZeros(20, avg.cov=100)
  input.data2 <- sim.data2[c("y", "n", "purity", "tcn", "m", "I", "S")]
  kToTest <- 10
  samps.list <- mclapply(kToTest,
                         function(k) runMCMC(input.data2, k, jags.file, inits, params,
                                             n.iter=1000, thin=1,
                                             n.burn=100),
                         mc.cores=8)
  # ---------------------------------------------------------------------------
  # check Z assignment for K=10
  chains <- ggs(samps.list[[which(kToTest == 10)]])
  
  # are there really 10 MAP clusters?
  map_z <- get.map.z(get.parameter.chain("z", chains))
  length(unique(map_z$value))
  expect_equivalent(length(unique(map_z$value)), 10)
  
  # compare cluster assignment (z) to truth
  # map cluster numbers back to truth (for K=10)
  w.z.relabeled.chains <- relabel.w.z.chains(sim.data2$z, chains)
  w.chain <- w.z.relabeled.chains[["w.chain"]]
  z.chain <- w.z.relabeled.chains[["z.chain"]]
  map_z_relabeled <- get.map.z(z.chain)
  
  expect_equivalent(map_z_relabeled$value, sim.data2$z)
  
  zs <- map_z_relabeled %>%
    mutate(truez=sim.data2$z)
  table(zs$value, zs$truez)
  ## this fails
  # ---------------------------------------------------------------------------
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
  samps.list2 <- mclapply(kToTest,
                         function(k) runMCMC(input.data2.clust49, k, jags.file, inits, params,
                                             n.iter=1000, thin=1,
                                             n.burn=100),
                         mc.cores=8)
  # need different jags model for K=1
  jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1.jags")
  samps.K1 <- runMCMC(input.data2.clust49, 1, jags.file.K1, inits, params, n.iter=1000, thin=1, n.burn=100)
  samps.list2 <- c(list(samps.K1), samps.list2)
  # check BIC
  BIC <- mapply(function(samps, k)
    calcBIC(input.data2.clust49$I*input.data2.clust49$S, k, calcChainLogLik(samps, input.data2.clust49, k)),
    samps = samps.list2, k = 1:maxK)
  min.BIC.k <- which(BIC == min(BIC))
  min.BIC.k # true K = 2, also fails
  
  # ---------------------------------------------------------------------------
  # # check Z assignment for K=2
  # chains2 <- ggs(samps.list2[[2]])
  # 
  # # are there really 2 MAP clusters?
  # map_z2 <- get.map.z(get.parameter.chain("z", chains2))
  # expect_equivalent(length(unique(map_z2$value)), 2)
  # 
  # # compare cluster assignment (z) to truth
  # zs2 <- map_z2 %>%
  #   mutate(truez=sim.data2$z[select.muts])
  # table(zs2$value, zs2$truez)
  # 
  # # check w 
  # w.map2 <- get.map.w(get.parameter.chain("w", chains2))
  # w.map2
  # ---------------------------------------------------------------------------
  
  # BIC_tb <- tibble(k = kToTest,
  #                  BIC = BIC)
  # ggplot(BIC_tb, aes(x = k, y = BIC)) +
  #   geom_line() +
  #   geom_hline(yintercept=min(BIC), lty="dashed", colour="darkgrey") +
  #   theme_light() +
  #   scale_x_continuous(breaks = kToTest)

  # compare cancer cell fraction (w) values to truth
  w.map.matrix <- get.map.w(w.chain)
  thresh <- 0.025
  w.diff <- abs(sim.data$w - w.map.matrix)
  for (k in 1:sim.data$K) {
    for (s in 1:sim.data$S) {
      expect_lte(w.diff[k,s], thresh)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Test case using IP30 CCFs for K=19
  sim.data3 <- simTestCaseIP30(30, avg.cov=100) 
  input.data3 <- sim.data3[c("y", "n", "purity", "tcn", "m", "I", "S")]
  
  inits <- list(".RNG.name" = "base::Wichmann-Hill",
                ".RNG.seed" = 123)
  params <- c("z", "w")
  kToTest <- 10:22
  #kToTest <- 19
  samps.list3 <- mclapply(kToTest,
                          function(k) runMCMC(input.data3, k, jags.file, inits, params,
                                              n.iter=1000, thin=1,
                                              n.burn=100),
                          mc.cores=8)
  # check if best K (min BIC) is same as truth (K=10)
  BIC3 <- mapply(function(samps, k)
    calcBIC(input.data3$I*input.data3$S, k, calcChainLogLik(samps, input.data3, k)),
    samps = samps.list3, k = kToTest)
  min.BIC.k <- kToTest[which(BIC3 == min(BIC3))]
  min.BIC.k
  
  
})
