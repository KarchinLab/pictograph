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
  input.data2 <- sim.data[c("y", "n", "purity", "tcn", "m", "I", "S")]
  kToTest <- 10
  samps.list <- mclapply(kToTest,
                         function(k) runMCMC(input.data2, k, jags.file, inits, params,
                                             n.iter=1000, thin=1,
                                             n.burn=100),
                         mc.cores=8)
  ## this fails
  ## try K 1 - smallnumber (4)  on just the mutations with 2 columns of zeros
  
  
  # BIC_tb <- tibble(k = kToTest,
  #                  BIC = BIC)
  # ggplot(BIC_tb, aes(x = k, y = BIC)) +
  #   geom_line() +
  #   geom_hline(yintercept=min(BIC), lty="dashed", colour="darkgrey") +
  #   theme_light() +
  #   scale_x_continuous(breaks = kToTest)
  
  # map cluster numbers back to truth (for K=10)
  chains <- ggs(samps.list[[which(kToTest == 10)]])
  
  # are there really 10 MAP clusters?
  z.chain <- get.parameter.chain("z", chains)
  it <- max(z.chain$Iteration)
  mcmc_z <- z.chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              maxiter=it) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(value=value[probability==max(probability)])
  expect_equivalent(length(unique(map_z$value)), 10)

  # compare cluster assignment (z) to truth
  w.z.relabeled.chains <- relabel.w.z.chains(sim.data$z, chains)
  w.chain <- w.z.relabeled.chains[["w.chain"]]
  
  mcmc_z_remapped <- w.z.relabeled.chains[["z.chain"]] %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              maxiter=max(Iteration)) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_z_remapped <- mcmc_z_remapped %>%
    group_by(Parameter) %>%
    summarize(value=value[probability==max(probability)])
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
