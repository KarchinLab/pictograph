context("Number of clusters")

test_that("Best K as min BIC", {
  library(rjags)
  library(ggmcmc)
  library(stringr)
  library(parallel)
  
  set.seed(123)
  extdir <- system.file("extdata", package="clone.tools")
  jags.file <- file.path(extdir, "spike_and_slab_purity_2.jags")
  sim.data <- simTestCase2(5)
  #sim.data <- simulateDataPurity() # this sim should pass all tests
  #sim.data <- simTestCase1(3) # 10 passes all tests, 3 fails
  input.data <- sim.data[c("y", "n", "purity", "tcn", "m", "I", "S")]
  
  inits <- list(".RNG.name" = "base::Wichmann-Hill",
                ".RNG.seed" = 123)
  params <- c("z", "w")
  kToTest <- 5:15
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
  
  # compare cluster assignment to truth
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
  expect_equivalent(length(unique(map_z_remapped$value)), 10)
  
  
})