context("simulation")


unit_test_data <- function(){
    set.seed(123)
    mcf <- matrix(c(0.98, 0.99, 0.97, 
                    0.55, 0.00, 0.80, 
                    0.30, 0.70, 0.00,
                    0.20, 0.22, 0.18),
                  byrow=TRUE,
                  nrow=4, ncol=3)
    dimnames(mcf) <- list(paste0("cluster", 1:4),
                          paste0("sample", 1:3))
    ## number variants per cluster
    nvarClust <- c(4, 10, 7, 15)
    dat <- simulateVAF(mcf, nvarClust)
    dat
}

test_that("simulateVAF", {
    ##
    ## mutation cell fraction
    ## - fraction of cells with mutation
    ## - affected by purity
    ##
    ## cancer cell fraction
    ## - fraction of cancer cells with mutation
    ##
    set.seed(123)
    mcf <- matrix(c(0.98, 0.99, 0.97, 
                    0.55, 0.00, 0.80, 
                    0.30, 0.70, 0.00,
                    0.20, 0.22, 0.18),
                  byrow=TRUE,
                  nrow=4, ncol=3)
    dimnames(mcf) <- list(paste0("cluster", 1:4),
                          paste0("sample", 1:3))
    ## number variants per cluster
    nvarClust <- c(4, 10, 7, 15)
    dat <- simulateVAF(mcf, nvarClust)
    nobs <- sum(nvarClust) * ncol(mcf)
    expect_true(nobs==nrow(dat))
})
        
test_that("jags_sampler", {
    library(rjags)
    library(ggmcmc)
    extdir <- system.file("extdata", package="clone.tools")
    jags.file <- file.path(extdir, "mcf_model.jag")
    dat <- unit_test_data()
    jags_inputs <- listJagInputs(dat)
    model <- jags.model(jags.file,
                        data=jags_inputs,
                        n.chains=3,
                        n.adapt=500)
    samples <- coda.samples(model, n.iter=1000, thin=1,
                            variable.names=c("w", "z", "ystar"))
    tib <- ggs(samples)

    w.chain <- ggs(samples, family="w") %>%
        mutate(cluster_index=clusterIndex(Parameter),
               sample_index=sampleIndex(Parameter))
    z.chain <- ggs(samples, family="z") %>%
        mutate(variant_index=mutationIndex(Parameter))
    expect_true(length(unique(z.chain$variant_index)) == jags_inputs$I)
})

