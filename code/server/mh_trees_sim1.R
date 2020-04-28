#!/usr/bin/env Rscript
suppressMessages(library(rjags))
suppressMessages(library(coda))
suppressMessages(library(ggmcmc))
suppressMessages(library(optparse))
suppressMessages(library(gridExtra))
suppressMessages(library(stringr))
suppressMessages(library(fs))
suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
suppressMessages(library(GGally))
suppressMessages(library(network))
suppressMessages(library(sna))

pictograph_path <- "~/GitHub/pictograph"

devtools::load_all(file.path(pictograph_path, "code", "clone.tools"))

data.dir <- file.path(pictograph_path, "data")

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-s", "--seed"), type="numeric", default=123, 
              help="random seed for tree MH", metavar="numeric"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory", metavar="character"),
  make_option(c("-i", "--numIter"), type="numeric", default=50000,
  			  help="number of iterations for tree MH", metavar="numeric"),
  make_option(c("-t", "--thin"), type="numeric", default=10, 
              help="the thinning interval between consecutive observations", metavar="numeric")
)
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory must be supplied", call.=FALSE)
}

dir_create(opt$outdir)
# ===================================================================================== #
# Data

I <- 120; K <- 10; S <- 3
input.data <- readRDS(file.path(pictograph_path, "output",
                                "cluster_variants_BIC_spikenslab.Rmd",
                                "input.data.rds"))
w.chain <- readRDS(file.path(pictograph_path, "output",
                            "cluster_variants_BIC_spikenslab.Rmd",
                            "K10_w_chain_relabeled.rds"))
z.chain <- readRDS(file.path(pictograph_path, "output",
                            "cluster_variants_BIC_spikenslab.Rmd",
                            "K10_z_chain_relabeled.rds"))

truth <- input.data$w
mcmc_z <- z.chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              maxiter=max(Iteration)) %>%
    mutate(probability=n/maxiter)

mcf_stats <- w.chain %>%
    group_by(Parameter) %>%
    summarize(sd=sd(value),
              mean=mean(value))
map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(value=value[probability==max(probability)])

# ===================================================================================== #
# Tree MH

set.seed(opt$seed)

sampled.w.chain <- list(sample.w(w.chain, K))
init_admat <- initializeAdjacencyMatrix(mcf_stats = mcf_stats, 
                                        mcf_matrix = sampled.w.chain[[1]])
admat.chain <- list(rand.admat(init_admat))
cpov.temp <- create.cpov(mcf_stats = mcf_stats, mcf_matrix = sampled.w.chain[[1]])
#cpov.chain <- list(create.cpov(mcf_stats = mcf_stats, mcf_matrix = sampled.w.chain[[1]]))
score.chain <- calc.tree.fitness(admat.chain[[1]], cpov.temp, sampled.w.chain[[1]])

numAccept  <-  0
ncol.to.mutate <- 1 

thin <- opt$thin

admat.prev <- admat.chain[[1]]
fit.prev <- score.chain[1]

for (i in seq_len(opt$numIter)) {
	#print(i)
	
	# propse new admat
    sampled.w.star <- sample.w(w.chain, K)
    admat.star <- mutate.admat.2(admat.prev, ncol.to.mutate)
	print(admat.star)
    cpov.star <- create.cpov(mcf_stats, mcf_matrix = sampled.w.star)
    fit.star <- calc.tree.fitness(admat.star, cpov.star, sampled.w.star)

    r <- fit.star / fit.prev
    u <- runif(1, 0, 1)

    if(u <= r) {
        admat.prev <- admat.star
        numAccept <- numAccept + 1
        fit.prev <- fit.star
        #cpov.chain[[i+1]] <- cpov.star
        #sampled.w.chain[[i+1]] <- sampled.w.star
    } # else stay the same

    # keep thinned obs
    if (i%%opt$thin == 0) {
        admat.chain <- append(admat.chain, list(admat.prev))
        score.chain <- c(score.chain, fit.prev)
    } 
}
results <- list(admat.chain=admat.chain,
                score.chain=score.chain,
                #cpov.chain=cpov.chain,
                #sampled.w.chain=sampled.w.chain,
                #numAccept=numAccept,
                acceptRate=numAccept/opt$numIter)
print(paste0("acceptance rate = ", results$acceptRate))

saveRDS(results, file.path(opt$outdir, paste0("trees_numIter", opt$numIter, "_thin", opt$thin, "_seed", opt$seed, ".rds")))

