library(ggnet)
library(network)
library(sna)
##library(rjags)
library(coda)
library(ggmcmc)
library(dplyr)
library(gridExtra)
library(stringr)
library(fs)
devtools::load_all("clone.tools")
outdir <- file.path("..", "output", "mh_trees.R")
dir_create(outdir)

I <- 100; K <- 10; S <- 3
set.seed(123)
test.data <- simulateData(I, K, S)
chains <- readRDS(file.path("..",
                            "output",
                            "cluster_variants.Rmd",
                            "chains.rds"))
w.chain <- chains[grep("w", chains$Parameter), ]
truth <- test.data$w
z.chain <- chains[grep("z", chains$Parameter), ]
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

set.seed(1234)
## construct adjacency matrix
answer <- initializeAdjacencyMatrix(mcf_stats, 0.01)
answer["root", "clone7"] <- answer["clone7", "clone1"] <-
    answer["clone1", "clone8"] <- answer["clone1", "clone9"] <-
    answer["clone8", "clone10"] <- answer["clone8", "clone2"] <-
    answer["clone10", "clone4"] <- answer["clone2", "clone3"] <-
    answer["clone9", "clone6"] <- answer["clone6", "clone5"] <- 1
truth <- plotDAG(answer)
##
## Each clone has one and only one parent. The parent can be the root.  Practically, this implies that the columns of the adjacency matrix sum to one.
##
stopifnot(all(colSums(answer, na.rm=TRUE) == 1))
## cluster precedence order violation matrix
cpov <- create.cpov(mcf_stats)
##trace(calc.tree.fitness, browser)
##calc.tree.fitness(answer, cpov, mcmc_w)
##trace(calc.topology.cost, browser)
TC <- calc.topology.cost(answer, cpov)
MC <- calc.mass.cost(answer, mcfMatrix(mcf_stats))
calc.tree.fitness(answer, cpov, mcf_stats)
##calc.tree.fitness(answer, cpov, w[c(1,2,9,5,10,8,6,3,7,4),])
##admat.chain <- list(init.admat(mcmc_w, zero.thresh=0.01))
## start at best
if(FALSE){
    best.admat.mcmc <- base.admat(mcmc_w)
    best.admat.mcmc[1,1] <- best.admat.mcmc[2,2] <-
        best.admat.mcmc[3,7] <- best.admat.mcmc[3,8] <-
        best.admat.mcmc[5,6] <- best.admat.mcmc[8,9] <- 
        best.admat.mcmc[9,4] <- best.admat.mcmc[9,10] <-
        best.admat.mcmc[10,5] <- best.admat.mcmc[11,3] <- 1
}
## shouldn't the truth be the best?'
best.admat.mcmc <- answer
admat.chain <- list(best.admat.mcmc)
numAccept  <-  0
ncol.to.mutate <- 1 
numIter <- 1000
score.chain <- c()
##N <- 10000
##THIN <- 20
N <- 10e3; THIN <- 250
probs <- rep(0, N)
accept <- rep(0, N)
for (i in seq_len(N)) {
    for(j in seq_len(THIN)){
        fit.prev <- calc.tree.fitness(admat.chain[[i]], cpov, mcf_stats)
        score.chain[i] <- fit.prev
        ##
        ## Propose edge
        ##
        admat.star <- mutate.admat(admat.chain[[i]], ncol.to.mutate)
        proposed <- plotDAG(admat.star)
        if(FALSE){
            fig <- arrangeGrob(truth, proposed)
            grid.draw(fig)
        }
        fit.star <- calc.tree.fitness(admat.star, cpov, mcf_stats)
        r <- fit.star / fit.prev
        u <- runif(1, 0, 1)
        if(u <= r)
            fit.prev <- fit.star
    }
    fit.prev <- calc.tree.fitness(admat.chain[[i]], cpov, mcf_stats)
    score.chain[i] <- fit.prev
    ##
    ## Propose edge
    ##
    admat.star <- mutate.admat(admat.chain[[i]], ncol.to.mutate)
    proposed <- plotDAG(admat.star)
    fit.star <- calc.tree.fitness(admat.star, cpov, mcf_stats)
    r <- fit.star / fit.prev
    u <- runif(1, 0, 1)    
    if(u <= r) {
        admat.chain[[i+1]] <- admat.star
        numAccept <- numAccept + 1
        accept[i] <- 1
        probs[i] <- fit.star
    } else {
        admat.chain[[i+1]] <- admat.chain[[i]]
        probs[i] <- fit.prev
    }
}
results <- list(admat.chain=admat.chain,
                probs=probs,
                numAccept=numAccept)
saveRDS(results, file.path("..", "output", "mh_trees.R",
                           "trees.rds"))
q('no')
library(ggmcmc)
results <- readRDS(file.path("..", "output", "mh_trees.R",
                             "trees.rds"))
probs <- as.mcmc(tibble(tree.lprob=log(results$probs)))
p <- mcmc.list(probs) %>%
    ggs()
ggs_traceplot(p) +
    scale_y_continuous(expand=expand_scale(mult=0.2)) +
    theme(panel.background=element_rect(fill="white",
                                        color="black"),
          strip.text=element_blank()) +
    ylab("log(p)") 
ggsave(file.path("..", "output", "mh_trees.R",
                 "tree_fitness.pdf"), width=10, height=4)

if(FALSE){
    ## From Thursday
    results <- readRDS(file.path("..", "output", "mh_trees.R",
                                 "tree_fitness.pdf"))
    admat.chain <- results[["admat.chain"]]
    chains <- sapply(admat.chain, function(x){
        x[is.na(x)] <- 0
        x <- as.numeric(x)
        paste(x, collapse="")
    })
    tab <- sort(table(chains), decreasing=TRUE)
    indices <- which(chains == names(tab)[1])
    ch <- admat.chain[[56]]
    plotDAG(ch)
}

