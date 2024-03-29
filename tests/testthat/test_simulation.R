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
    library(stringr)
    extdir <- system.file("extdata", package="clone.tools")
    jags.file <- file.path(extdir, "mcf_model.jag")
    dat <- unit_test_data()
    jags_inputs <- listJagInputs(dat)
    model <- jags.model(jags.file,
                        data=jags_inputs,
                        ## confusing to have more than 1 chain
                        ## because of label-switching
                        n.chains=1, 
                        n.adapt=500)
    samples <- coda.samples(model, n.iter=1000, thin=1,
                            variable.names=c("w", "z", "ystar"))
    tib <- ggs(samples)

    mcf.chain <- ggs(samples, family="w") %>%
        mutate(cluster_id=clusterIndex(Parameter),
               sample_id=sampleIndex(Parameter))
    z.chain <- ggs(samples, family="z") %>%
        mutate(variant_id=mutationIndex(Parameter)) %>%
        rename(cluster_id=value)
    ystar.chain <- ggs(samples, family="ystar") %>%
        mutate(variant_id=clusterIndex(Parameter),
               sample_id=sampleIndex(Parameter))
    expect_true(length(unique(z.chain$variant_id)) == jags_inputs$I)

    ##
    ## plot cluster assignments
    ##
    z.chain %>%
        group_by(variant_id, cluster_id) %>%
        tally() %>%
        mutate(prob=n/1000) %>%
        ggplot(aes(variant_id, cluster_id)) +
        geom_point(aes(size=prob)) +
        theme_bw()
    ##
    ## Summarize MCFs
    ##
    mcf.chain %>%
        group_by(cluster_id, sample_id) %>%
        summarize(mcf=mean(value)) %>%
        ggplot(aes(cluster_id, mcf)) +
        geom_point() +
        theme_bw(base_size=13)

    mcf <- mcf.chain %>%
        group_by(cluster_id, sample_id) %>%
        summarize(mean=mean(value))
    ##trace(constrainedEdges, browser)
    A <- constrainedEdges(mcf, zero.thresh=0.01)
    wmat <- spread(mcf, sample_id, mean) %>%
        ungroup() %>%
        select(-cluster_id) %>%
        as.matrix()    
    A2 <- base.admat(wmat, zero.thresh = 0.01)

    tmp <- A %>%
        select(parent, child, connected) %>%
        spread(child, connected) %>%
        select(-parent) %>%
        as.matrix()
    expect_equivalent(tmp, A2)
})



test_that("initializing a graph", {
    ##
    ## For now, proceed with graph sampler
    ## based on posterior means of the MCFs
    ##
    set.seed(123)
    mcf <- matrix(c(0.98, 0.99, 0.97, 
                    0.55, 0.00, 0.80, 
                    0.30, 0.70, 0.00,
                    0.20, 0.22, 0.18),
                  byrow=TRUE,
                  nrow=4, ncol=3)
    mcf.long <- tibble(cluster_id=as.character(rep(1:4, 3)),
                       sample_id=as.character(rep(1:3, each=4)),
                       mean=as.numeric(mcf))
    set.seed(2)
    am <- init.admat(mcf, zero.thresh=0.01)
    am2.long <- constrainedEdges(mcf.long)
    ##
    ## Below, we establish that the edges in am3.long are equivalent
    ## to am3. This means we could replace the first part of
    ## rand.admat function with the long formatted data
    ##
    set.seed(2)
    am3.long <- randAdmat(am2.long)
    am3.wide <- am3.long %>%
        select(parent, child, connected) %>%
        spread(child, connected)
    tmp <- am3.wide %>% select(-parent) %>%
        as.matrix()
    ## remove row with all NAs in am
    am <- am[rowSums(is.na(am)) < 4, ]
    expect_equivalent(tmp, am)
    expect_true(validGraph(am3.long))

    ## constrainedEdges and randAdmat are combined in the function initializeGraph
    set.seed(2)
    am2 <- initializeGraph(mcf)
    expect_equivalent(am2, am3.long)
    ##
    ## Note, only the root and parents 1 and 4 are allowed to have
    ## children as clusters 2 and 3 have samples with undetectable
    ## mutations
    ##
    ##
    ## The next part of rand.admat requires that there is at
    ## least one edge from the root.
    ##
    expect_true(isRootConnected(am3.long))
    ##
    ## Disconnect root
    ##
    am3.long$connected[[1]] <- 0
    expect_true(!isRootConnected(am3.long))
    expect_true(!validGraph(am3.long))
    ##
    ## connecting the root
    ##
    ##  1. randomly select child to connect to root
    ##     
    ##  2. For the child selected, disconnect any other edges to that child
    ##
    ##  3. Check that we have a valid directed acyclic graph (DAG)
    ##   
    edge <- am3.long %>%
        filter(parent == "root") %>%
        sample_n(1)
    am3.long2 <- addEdge(am3.long, edge)
    expect_true(validGraph(am3.long2))
    expect_true(isRootConnected(am3.long2))
    expect_true(isDirected(am3.long2))
})

test_that("mutateA", {
    set.seed(123)
    mcf <- matrix(c(0.98, 0.99, 0.97, 
                    0.55, 0.00, 0.80, 
                    0.30, 0.70, 0.00,
                    0.20, 0.22, 0.18),
                  byrow=TRUE,
                  nrow=4, ncol=3)
    mcf.long <- tibble(cluster_id=as.character(rep(1:4, 3)),
                       sample_id=as.character(rep(1:3, each=4)),
                       mean=as.numeric(mcf))
    set.seed(2)
    A <- init.admat(mcf, zero.thresh=0.01)

    set.seed(2)
    a <- initializeGraph(mcf)
    expect_equivalent(toWide(a), A[rowSums(is.na(A)) < 4, ])
    expect_true(validGraph(a))
    
    set.seed(3)
    astar <- .mutateA(A)
    K <- ncol(A)
    npossible <- colSums(!is.na(A))
    columns <- which(npossible > 1)
    set.seed(3)
    rand.k <- sample(columns, size=1)

    set.seed(3)
    possible_moves <- a %>% group_by(child) %>%
        summarize(n=sum(connected==0)) %>%
        mutate(possible=which(n >= 1)) %>%
        sample_n(1)
    ## rand.k and child are both 1
    expect_true(rand.k == pull(possible_moves, child))

    ## list the possible edges to child 1
    new_edge <- filter(a, child==pull(possible_moves, child),
                       connected==0)
    a2 <- addEdge(a, new_edge)
    astar2 <- astar[rowSums(is.na(astar))< 4, ]
    toWide(a2)
    expect_equivalent(toWide(a2), astar2)
    ##
    ## The test below fails.
    ## graph is not valid since root is not connected
    ##
    ## TODO: I think we should evaluate the possible valid move sets
    ## and never move to an invalid graph.
    expect_true(validGraph(a2))
    ## 
    ## assumptions:
    ##  The current graph is valid, which means the following
    ##    - is directed
    ##    - each child has a single 1
    ##    - the root is connected to at least 1 child
    ##
    ## condition 1:
    ##     child node must have 1 or more zeros
    condition1 <- a %>% group_by(child) %>%
        summarize(n=sum(connected==0)) %>%
        mutate(possible=which(n >= 1))
    ## condition 2:
    ##   for each child and for each zero of that child
    ##     determine whether changing the 1 to a zero and the zero to a 1 would be valid
    possible_moves <- filter(a, connected==0, child %in% condition1$child)
    is_valid <- rep(NA, nrow(possible_moves))
    for(i in seq_len(nrow(possible_moves))){
        astar <- addEdge(a, possible_moves[i, ])
        is_valid[i] <- validGraph(astar)
    }
    move_set <- possible_moves[is_valid, ]
    ix <- sample(seq_len(nrow(move_set)), 1)
    astar <- addEdge(a, move_set[ix, ])
    expect_true(validGraph(astar))

    ## sampleNewEdge is wrapper for above
    astar2 <- sampleNewEdge(a)
    expect_true(validGraph(astar2))
})



test_that("skip", {
    skip(" misc functions ")
    addRootEdge <- function(A){
        K <- ncol(A)
        j <- sample(seq_len(K), 1)
        ind.1 <- which(A[, j] == 1)
        A[1, j] <- 1
        A[ind.1, j] <- 0
        A
    }

    isFullyConnected <- function(A){
        numClusters <- nrow(A)
        nodesInMainTree <- bfs(A)
        numNodesInMainTree <- length(nodesInMainTree)
        numClusters == numNodesInMainTree
    }


    ## refactored mutateA
    edgeSampler <- function(A) {
        ## choose a column to mutate
        K <- ncol(A)
        npossible <- colSums(!is.na(A))
        columns <- which(npossible > 1)
        rand.k <- sample(columns, size=1)
        possible_edges <- which(!is.na(A[, rand.k]) & A[, rand.k] != 1)
        ## current position with 1
        ind.1 <- which(A[, rand.k] == 1)
        ## select new position
        if (length(possible_edges) == 1) {
            new.1 <- possible_edges
        } else {
            new.1 <- sample(possible_edges, size=1)
        }
        A2 <- A
        if(length(ind.1) > 0)
            A2[ind.1, rand.k] <- 0
        A2[new.1, rand.k] <- 1
        if(sum(A2[1, ]) == 0) {
            A2 <- addRootEdge(A2)
        }
        if(!isFullyConnected(A2)) return(A)
        A2
    }

    sampleEdges <- function(A){
        A2 <- edgeSampler(A)
        while(!isFullyConnected(A2)){
            j <- sample(seq_len(ncol(A)), 1)
            A2 <- mutate.column(A2, j)
        }
        ## fix bi-directional
        ##
        A2
    }


##    trace(sampleEdges, browser)
    Anext <- sampleEdges(A)
})

