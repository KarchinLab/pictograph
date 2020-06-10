samplePurityFromUnif <- function(min.purity, max.purity, S) {
  round(runif(S, min=min.purity, max=max.purity), 2)
}

generateRandomCCFsFromGraph <- function(rand.am, S, K) {
  edges <- rand.am[rand.am$connected == 1, ]
  w <- matrix(NA, nrow=K, ncol=S)
  
  # Assign CCFs for root children
  root_children <- edges$child[edges$parent == "root"]
  w[as.numeric(root_children), ] <- sampleCCF(length(root_children), S, rep(1, S))
  
  # move down tree
  all_parents <- as.character(edges$parent[which(edges$parent != "root")])
  next_parents <- all_parents[all_parents %in% root_children]
  
  # assign rest of CCFs
  while (length(next_parents) > 0) {
    temp_parent <- next_parents[1]
    temp_children <- edges[edges$parent == temp_parent, ]$child
    w[as.numeric(temp_children), ] <- sampleCCF(length(temp_children), S, w[as.numeric(temp_parent), ])
    
    next_parents <- c(next_parents[-1], all_parents[all_parents %in% temp_children])
  }
  
  w
}

sampleCCF <- function(numClusters, numSamples, parentCCF) {
  if (numClusters == 1) {
    x <- runif(numSamples, 0, 1)
  } else {
    x <- t(MCMCpack::rdirichlet(numSamples, rep(1, numClusters)))
  }
  x * parentCCF
}

simDataFromScratch <- function(S, K, minVarPerClust, maxVarPerClust, avg.cov=100) {
  purity <- samplePurityFromUnif(min.purity = 0.5, max.purity = 0.9, S)
  rand.am <- generateRandomGraphFromK(K)
  w <- round(generateRandomCCFsFromGraph(rand.am, S, K), 2)
  
  z <- unlist(mapply(function(clust, numVar) rep(clust, numVar), 
                     1:K, 
                     sample(minVarPerClust:maxVarPerClust, K)))
  I <- length(z)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  
  # variant x sample matrices
  W <- w[z, ]
  P <- matrix(rep(purity, each = I), nrow = I, ncol = S)
  
  theta <- calcTheta2(m, tcn, W, P)
  ##
  ## Simulate altered reads
  ##
  n <- replicate(S, rpois(I, avg.cov))
  y <- matrix(NA, nrow=I, ncol=S)
  for (i in 1:I) {
    for (s in 1:S) {
      y[i, s] <- rbinom(1, n[i, s], theta[i,s])
    }
  }
  
  colnames(w) <- colnames(m) <- paste0("sample", 1:S)
  test.data <- list("I" = I, "S" = S, "K" = K, 
                    "y" = y, "n" = n,
                    "m" = m, "tcn" = tcn,
                    z=z,
                    w=w,
                    theta=theta,
                    purity=purity)
  test.data
}
