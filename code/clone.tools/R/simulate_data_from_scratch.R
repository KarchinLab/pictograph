samplePurityFromUnif <- function(min.purity, max.purity, S) {
  round(runif(S, min=min.purity, max=max.purity), 2)
}

generateRandomCCFsFromGraph <- function(rand.am, S, K, purity) {
  edges <- rand.am %>%
    filter(connected == 1)
  w <- matrix(NA, nrow=K, ncol=S)
  
  # Assign CCFs for root children
  root_children <- edges$child[edges$parent == "root"]
  # if root has only one child, child CCFs should be equal to 1
  if (length(root_children) == 1) {
    w[as.numeric(root_children), ] <- rep(1, S)
  } else {
    w[as.numeric(root_children), ] <- sampleCCF(length(root_children), S, rep(1, S))
  }
  
  # move down tree
  all_parents <- unique(as.character(edges$parent[which(edges$parent != "root")]))
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
    children.CCF <- round(x * parentCCF, 2)
  } else {
    x <- t(MCMCpack::rdirichlet(numSamples, rep(1, numClusters)))
    y <- matrix(rep(parentCCF, numClusters), nrow=numClusters, byrow = T)
    children.CCF <- round(x * y, 2)
  }
  
  children.CCF
}

.noAllZeroes <- function(w) {
  !any(rowSums(w) == 0)
}
.allDistinct <- function(w) {
  w.distinct <- w %>% 
    as_tibble() %>% 
    distinct()
  nrow(w) == nrow(w.distinct)
}

.isValidW <- function(w) {
  .noAllZeroes(w) & .allDistinct(w)
}

simDataFromScratch <- function(S, K, varPerClustMode, 
                               minVarPerClust, maxVarPerClust, 
                               varPerClust, avg.cov=100,
                               max.num.root.children) {
  purity <- samplePurityFromUnif(min.purity = 0.5, max.purity = 0.9, S)
  rand.am <- generateRandomGraphFromK(K, max.num.root.children)
  w <- generateRandomCCFsFromGraph(rand.am, S, K, purity)
  
  # make sure w is valid (all clusters have distinct CCFs and no all 0's)
  while(!.isValidW(w)) w <- generateRandomCCFsFromGraph(rand.am, S, K, purity)
  
  if (varPerClustMode == "variable") {
    z <- unlist(mapply(function(clust, numVar) rep(clust, numVar), 
                       1:K, 
                       sample(minVarPerClust:maxVarPerClust, K, replace = T)))
  } else if (varPerClustMode == "constant") {
    z <- rep(1:K, each = varPerClust)
  } else stop("mutPerClust must be either variable or constant")
  
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
                    purity=purity,
                    am.long=rand.am)
  test.data
}

inputDataFromSim <- function(sim.data) {
  sim.data[c("y", "n", "purity", "tcn", "m", "I", "S")]
}
