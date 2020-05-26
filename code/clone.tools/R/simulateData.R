simulateData <- function(I, K, S, avg.cov=100){
  pi <- rep(1/K, K)
  ##z <- sample(1:K, size = I, replace = T, prob = pi)
  ## True cancer cell fraction
  z <- rep(1:K, each=I/K)
  w <- matrix(c(0.98, 0.99, 0.97, 
                0.98, 0.90, 0.82,
                0.55, 0.00, 0.80, 
                0.20, 0.00, 0.50,
                0.30, 0.00, 0.30, 
                0.43, 0.90, 0.00,
                0.30, 0.70, 0.00,
                0.20, 0.00, 0.00,
                0.00, 0.00, 0.30,
                0.00, 0.50, 0.00),
              byrow=T,
              nrow=K, ncol=S)
  
  colnames(w) <-  paste0("sample", 1:S)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  W <- w[z, ]
  
  theta <- calcTheta(m, tcn, W)
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
  
  p <- y/n
  colnames(w) <- colnames(p) <- paste0("sample", 1:S)
  colnames(m) <- colnames(p)
  test.data <- list("I" = I, "S" = S, "K" = K, 
                    "y" = y, "n" = n,
                    "m" = m, "tcn" = tcn,
                    z=z,
                    w=w,
                    theta=theta,
                    p=p)
  test.data
}

simulateData3 <- function(I, K, S, avg.cov=100){
  pi <- rep(1/K, K)
  ##z <- sample(1:K, size = I, replace = T, prob = pi)
  ## True cancer cell fraction
  z <- rep(1:K, each=I/K)
  w <- matrix(c(0.98, 0.99, 0.97, 
                0.98, 0.90, 0.82,
                0.55, 0.00, 0.80, 
                0.20, 0.00, 0.50,
                0.30, 0.00, 0.30, 
                0.43, 0.90, 0.00,
                0.30, 0.70, 0.00,
                0.06, 0.00, 0.00,
                0.00, 0.00, 0.08,
                0.00, 0.07, 0.00),
              byrow=T,
              nrow=K, ncol=S)
  
  colnames(w) <-  paste0("sample", 1:S)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  W <- w[z, ]
  
  theta <- calcTheta(m, tcn, W)
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
  
  p <- y/n
  colnames(w) <- colnames(p) <- paste0("sample", 1:S)
  colnames(m) <- colnames(p)
  test.data <- list("I" = I, "S" = S, "K" = K, 
                    "y" = y, "n" = n,
                    "m" = m, "tcn" = tcn,
                    z=z,
                    w=w,
                    theta=theta,
                    p=p)
  test.data
}

simulateData2 <- function(I, K, S){
  pi <- rep(0.1, 10)
  ##z <- sample(1:K, size = I, replace = T, prob = pi)
  ## True cancer cell fraction
  z <- rep(1:10, each=10)
  w <- matrix(c(0.98, 0.99, 0.97, 
                0.98, 0.90, 0.82,
                0.55, 0.00, 0.80, 
                0.20, 0.00, 0.50,
                0.30, 0.00, 0.30, 
                0.43, 0.90, 0.00,
                0.30, 0.70, 0.00,
                0.20, 0.00, 0.00,
                0.00, 0.00, 0.30,
                0.00, 0.50, 0.00),
              byrow=T,
              nrow=K, ncol=S)
  
  colnames(w) <-  paste0("sample", 1:S)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  W <- w[z, ]
  m[w==0] <- 0
  
  theta <- calcTheta2(m, tcn, W, )
  ##
  ## Simulate altered reads
  ##
  n <- replicate(S, rpois(I, 100))
  y <- matrix(NA, nrow=I, ncol=S)
  for (i in 1:I) {
    for (s in 1:S) {
      y[i, s] <- rbinom(1, n[i, s], theta[i,s])
    }
  }
  
  p <- y/n
  colnames(w) <- colnames(p) <- paste0("sample", 1:S)
  colnames(m) <- colnames(p)
  test.data <- list("I" = I, "S" = S, "K" = K, 
                    "y" = y, "n" = n,
                    "m" = m, "tcn" = tcn,
                    z=z,
                    w=w,
                    theta=theta,
                    p=p)
  test.data
}

simulateDataPurity <- function(I=120, K=10, S=3, avg.cov=100,
                               purity=c(0.80, 0.85, 0.90)){
  set.seed(1234)
  pi <- rep(1/K, K)
  z <- rep(1:K, each=I/K)
  ## True cancer cell fraction
  w <- matrix(c(1, 1, 1, 
                0.98, 0.90, 0.82,
                0.55, 0.00, 0.80, 
                0.20, 0.00, 0.50,
                0.30, 0.00, 0.30, 
                0.43, 0.90, 0.00,
                0.30, 0.70, 0.00,
                0.13, 0.00, 0.00,
                0.00, 0.00, 0.16,
                0.00, 0.14, 0.00),
              byrow=T,
              nrow=K, ncol=S)
  
  colnames(w) <-  paste0("sample", 1:S)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  
  # variant x sample matrix
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
  
  colnames(m)  <- colnames(w) <- colnames(P) <- paste0("sample", 1:S)
  
  test.data <- list("I" = I, "S" = S, "K" = K, 
                    "y" = y, "n" = n,
                    "m" = m, "tcn" = tcn,
                    z=z,
                    w=w,
                    theta=theta,
                    purity=purity)
  test.data
}

simTestCase1 <- function(){
  # single sample clusters with few variants
  set.seed(1234)
  K=10; S=3;
  purity=c(0.80, 0.85, 0.90)
  avg.cov=100
  
  z <- c(rep(1:7, each=20), 
         8, 9, 9, 10, 10, 10)
  I <- length(z)
  ## True cancer cell fraction
  w <- matrix(c(1.00, 1.00, 1.00,
                0.98, 0.90, 0.82, 
                0.55, 0.00, 0.80, 
                0.10, 0.00, 0.20, 
                0.43, 0.90, 0.00,
                0.30, 0.70, 0.00,
                0.30, 0.10, 0.00,
                0.10, 0.00, 0.00,
                0.00, 0.00, 0.10,
                0.00, 0.07, 0.00),
              byrow=T,
              nrow=K, ncol=S)
  
  colnames(w) <-  paste0("sample", 1:S)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  
  # variant x sample matrix
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

simTestCaseIP30 <- function(){
  # single sample clusters with few variants
  set.seed(1234)
  K=19; S=3;
  purity=c(0.80, 0.85, 0.90)
  avg.cov=100
  
  z <- rep(1:19, each=5)
  I <- length(z)
  ## True cancer cell fraction
  w <- matrix(c(0.53, 0.61, 0.00, 0.00, 0.00,
                0.00, 0.17, 0.00, 0.00, 0.00,
                0.48, 0.47, 0.17, 0.00, 0.00,
                0.98, 0.99, 0.00, 0.00, 0.00,
                0.65, 0.58, 0.97, 0.98, 0.00,
                0.19, 0.00, 0.00, 0.00, 0.00,
                0.99, 0.94, 0.98, 0.98, 0.00,
                1.00, 0.99, 0.23, 0.00, 0.00,
                0.00, 0.00, 0.99, 0.58, 0.00,
                0.00, 0.00, 0.98, 0.99, 0.00,
                0.97, 0.99, 0.92, 0.63, 0.63,
                0.99, 0.59, 0.25, 0.00, 0.00,
                0.63, 0.97, 0.95, 0.98, 0.26,
                1.00, 0.99, 1.00, 0.97, 0.72,
                0.49, 0.24, 0.00, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.14, 0.00,
                0.00, 0.00, 0.18, 0.00, 0.00,
                0.00, 0.00, 0.00, 0.00, 0.40,
                0.00, 0.00, 0.33, 0.43, 0.00),
              byrow=T,
              nrow=K, ncol=S)
  
  colnames(w) <-  paste0("sample", 1:S)
  
  tcn <- matrix(2, nrow=I, ncol=S)
  m <- matrix(rep(sample(1:2, size = I, replace = T), S), 
              nrow=I, ncol=S)
  
  # variant x sample matrix
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