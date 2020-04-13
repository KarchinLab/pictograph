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
  
  theta <- calcTheta(m, tcn, W)
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