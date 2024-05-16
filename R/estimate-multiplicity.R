get95CI <- function(y, n) {
  test <- prop.test(y, n, correct = T)
  ci_lower <- test$conf.int[1]
  ci_upper <- test$conf.int[2]
  return(c(ci_lower, ci_upper))
}

calcmC <- function(VAF, purity, tcn) {
  mC <- (VAF * (purity * tcn + (1-purity) * 2)) / purity
  return(mC)
}

getmCCI <- function(y, n, purity, tcn) {
  VAF_CI <- get95CI(y, n)
  lower_mC <- calcmC(VAF_CI[1], purity, tcn)
  upper_mC <- calcmC(VAF_CI[2], purity, tcn)
  return(c(lower_mC, upper_mC))
}

assignMultiplicity <- function(lower_mC, upper_mC, tcn) {
  # (1) If the CI for mC overlaps an integer value, 
  # that value is estimated to indicate the multiplicity 
  # of the mutation and the mutation is clonal (C=1)
  if (ceiling(lower_mC) == floor(upper_mC)) {
    test_m <- ceiling(lower_mC)
    # cap m at tcn
    if (test_m > tcn) return (tcn)
    return(test_m)
  }
  
  # (2) If the upper bound of the CI for mC is below 1, 
  # the multiplicity is set to 1, and the mutation is 
  # subclonal, unless the resulting estimate for C is 
  # within a tolerance threshold (0.25) of 1
  if (upper_mC < 1) return(1)
  
  # (3) If the CI for mC is above 1 and does not overlap 
  # any integer values, multiplicity is greater than 1 
  # and m is set such that the confidence interval for C 
  # falls within the expected intervals of [0,1]
  if (lower_mC > 1) {
    keep_going <- TRUE
    m <- 1
    while(keep_going) {
      # cap m at tcn 
      if (m == tcn) return(m)
      
      m <- m + 1
      C_lower <- lower_mC / m
      C_upper <- upper_mC / m
      if (C_lower > 0 & C_upper < 1) {
        keep_going <- FALSE
      }
    }
    return(m)
  }
  
  # CI for mC spans more than 1 integer value including 1
  if (lower_mC <= 1 && upper_mC >= 1) return(1)
  
  return(NA)
}

estimateMultiplicity1 <- function(y, n, tcn) {
  # if no variant reads, assigning multiplicity of 1
  if (y == 0) return(1)
  mC_CI <- getmCCI(y, n, 0.9, tcn)
  m <- max(1, assignMultiplicity(mC_CI[1], mC_CI[2], tcn))
  return(m)
}

#' @export
estimateMultiplicityMatrix <- function(data) {
  Y = data$y
  N = data$n
  Tcn = data$tcn
  I = data$I
  S = data$S
  M = matrix(100, I, S)
  
  for (i in 1:I) {
    for (s in 1:S) {
      y = Y[[i, s]]
      n = N[[i, s]]
      tcn = Tcn[[i, s]]
      m = estimateMultiplicity1(y, n, tcn)
      M[i, s] = m
    }
  }
  return(M)
}