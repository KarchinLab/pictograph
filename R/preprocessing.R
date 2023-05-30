#' read input data file and store in PICTOGRAPH input format
#' 
#' @export
#' @param input_file input data file; 
importCSV <- function(inputFile) {
  data <- read_csv(inputFile, show_col_types = FALSE)
  output_data <- list()
  
  output_data$y <- as.matrix(data[c("mutation", "sample", "alt_reads")] %>% pivot_wider(names_from = sample, values_from = alt_reads, values_fill = 0))
  rownames(output_data$y) <- output_data$y[,'mutation']
  output_data$y <- output_data$y[,-1, drop=FALSE]
  rowname = rownames(output_data$y)
  colname = colnames(output_data$y)
  output_data$y <- matrix(as.numeric(output_data$y), ncol = ncol(output_data$y))
  rownames(output_data$y) = rowname
  colnames(output_data$y) = colname
  
  output_data$MutID <- rowname
  
  output_data$n <- as.matrix(data[c("mutation", "sample", "total_reads")] %>% pivot_wider(names_from = sample, values_from = total_reads, values_fill = 0))
  rownames(output_data$n) <- output_data$n[,'mutation']
  output_data$n <- output_data$n[,-1, drop=FALSE]
  rowname = rownames(output_data$n)
  colname = colnames(output_data$n)
  output_data$n <- matrix(as.numeric(output_data$n), ncol = ncol(output_data$n))
  rownames(output_data$n) = rowname
  colnames(output_data$n) = colname
  
  if (any((output_data$y - output_data$n) > 0)) {
    warning("Total read count must be equal or bigger than alt read count. Please check input data before proceeding!")
    stop()
  }
  
  if (any(output_data$n==0)) {
    warning("Total read counts of 0 encoutered. Replaced 0 with mean total read count.")
    output_data$n[output_data$n==0] <- round(mean(output_data$n))
  }
  
  output_data$tcn <- as.matrix(data[c("mutation", "sample", "tumor_integer_copy_number")] %>% pivot_wider(names_from = sample, values_from = tumor_integer_copy_number, values_fill = 0))
  rownames(output_data$tcn) <- output_data$tcn[,'mutation']
  output_data$tcn <- output_data$tcn[,-1, drop=FALSE]
  rowname = rownames(output_data$tcn)
  colname = colnames(output_data$tcn)
  output_data$tcn <- matrix(as.numeric(output_data$tcn), ncol = ncol(output_data$tcn))
  rownames(output_data$tcn) = rowname
  colnames(output_data$tcn) = colname
  
  output_data$purity <- as.matrix(data[c("mutation", "sample", "purity")] %>% pivot_wider(names_from = sample, values_from = purity))
  rownames(output_data$purity) <- output_data$purity[,'mutation']
  output_data$purity <- output_data$purity[,-1, drop=FALSE]
  rowname = rownames(output_data$purity)
  colname = colnames(output_data$purity)
  output_data$purity <- matrix(as.numeric(output_data$purity), ncol = ncol(output_data$purity))
  rownames(output_data$purity) = rowname
  colnames(output_data$purity) = colname
  if (ncol(output_data$purity) == 1) {
    if (length(unique(output_data$purity[,])) != 1) {
      warning("purity not consistent for the same sample; taking the mean purity")
      output_data$purity <- round(colMeans(output_data$purity), digit=2)
    } else {
      output_data$purity <- unique(output_data$purity[,])
    }
  } else {
    if (nrow(unique(output_data$purity[,])) != 1) {
      warning("Purity not consistent for the same sample; taking the mean purity")
      output_data$purity <- round(colMeans(output_data$purity), digit=2)
    } else {
      output_data$purity <- unique(output_data$purity[,])[1,]
    }
  }
  
  output_data$S = ncol(output_data$y)
  output_data$I = nrow(output_data$y)
  
  if ("multiplicity" %in% colnames(data)) {
    output_data$m <- as.matrix(data[c("mutation", "sample", "multiplicity")] %>% pivot_wider(names_from = sample, values_from = multiplicity))
    rownames(output_data$m) <- output_data$m[,'mutation']
    output_data$m <- output_data$m[,-1, drop=FALSE]
    rowname = rownames(output_data$m)
    colname = colnames(output_data$m)
    output_data$m <- matrix(as.numeric(output_data$m), ncol = ncol(output_data$m))
    rownames(output_data$m) = rowname
    colnames(output_data$m) = colname
    if(any(data$multiplicity == 0)){
      warning("Multiplicity of 0 encoutered in the input data file which may cause issue for jags model!")
    }
  } else {
    output_data$m <- estimateMultiplicityMatrix(output_data)
  }
  
  return(output_data)
}

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

estimateMultiplicity <- function(y, n, purity, tcn) {
  # if no variant reads, assigning multiplicity of 1
  if (y == 0) return(1)
  mC_CI <- getmCCI(y, n, purity, tcn)
  m <- max(1, assignMultiplicity(mC_CI[1], mC_CI[2], tcn))
  return(m)
}

estimateMultiplicityMatrix <- function(data) {
  Y = data$y
  N = data$n
  Tcn = data$tcn
  Purity = data$purity
  I = data$I
  S = data$S
  M = matrix(100, I, S)
  
  for (i in 1:I) {
    for (s in 1:S) {
      y = Y[[i, s]]
      n = N[[i, s]]
      tcn = Tcn[[i, s]]
      purity = Purity[s]
      m = estimateMultiplicity(y, n, purity, tcn)
      M[i, s] = m
    }
  }
  return(M)
}