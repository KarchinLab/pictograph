calcLogLik <- function(z.iter, w.iter, input.data) {
  W <- w.iter[z.iter, ]
  
  if (is.null(input.data$purity)) {
    theta <- calcTheta(input.data$m, input.data$tcn, W)
  } else {
    purity <- input.data$purity
    P <- matrix(rep(purity, each = input.data$I), nrow = input.data$I, ncol = input.data$S)
    theta <- calcTheta2(input.data$m, input.data$tcn, W, P)
  }
  
  sum(dbinom(as.matrix(input.data$y), as.matrix(input.data$n), as.matrix(theta), log=T))
}

calcChainLogLik <- function(chains, input.data, est_K) {
  num_iter <- max(chains$z_chain$Iteration)
  
  lik <- c()
  for(iter in 1:num_iter) {
    z.iter <- chains$z_chain %>%
      filter(Iteration == iter) %>%
      pull(value)
    w.iter <- filter(chains$w_chain, Iteration == iter) %>% 
      reshapeW()
    lik <- c(lik, calcLogLik(z.iter, w.iter, input.data))
  }
  return(mean(lik))
}


#' @importFrom stringr str_replace_all
reshapeW <- function(w.chain.iter) {
  w.mat <- w.chain.iter %>%
    mutate(sample=stringr::str_replace_all(Parameter, "w\\[[:digit:]+,", ""),
           sample=as.numeric(stringr::str_replace_all(sample, "\\]", "")),
           cluster=stringr::str_replace_all(Parameter, "w\\[", ""),
           cluster=as.numeric(stringr::str_replace_all(cluster, ",[:digit:]\\]", ""))) %>%
    select(cluster, sample, value) 
  S <- max(w.mat$sample)
  w.mat <- w.mat %>%
    pivot_wider(names_from = sample, 
                values_from = value)
  w.mat$cluster <- NULL
  w.mat <- as.matrix(w.mat)
  colnames(w.mat) <- paste0("sample", 1:S)
  return(w.mat)
}

calcBIC <- function(n, k, ll) log(n)*k - 2*ll

#' @import magrittr
calcChainBIC <- function(chains, input.data) {
  n <- input.data$I * input.data$S
  est_K <- estimateCCFs(chains$w_chain) %>%
    nrow(.)
  ll <- calcChainLogLik(chains, input.data, est_K)
  
  BIC <- calcBIC(n, est_K, ll)
  return(BIC)
}

calcBICForRangeK <- function(samps.list, kToTest, input.data) {
  mapply(function(samps, k)
    calcBIC(input.data$I*input.data$S, k, calcChainLogLik(samps, input.data, k)),
    samps = samps.list, k = kToTest)
}

calcTheta <- function(m, tcn, w) {
  (m * w) / (tcn * w + 2*(1-w))
}

calcTheta2 <- function(m, tcn, w, p) {
  (m * w * p) / (tcn * p + 2*(1-p))
}

#' Make table listing possible choices of K (minimum BIC and elbow of BIC plot) for each mutation set
#' 
#' @export
#' @import dplyr
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param sample_names (Optional) Vector of sample IDs, same order as provided as input data (e.g. indata$Sample_ID)
writeSetKTable <- function(all_set_results, sample_names = NULL) {
  min_bic_k <- sapply(all_set_results, function(x) x$best_K)
  elbow_k <- sapply(all_set_results, function(x) ifelse(is.logical(x$BIC), 1, findElbow1(x$BIC$BIC)))
  knee_k <- sapply(all_set_results, function(x) ifelse(is.logical(x$BIC), 1, findKnee(x$BIC$BIC)))
  min_bic_k_tb <- tibble(set_name_bin = names(all_set_results),
                         min_BIC = min_bic_k,
                         elbow = elbow_k,
                         knee = knee_k)
  
  if (!is.null(sample_names)) {
    min_bic_k_tb <- min_bic_k_tb %>%
      mutate(set_name_full = sapply(min_bic_k_tb$set_name_bin, 
                                    function(x) getSetName(x, sample_names, collapse_string = ","))) %>%
      select(set_name_bin, set_name_full, min_BIC, elbow)
  }
  
  # write chosen K if elbow, knee, and minimum BIC all agree 
  chosen_K <- rep(NA, length(all_set_results))
  for (i in seq_len(length(all_set_results))) {
    chosen_K[i] <- ifelse(min_bic_k_tb$min_BIC[i] == min_bic_k_tb$elbow[i] & 
                            min_bic_k_tb$min_BIC[i] == min_bic_k_tb$knee[i],
                          min_bic_k_tb$min_BIC[i],
                          NA)
  }
  
  min_bic_k_tb <- min_bic_k_tb %>%
    mutate(chosen_K = chosen_K)
  
  return(min_bic_k_tb)
}

# find elbow of bic plot
findElbow1 <- function(BIC) {
  
  # return first ind if all increasing 
  if(all (BIC >= BIC[1])) return(1)
  
  delta1 <- diff(BIC)
  delta2 <- diff(delta1)
  elbow_ind <- which.min(delta2) + 2
  return(elbow_ind)
}

# angle-based method for knee point detection of BIC 
# Zhao et al 2008
# returns index of knee point
findKnee <- function(BIC, n = 5) {
  if (length(BIC) == 1) return(1)
  
  # return first ind if all increasing 
  if(all (BIC >= BIC[1])) return(1)
  
  # initialize
  curr_val <- BIC[1]
  prev_val <- BIC[1]
  next_val <- BIC[1]
  
  # begin
  diff_fun <- rep(NA, length(BIC))
  for (m in seq_len(length(BIC))) {
    curr_val <- BIC[m]
    next_val <- BIC[min(m+1, length(BIC))]
    diff_fun[m] <- DiffFun(prev_val, next_val, curr_val)
    prev_val <- curr_val
  }
  
  # find first n local maximas in diff_fun
  local_max <- localMaxima(diff_fun)
  local_max <- local_max[1:min(length(local_max), n)]
  if (length(local_max) == 1) return(local_max)
  
  # for each n with decreasing order of LocalMax value,
  angle <- c()
  for (n in local_max) {
    angle <- c(angle, AngleFun(BIC[max(n-1, 1)], BIC[min(n+1, length(BIC))], BIC[n]))
  }
  # return m with the first minima angle
  return(local_max[localMinima(angle)[1]])
}
DiffFun <- function(prev_val, next_val, curr_val) {
  prev_val + next_val - 2*curr_val
}
AngleFun <- function(prev_val, next_val, curr_val) {
  atan( 1 / abs(curr_val - prev_val) ) + atan( 1 / abs(next_val - curr_val) )
}

# detect local maxima
localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L # for maxima
  #y <- diff(c(.Machine$integer.max, x)) < 0L # for mimuma
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}
localMinima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  #y <- diff(c(-.Machine$integer.max, x)) > 0L # for maxima
  y <- diff(c(.Machine$integer.max, x)) < 0L # for mimuma
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}
