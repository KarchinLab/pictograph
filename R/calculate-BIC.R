#' @import magrittr
calcChainBIC <- function(chains, input.data, pattern) {
  options(dplyr.summarise.inform = FALSE)
  mcf <- chains$mcf_chain%>%
    mutate(value = round(value,5),
           Cluster = as.numeric(gsub("mcf\\[(.*),.*","\\1", Parameter)),
           Sample = as.numeric(gsub(".*,(.*)\\]","\\1", Parameter)))%>%
    group_by(Cluster,Sample)%>%
    reframe(mcf = mean(value))%>%
    ungroup()%>%
    spread(key = Sample, value = mcf)
  
  ww <- writeClusterAssignmentsTable(chains$z_chain)%>%
    mutate(Mut_ID = as.numeric(gsub("Mut","",Mut_ID)))%>%
    arrange(Mut_ID)%>%
    left_join(mcf, by = "Cluster")%>%
    select(-c("Mut_ID","Cluster"))%>%
    as.matrix()
  
  ww <- ww[,which(strsplit(pattern, split="")[[1]]=="1"),drop=FALSE]
  
  mm <- writeMultiplicityTable(chains$m_chain)%>%
    mutate(Mut_ID = as.numeric(gsub("Mut","",Mut_ID)))%>%
    arrange(Mut_ID)%>%
    select(c("Multiplicity")) %>%
    as.matrix()
  
  mm <- replicate(input.data$S, mm[, 1])
  
  is_cn <- replicate(input.data$S, input.data$is_cn)
  vaf <- ifelse(is_cn==0, (ww + (mm-1) * input.data$cncf)/input.data$tcn, (ww * mm + 1 - ww) / input.data$tcn)
  vaf <- ifelse(vaf<=0, 0.001, vaf)
  vaf <- ifelse(vaf>=1, 0.999, vaf)
  
  lik <- sum(dbinom(input.data$y,input.data$n,vaf,log = T))
  
  est_K <- estimateMCFs(chains$mcf_chain) %>% nrow(.)
  
  BIC <- log(input.data$I*input.data$S)*est_K-2*lik
  return(BIC)
}

#' @import magrittr
calcChainSilhouette <- function(chains, input.data, pattern) {
  
  vaf <- input.data$y / input.data$n
  
  zz <- writeClusterAssignmentsTable(chains$z_chain) %>%
    mutate(Mut_ID = as.numeric(gsub("Mut","",Mut_ID))) %>%
    arrange(Mut_ID)%>%
    select(c("Cluster")) %>%
    as.matrix()
  zz <- zz[,1]
  
  mcf1 <- chains$mcf_chain%>%
    mutate(value = round(value,5),
           Cluster = as.numeric(gsub("mcf\\[(.*),.*","\\1", Parameter)),
           Sample = as.numeric(gsub(".*,(.*)\\]","\\1", Parameter)))%>%
    group_by(Cluster,Sample)%>%
    reframe(mcf = mean(value))%>%
    ungroup()%>%
    spread(key = Sample, value = mcf)
  
  ww <- writeClusterAssignmentsTable(chains$z_chain)%>%
    mutate(Mut_ID = as.numeric(gsub("Mut","",Mut_ID)))%>%
    arrange(Mut_ID)%>%
    left_join(mcf1, by = "Cluster")%>%
    select(-c("Mut_ID","Cluster"))%>%
    as.matrix()
  
  ww <- ww[,which(strsplit(pattern, split="")[[1]]=="1"),drop=FALSE]
  
  mm <- writeMultiplicityTable(chains$m_chain)%>%
    mutate(Mut_ID = as.numeric(gsub("Mut","",Mut_ID)))%>%
    arrange(Mut_ID)%>%
    select(c("Multiplicity")) %>%
    as.matrix()
  mm <- replicate(input.data$S, mm[, 1])
  
  is_cn <- replicate(input.data$S, input.data$is_cn)
  
  mcf <- ifelse(is_cn==0, input.data$tcn * vaf - (mm-1) * input.data$cncf, (input.data$tcn * vaf - 1) / (mm - 1))
  mcf <- ifelse(mcf<0, ww, mcf)    
  mcf <- ifelse(mcf>1, ww, mcf)
  
  sil_widths <- silhouette(zz, dist(mcf))
  if (length(sil_widths)==1) {
    mean_silhouette_score <- 0
  } else {
    mean_silhouette_score <- mean(sil_widths[, "sil_width"])
  }
  
  return(mean_silhouette_score)
}

#' Make table listing possible choices of K (minimum BIC and elbow of BIC plot) for each mutation set
#' 
#' @export
#' @import dplyr
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param sample_names (Optional) Vector of sample IDs, same order as provided as input data (e.g. indata$Sample_ID)
writeSetKTable <- function(all_set_results, sample_names = NULL) {
  
min_bic_k <- sapply(all_set_results, function(x) x$BIC_best_K)
    elbow_k <- sapply(all_set_results, function(x) ifelse(is.logical(x$BIC), 1, findElbow(x$BIC$BIC)))
    knee_k <- sapply(all_set_results, function(x) ifelse(is.logical(x$BIC), 1, findKnee(x$BIC$BIC)))
    min_bic_k_tb <- tibble(set_name_bin = names(all_set_results),
                           min_BIC = min_bic_k,
                           elbow = elbow_k,
                           knee = knee_k)
    
    if (!is.null(sample_names)) {
      min_bic_k_tb <- min_bic_k_tb %>%
        mutate(set_name_full = sapply(min_bic_k_tb$set_name_bin, 
                                      function(x) getSetName(x, sample_names, collapse_string = ","))) %>%
        select(set_name_bin, set_name_full, min_BIC, elbow, knee)
    }
    
    min_bic_k_tb <- min_bic_k_tb %>%
      mutate(BIC_K = round((min_BIC + elbow + knee) / 3))
    
    silhouette_k <- sapply(all_set_results, function(x) x$silhouette_best_K)
    min_bic_k_tb <- min_bic_k_tb %>%
      mutate(silhouette_K = silhouette_k)

  return(min_bic_k_tb)
}

findElbow <- function(BIC) {
  
  # return first ind if all increasing 
  if(all (BIC >= BIC[1])) return(1)
  
  K_vec <- seq_len(length(BIC))
  # line defined by two end-points of BIC plot P1 = (x1, y1) and P2 = (x2, y2)
  x1 <- K_vec[1]
  y1 <- BIC[1]
  x2 <- length(BIC)
  y2 <- BIC[length(BIC)]
  
  # calculate distance of each point (x0, y0) to line
  perp_dist <- sapply(K_vec, function(i) calcDistFromPointToLine(K_vec[i], BIC[i],
                                                                 x1, y1, 
                                                                 x2, y2))
  return(which.max(perp_dist))
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

calcDistFromPointToLine <- function (x0, y0,
                                     x1, y1,
                                     x2, y2) {
  numerator <- abs( (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1) )
  denominator <- sqrt( (x2-x1)^2 + (y2-y1)^2 )
  distance <- numerator / denominator
  return(distance)
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
