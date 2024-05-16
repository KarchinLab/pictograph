#' run MCMC using JAGS
#' 
#' @export
runMCMCForAllBoxes <- function(sep_list,
                               sample_presence=TRUE,
                               ploidy=2,
                               max_K = 5,
                               min_mutation_per_cluster = 5,
                               cluster_diff_thresh=0.05,
                               n.iter = 5000,
                               n.burn = 1000, 
                               thin = 10,
                               mc.cores = 4, 
                               inits = list(".RNG.name" = "base::Wichmann-Hill",
                                            ".RNG.seed" = 123)){
  
  if (!sample_presence) {
    all_set_results <- vector("list", 1)
    names(all_set_results) <- paste0(rep("1", ncol(sep_list$y)), collapse = "")
    params = c("z", "mcf", "icn", "m", "ystar")
    
    temp_box <- sep_list
    temp_box$pattern <- paste0(rep("1", ncol(sep_list$y)), collapse = "")
    temp_max_K <- min(max_K, floor(nrow(temp_box$y)/min_mutation_per_cluster))
    temp_max_K <- max(temp_max_K, 1)
    
    temp_samps_list <- runMutSetMCMC(temp_box, 
                                     ploidy=ploidy,
                                     n.iter = n.iter, 
                                     n.burn = n.burn, 
                                     thin = thin, 
                                     mc.cores = mc.cores,
                                     inits = inits,
                                     temp_max_K = temp_max_K,
                                     params = params,
                                     min_mutation_per_cluster = min_mutation_per_cluster, 
                                     cluster_diff_thresh = cluster_diff_thresh,
                                     sample_presence=sample_presence)
    
    all_set_results[[1]] <- temp_samps_list
    
  } else {
    all_set_results <- vector("list", length(sep_list))
    names(all_set_results) <- names(sep_list)
    params = c("z", "mcf", "icn", "m", "ystar")
    
    for (i in seq_len(length(sep_list))) {
      temp_box <- sep_list[[i]]
      # Max number of clusters cannot be more than number of mutations/min_mutation_per_cluster
      temp_max_K <- min(max_K, floor(length(temp_box$mutation_indices)/min_mutation_per_cluster))
      temp_max_K <- max(temp_max_K, 1)
      temp_samps_list <- runMutSetMCMC(temp_box,
                                       ploidy=ploidy,
                                       n.iter = n.iter, 
                                       n.burn = n.burn, 
                                       thin = thin, 
                                       mc.cores = mc.cores,
                                       inits = inits,
                                       temp_max_K = temp_max_K,
                                       params = params,
                                       min_mutation_per_cluster = min_mutation_per_cluster, 
                                       cluster_diff_thresh = cluster_diff_thresh,
                                       sample_presence=sample_presence)
      
      all_set_results[[i]] <- temp_samps_list
    }
  }
  return(all_set_results)
}

runMutSetMCMC <- function(temp_box, 
                          ploidy=2,
                          n.iter = 10000, 
                          n.burn = 1000, 
                          thin = 10, 
                          mc.cores = 1,
                          inits = list(".RNG.name" = "base::Wichmann-Hill",
                                       ".RNG.seed" = 123),
                          temp_max_K = 5,
                          params = c("z", "mcf", "icn", "m", "ystar"),
                          min_mutation_per_cluster = 1,
                          cluster_diff_thresh=0.05,
                          sample_presence=FALSE) {

  temp_samps_list <- runMCMCForABox(temp_box,
                                    ploidy=ploidy,
                                    n.iter = n.iter, 
                                    n.burn = n.burn, 
                                    thin = thin, 
                                    mc.cores = mc.cores,
                                    inits = inits,
                                    params = params,
                                    max_K = temp_max_K,
                                    sample_presence=sample_presence)
  
  # Format chains
  if (length(temp_samps_list) == 1) {
    samps_list <- list(formatChains(temp_samps_list))
    names(samps_list) <- "K1"
  } else {
    samps_list <- parallel::mclapply(temp_samps_list, formatChains,
                                     mc.cores = mc.cores)
  }
  
  # check whether: 
  # 1) number of mutations per cluster is at least min_mutation_per_cluster 
  # 2) difference between any two cluster less than cluster_diff_thresh 
  filtered_samps_list <- filterK(samps_list, min_mutation_per_cluster = min_mutation_per_cluster,
                                 cluster_diff_thresh = cluster_diff_thresh)
  
  # Calculate BIC and silhouette
  K_tested <- seq_len(length(filtered_samps_list))
  if (temp_max_K > 1) {
    box_indata <- getBoxInputData(temp_box, ploidy)
    
    bic_vec <- unname(unlist(parallel::mclapply(filtered_samps_list,
                                                  function(chains) calcChainBIC(chains=chains, input.data=box_indata, pattern=temp_box$pattern),
                                                  mc.cores = mc.cores)))
    bic_tb <- tibble(K_tested = K_tested,
                     BIC = bic_vec)
    BIC_best_chains <- samps_list[[which.min(bic_vec)]]
    sc_vec <- unname(unlist(parallel::mclapply(filtered_samps_list,
                                                function(chains) calcChainSilhouette(chains=chains, input.data=box_indata, pattern=temp_box$pattern),
                                                mc.cores = mc.cores)))
    
    sc_tb <- tibble(K_tested = K_tested,
                     silhouette = sc_vec)
    sc_best_chains <- samps_list[[which.max(sc_vec)]]
    res_list <- list(all_chains = samps_list,
                     silhouette = sc_tb,
                     BIC = bic_tb,
                     BIC_best_chains = BIC_best_chains,
                     sc_best_chains = sc_best_chains,
                     BIC_best_K = which.min(bic_vec),
                     silhouette_best_K = which.max(sc_vec))
  } else {
    # only 1 variant, so must be 1 cluster and don't need to check BIC
    res_list <- list(all_chains = filtered_samps_list,
                     silhouette = NA,
                     BIC = NA,
                     BIC_best_chains = filtered_samps_list[[1]],
                     sc_best_chains = filtered_samps_list[[1]],
                     BIC_best_K = 1,
                     silhouette_best_K = 1)
  }
  return(res_list)
}

runMCMCForABox <- function(box, 
                           ploidy=2,
                           n.iter = 10000, 
                           n.burn = 1000, 
                           thin = 10, 
                           mc.cores = 1,
                           inits = list(".RNG.name" = "base::Wichmann-Hill",
                                        ".RNG.seed" = 123),
                           params = c("z", "mcf", "icn", "m", "ystar"),
                           max_K = 5, 
                           sample_presence=FALSE) {

  # select columns if the presence pattern is 1
  box_input_data <- getBoxInputData(box, ploidy)
  
  extdir <- system.file("extdata", package="pictograph")
  
  # choose sample in which mutations are present
  sample_to_sort <- which(colSums(box_input_data$y) > 0)[1]
  
  jags.file <- file.path(extdir, "model.jags")
  jags.file.K1 <- file.path(extdir, "model_K1.jags")
  
  
  samps_K1 <- runMCMC(box_input_data, 
                      1, 
                      jags.file.K1,
                      inits, 
                      params, 
                      n.iter=n.iter, 
                      thin=thin, 
                      n.burn=n.burn)
  
  if(box_input_data$S == 1) {
    colnames(samps_K1[[1]])[which(colnames(samps_K1[[1]]) == "mcf")] <- "mcf[1,1]"
  }
  
  if (sample_presence) {
    samps_K1 <- reverseDrop(samps_K1, box$pattern, n.iter)
  }

  if (max_K > 1) {
    
    box_input_data$sample_to_sort <- sample_to_sort
    
    samps_2 <- parallel::mclapply(2:max_K,
                                  function(k) runMCMC(box_input_data, k,
                                                      jags.file, inits, params,
                                                      n.iter=n.iter, thin=thin,
                                                      n.burn=n.burn),
                                  mc.cores=mc.cores)
    
    if (sample_presence) {
      for (i in seq_len(length(samps_2))) {
        samps_2[[i]] <- reverseDrop(samps_2[[i]], box$pattern, n.iter)
      }
    }
    samps_list <- c(list(samps_K1), samps_2)
    names(samps_list) <- paste0("K", 1:max_K)
    return(samps_list)
    
  } else {
    names(samps_K1) <- "K1"
    return(samps_K1)
  }
  
}

getBoxInputData <- function(box, ploidy=2) {
  sample_list = vector()
  
  # include samples if the pattern is 1; i.e. presence of mutations in the sample
  for (j in 1:ncol(box$y)) {
    if (strsplit(box$pattern, "")[[1]][j] == "1") {
      sample_list <- append(sample_list, j)
    }
  }
  box_input_data <- list(I = nrow(box$y),
                         S = length(sample_list),
                         y = box$y[,sample_list,drop=FALSE],
                         n = box$n[,sample_list,drop=FALSE],
                         tcn = box$tcn[,sample_list,drop=FALSE],
                         is_cn = box$is_cn,
                         mtp = box$mtp,
                         cncf = box$cncf[,sample_list,drop=FALSE],
                         icn = box$icn,
                         purity=box$purity,
                         ploidy=ploidy)
  
  # set tcn to 2 if 0
  box_input_data$tcn[box_input_data$tcn==0] <- 2
  return(box_input_data)
}

runMCMC <- function(box_input_data, 
                    K, 
                    jags.file, 
                    inits, 
                    params,
                    n.iter=10000, 
                    thin=10, 
                    n.chains=1,
                    n.adapt=1000, 
                    n.burn=1000) {
  if (K > 1) box_input_data$K <- K
  jags.m <- jags.model(jags.file,
                       box_input_data,
                       n.chains = n.chains,
                       inits = inits,
                       n.adapt = n.adapt)
  if (n.burn > 0) update(jags.m, n.burn)
  samps <- coda.samples(jags.m, params, n.iter=n.iter, thin=thin)
  
  return(samps)
}

reverseDrop <- function(samps, pattern, n.iter) {
  total_sample = nchar(pattern)
  sample_list = vector()
  for (j in seq_len(nchar(pattern))) {
    if (strsplit(pattern, "")[[1]][j] == "1") {
      sample_list <- append(sample_list, j)
    }
  }
  k_list = vector()
  ystar_list = vector()
  # replace current sample id by true sample id from pattern
  for (i in seq_len(length(colnames(samps[[1]])))) {
    if (startsWith(colnames(samps[[1]])[i], "mcf")) {
      para <- str_extract_all(colnames(samps[[1]])[i], "[0-9]+")[[1]]
      colnames(samps[[1]])[i] <- paste("mcf[", para[1], ",", sample_list[strtoi(para[2])], "]", sep = "")
      k_list <- c(k_list, para[1])
    } else if (startsWith(colnames(samps[[1]])[i], "ystar")) {
      para <- str_extract_all(colnames(samps[[1]])[i], "[0-9]+")[[1]]
      colnames(samps[[1]])[i] <- paste("ystar[", para[1], ",", sample_list[strtoi(para[2])], "]", sep = "")
      ystar_list <- c(ystar_list, para[1])
    }
  }
  k_list <- unique(k_list)
  ystar_list <- unique(ystar_list)
  
  # add back dropped samples
  absent_sample <- vector()
  for (sample in seq_len(total_sample)) {
    if (! sample %in% sample_list) {
      absent_sample <- append(absent_sample, sample)
    }
  }
  for (k in seq_len(length(k_list))) {
    for (j in seq_len(length(absent_sample))) {
      col = paste("mcf[", k_list[k], ",", absent_sample[j], "]", sep = "")
      samps[[1]] <- cbind(samps[[1]], col=0)
      colnames(samps[[1]])[colnames(samps[[1]]) == 'col'] <- col
    }
  }
  for (ystar in seq_len(length(ystar_list))) {
    for (j in seq_len(length(absent_sample))) {
      col = paste("ystar[", ystar_list[ystar], ",", absent_sample[j], "]", sep = "")
      samps[[1]] <- cbind(samps[[1]], col=0)
      colnames(samps[[1]])[colnames(samps[[1]]) == 'col'] <- col
    }
  }
  
  samps[[1]] <- samps[[1]][,order(colnames(samps[[1]]))]
  
  return(samps)
}


filterK <- function(samps_list, min_mutation_per_cluster=1, cluster_diff_thresh=0.05) {
  filtered_samps_list <- list()
  toBreak = F
  for (i in seq_len(length(samps_list))) {
    k = as.numeric(gsub("\\D", "", names(samps_list)[i]))
    if (k > 1) {
      mcf_chain = samps_list[[i]]$mcf_chain
      z_chain = samps_list[[i]]$z_chain
      clusterTable = writeClusterAssignmentsTable(z_chain)
      # check whether all cluster contains at least one mutation
      if (length(unique(clusterTable$Cluster))==k) {
        # check whether all cluster contains at least min_mutation_per_cluster mutations
        if (any(table(clusterTable$Cluster) < min_mutation_per_cluster)) {
          break
        }
      } else {
        break
      }
      mcfTable = writeClusterMCFsTable(mcf_chain)
      # check whether mcf for any cluster is less than cluster_diff_thresh in all samples
      for (j1 in seq_len(k)) {
        if (all(mcfTable[j1,2:ncol(mcfTable)] < cluster_diff_thresh)) {
          toBreak = T
        }
      }
      # check whether mcf difference between any two clusters less than cluster_diff_thresh in all samples
      for (j1 in seq_len(k-1)) {
        for (j2 in seq(j1+1, k)) {
          diff = abs(mcfTable[j1,2:ncol(mcfTable)] - mcfTable[j2,2:ncol(mcfTable)])
          if (all(diff < cluster_diff_thresh)) {
            toBreak = T
          }
        }
      }
    }
    if (toBreak) { break }
    filtered_samps_list[[names(samps_list)[i]]] <- samps_list[[i]]
  }
  return(filtered_samps_list)
}

formatChains <- function(samps) {
  temp_z <- get.parameter.chain("z", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  if (nrow(temp_z) == 0) {
    temp_z <- get.parameter.chain("z", ggmcmc::ggs(samps) %>% mutate(Parameter = gsub("z","z[1]",Parameter))) %>%
      mutate(Parameter = as.character(Parameter))
  }
  temp_mcf <- get.parameter.chain("mcf", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  if (nrow(temp_mcf) == 0) {
    temp_mcf <- get.parameter.chain("mcf", ggmcmc::ggs(samps) %>% mutate(Parameter = gsub("mcf","mcf[1,1]",Parameter))) %>%
      mutate(Parameter = as.character(Parameter))
  }
  temp_ystar <- get.parameter.chain("ystar", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  temp_icn <- get.parameter.chain("icn", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  if (nrow(temp_icn) == 0) {
    temp_icn <- get.parameter.chain("icn", ggmcmc::ggs(samps) %>% mutate(Parameter = gsub("icn","icn[1]",Parameter))) %>%
      mutate(Parameter = as.character(Parameter))
  }
  temp_m <- get.parameter.chain("m", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  if (nrow(temp_m) == 0) {
    temp_m <- get.parameter.chain("m", ggmcmc::ggs(samps) %>% mutate(Parameter = gsub("\\bm\\b","m[1]",Parameter))) %>%
      mutate(Parameter = as.character(Parameter))
  }
  samps_list_formatted <- list(mcf_chain = temp_mcf,
                               z_chain = temp_z,
                               icn_chain = temp_icn,
                               m_chain = temp_m,
                               ystar_chain = temp_ystar)
  return(samps_list_formatted)
}

get.parameter.chain <- function(param, chains) {
  chains[grep(paste0(param, "\\["), chains$Parameter), ]
}
