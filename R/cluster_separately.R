#' Run MCMC to cluster mutations and estimate CCFs
#' 
#' @export
#' @importFrom ggmcmc ggs
#' @param input_data list of input data objects; 
#' @param n.iter number of iterations to run MCMC
#' @param n.burn number of iterations for burn in 
#' @param thin thinning parameter
#' @param mc.cores number of cores for parallelization
#' @param max_K maximum number of clusters to assess for each mutation set
#' @param model_type hierarchical model type for ("spike_and_slab" or "simple)
#' @param beta.prior option to run an initial MCMC chain and use results to specify beta priors for a second MCMC chain 
#' @param one_box option to run the MCMC chain without using sample presence
clusterSep <- function(input_data,
                       n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                       inits = list(".RNG.name" = "base::Wichmann-Hill",
                                    ".RNG.seed" = 123),
                       max_K = 5,
                       model_type = "spike_and_slab",
                       beta.prior = FALSE,
                       drop_zero = FALSE,
                       one_box = TRUE) {
  # 1. separate mutations by sample presence
  if (one_box) {
    input_data$mutation_indices <- seq_len(input_data$I)
    sep_list <- vector("list", 1)
    sep_list[[1]] <- input_data
    names(sep_list) <- "one_box"
  } else {
    sep_list <- separateMutationsBySamplePresence(input_data)
  }
  
  # 2a. For each presence set, run clustering MCMC, calc BIC and choose best K (min BIC)
  all_set_results <- vector("list", length(sep_list))
  names(all_set_results) <- names(sep_list)
  params = c("z", "w", "ystar")
  
  for (i in seq_len(length(sep_list))) {
    temp_box <- sep_list[[i]]
    # Max number of clusters cannot be more than number of mutations
    temp_max_K <- min(max_K, length(temp_box$mutation_indices))
    
    temp_samps_list <- runMutSetMCMC(temp_box, 
                                     n.iter = n.iter, n.burn = n.burn, thin = thin, 
                                     mc.cores = mc.cores,
                                     inits = inits,
                                     temp_max_K = temp_max_K,
                                     model_type = model_type,
                                     params = params,
                                     beta.prior = beta.prior,
                                     drop_zero = drop_zero)
    all_set_results[[i]] <- temp_samps_list
  }
  
  return(all_set_results)
}

separateMutationsBySamplePresence <- function(input_data) {
  # returns list of lists -- 
  # each item of list contains input data for a mutation sample presence set 
  # original mutation indices from input_data are recorded in $mutation_indices
  pres <- ifelse(input_data$y > 0, 1, 0)
  pat <- apply(pres, 1, function(x) paste0(x, collapse=""))
  types <- sort(names(table(pat)), decreasing=TRUE)
  if (length(types) == 1) {
    type_indices <- list()
    type_indices[[types]] <- seq_len(input_data$I)
  } else {
    type_indices <- lapply(types, function(x) which(pat == x))
    names(type_indices) <- types
  }
  
  sep_list <- list()
  for (t in seq_len(length(types))) {
    sep_list[[types[t]]] <- list(pattern = types[t],
                                 mutation_indices = type_indices[[types[t]]],
                                 purity = input_data$purity,
                                 I = length(type_indices[[types[t]]]),
                                 S = input_data$S,
                                 y = input_data$y[type_indices[[types[t]]], ,drop=FALSE],
                                 n = input_data$n[type_indices[[types[t]]], ,drop=FALSE],
                                 tcn = input_data$tcn[type_indices[[types[t]]], ,drop=FALSE],
                                 m = input_data$m[type_indices[[types[t]]], ,drop=FALSE])
    if (input_data$S == 1) {
      break
    }
  }
  return(sep_list)
}

getBoxInputData <- function(box) {
  list(purity = box$purity,
       I = box$I, 
       S = box$S,
       y = box$y,
       n = box$n,
       m = box$m,
       tcn = box$tcn)
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
    # print(colnames(samps[[1]])[i])
    if (startsWith(colnames(samps[[1]])[i], "w")) {
      para <- str_extract_all(colnames(samps[[1]])[i], "[0-9]+")[[1]]
      colnames(samps[[1]])[i] <- paste("w[", para[1], ",", sample_list[strtoi(para[2])], "]", sep = "")
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
      col = paste("w[", k_list[k], ",", absent_sample[j], "]", sep = "")
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

runMCMCForABox <- function(box, 
                           n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                           inits = list(".RNG.name" = "base::Wichmann-Hill",
                                        ".RNG.seed" = 123),
                           params = c("z", "w", "ystar"),
                           max_K = 5, model_type = "simple",
                           beta.prior = FALSE,
                           drop_zero = FALSE) {
  # returns samps_list 
  box_input_data <- getBoxInputData(box)
  
  # modify box_input_data so it only contain non-zero samples
  if (drop_zero) {
    sample_list = vector()
    for (j in 1:box$S) {
      if (strsplit(box$pattern, "")[[1]][j] == "1") {
        sample_list <- append(sample_list, j)
      }
    }
    # print(temp_box)
    box_input_data$purity <- box_input_data$purity[sample_list]
    box_input_data$y <- box_input_data$y[,sample_list,drop=FALSE]
    box_input_data$n <- box_input_data$n[,sample_list,drop=FALSE]
    box_input_data$tcn <- box_input_data$tcn[,sample_list,drop=FALSE]
    box_input_data$m <- box_input_data$m[,sample_list,drop=FALSE]
    box_input_data$S <- length(sample_list)
  }
  
  extdir <- system.file("extdata", package="pictograph")
  if (box$I == 1) {
    jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1_I1.jags")
    box_input_data$I <- NULL
  } else {
    jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1.jags")
  }
  
  if (model_type == "simple") {
    jags.file <- file.path(extdir, "model-test.jags") # fixes order of CCFs in one sample, not spike and slab 
  } else if (model_type == "spike_and_slab") {
    jags.file <- file.path(extdir, "spike_and_slab_purity_ident.jags") # fixing order of CCFs in one sample
  } else stop("provide model_type either 'spike_and_slab' or 'simple'")
  
  # choose sample in which mutations are present
  sample_to_sort <- which(colSums(box_input_data$y) > 0)[1]
  
  samps_K1 <- runMCMC(box_input_data, 1, jags.file.K1, 
                      inits, params, n.iter=n.iter, thin=thin, n.burn=n.burn)
  
  if(box_input_data$S == 1) {
    colnames(samps_K1[[1]])[which(colnames(samps_K1[[1]]) == "w")] <- "w[1,1]"
  }
  
  if (drop_zero) {
    samps_K1 <- reverseDrop(samps_K1, box$pattern, n.iter)
  }
  
  if(box$I == 1) {
    colnames(samps_K1[[1]])[which(colnames(samps_K1[[1]]) == "z")] <- "z[1]"
  }
  
  # Max number of clusters cannot be more than number of mutations
  max_K <- min(max_K, length(box$mutation_indices)) 
  if (max_K > 1) {
  
    box_input_data$sample_to_sort <- sample_to_sort
    
    samps_2 <- parallel::mclapply(2:max_K,
                                  function(k) runMCMC(box_input_data, k,
                                                      jags.file, inits, params,
                                                      n.iter=n.iter, thin=thin,
                                                      n.burn=n.burn,
                                                      beta.prior=beta.prior),
                                  mc.cores=mc.cores)
    
    if (drop_zero) {
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

formatChains <- function(samps) {
  temp_z <- get.parameter.chain("z", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  temp_w <- get.parameter.chain("w", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  if (nrow(temp_w) == 0) {
    temp_w <- get.parameter.chain("w", ggmcmc::ggs(samps) %>% mutate(Parameter = gsub("w","w[1,1]",Parameter))) %>%
      mutate(Parameter = as.character(Parameter))
  }
  temp_ystar <- get.parameter.chain("ystar", ggmcmc::ggs(samps)) %>%
    mutate(Parameter = as.character(Parameter))
  samps_list_formatted <- list(w_chain = temp_w,
                               z_chain = temp_z,
                               ystar_chain = temp_ystar)
  return(samps_list_formatted)
}

runMutSetMCMC <- function(temp_box, 
                          n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                          inits = list(".RNG.name" = "base::Wichmann-Hill",
                                       ".RNG.seed" = 123),
                          temp_max_K = 5,
                          model_type = "spike_and_slab",
                          params = c("z", "w", "ystar"),
                          beta.prior = FALSE,
                          drop_zero = FALSE) {
  
  # Run MCMC
  if (temp_max_K == 1) {
    # only 1 possible cluster
    temp_samps_list <- runMCMCForABox(temp_box,
                                      n.iter = n.iter, n.burn = n.burn, 
                                      thin = thin, mc.cores = mc.cores,
                                      inits = inits,
                                      params = params,
                                      max_K = temp_max_K,
                                      drop_zero = drop_zero)
  } else {
    # assess range of K: [1, temp_max_K]
    temp_samps_list <- runMCMCForABox(temp_box,
                                      n.iter = n.iter, n.burn = n.burn, 
                                      thin = thin, mc.cores = mc.cores,
                                      inits = inits,
                                      params = params,
                                      max_K = temp_max_K,
                                      model = model_type,
                                      beta.prior = beta.prior,
                                      drop_zero = drop_zero)
  }
  
  # Format chains
  if (drop_zero && temp_box$I == 1) {
    samps_list <- list(formatChains(temp_samps_list))
    names(samps_list) <- "K1"
  } else {
    samps_list <- parallel::mclapply(temp_samps_list, formatChains,
                                     mc.cores = mc.cores)
  }
  
  # Calculate BIC
  K_tested <- seq_len(temp_max_K)
  if (temp_max_K > 1) {
    box_indata <- getBoxInputData(temp_box)
    bic_vec <- unname(unlist(parallel::mclapply(samps_list, 
                                         function(chains) calcChainBIC(chains, box_indata),
                                         mc.cores = mc.cores)))
    bic_tb <- tibble(K_tested = K_tested,
                     BIC = bic_vec)
    best_chains <- samps_list[[which.min(bic_vec)]]
    res_list <- list(all_chains = samps_list,
                     BIC = bic_tb,
                     best_chains = best_chains,
                     best_K = which.min(bic_vec))
  } else {
    # only 1 variant, so must be 1 cluster and don't need to check BIC
    res_list <- list(all_chains = samps_list,
                     BIC = NA,
                     best_chains = samps_list[[1]],
                     best_K = 1)
  }
  
  return(res_list)
}






#' Collect chains for best K of each mutation set 
#' 
#' @export
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param chosen_K (Optional) Vector of K to choose for each mutation set, in the same order as all_set_results. If left blank, function will select best K automatically selected by \code{clusterSep}
collectBestKChains <- function(all_set_results, chosen_K = NULL) {
  if (is.null(chosen_K)) {
    best_set_chains <- lapply(all_set_results, function(x) x$best_chains)
  } else {
    best_set_chains <- mapply(function(set_res, choose_K) set_res$all_chains[[choose_K]],
                              set_res = all_set_results,
                              chosen_K,
                              SIMPLIFY = FALSE)
  }
  return(best_set_chains)
}

#' Relabel chains for all sets and merge 
#' 
#' @export
#' @import dplyr
#' @param best_set_chains List of lists of MCMC chains (w_chain, z_chain, ystar_chain) for each mutation set
#' @param indata List of input data objects (same as provided to clusterSep)
mergeSetChains <- function(best_set_chains, indata) {
  best_K_vals <- unname(sapply(best_set_chains, function(x) max(x$z_chain$value)))
  sep_list <- separateMutationsBySamplePresence(indata)
  
  # first set doesn't need to change cluster labels
  w_chain <- best_set_chains[[1]]$w_chain
  temp_z_chain <- best_set_chains[[1]]$z_chain
  temp_ystar_chain <- best_set_chains[[1]]$ystar_chain
  
  if (length(best_set_chains) > 1) {
    # still need to change mutation indices if more than 1 box
    z_chain <- relabel_z_chain_mut_only(temp_z_chain, sep_list[[1]]$mutation_indices)
    ystar_chain <- relabel_ystar_chain(temp_ystar_chain,
                                       sep_list[[1]]$mutation_indices)
    for (i in 2:length(best_set_chains)) {
      temp_w_chain <- best_set_chains[[i]]$w_chain
      temp_z_chain <- best_set_chains[[i]]$z_chain
      temp_ystar_chain <- best_set_chains[[i]]$ystar_chain
      new_cluster_labels <- seq_len(best_K_vals[i]) + sum(best_K_vals[1:(i-1)])
      
      temp_relabeled_w_chain <- relabel_w_chain(temp_w_chain, new_cluster_labels)
      temp_relabeled_z_chain <- relabel_z_chain(temp_z_chain, new_cluster_labels, 
                                                sep_list[[i]]$mutation_indices)
      temp_relabeled_ystar_chain <- relabel_ystar_chain(temp_ystar_chain,
                                                        sep_list[[i]]$mutation_indices)
      
      w_chain <- rbind(w_chain, temp_relabeled_w_chain)
      z_chain <- rbind(z_chain, temp_relabeled_z_chain)
      ystar_chain <- rbind(ystar_chain, temp_relabeled_ystar_chain)
    }
  } else {
    z_chain <- temp_z_chain
    ystar_chain <- temp_ystar_chain
  }
  
  # set levels for Parameter
  w_chain <- w_chain %>% 
    mutate(k = as.numeric(gsub("w\\[", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][2])))) %>%
    arrange(k, s) %>%
    mutate(Parameter = factor(Parameter, levels = unique(w_chain$Parameter))) %>%
    select(Iteration, Chain, Parameter, value)
  
  z_chain_param_order <- tibble(Parameter = unique(z_chain$Parameter)) %>%
    mutate(Variant = as.numeric(gsub("z\\[", "", 
                                     gsub("\\]", "", 
                                          unique(z_chain$Parameter))))) %>%
    arrange(Variant)
  z_chain <- z_chain %>%
    mutate(Parameter = factor(Parameter, levels = z_chain_param_order$Parameter))
  
  ystar_chain <- ystar_chain %>%
    mutate(Mutation_index = as.numeric(gsub("ystar\\[", "",
                                            sapply(ystar_chain$Parameter,
                                                   function(x) strsplit(as.character(x), ",")[[1]][1]))),
           s = as.numeric(gsub("\\]", "",
                               sapply(ystar_chain$Parameter,
                                      function(x) strsplit(as.character(x), ",")[[1]][2]))))
  ystar_chain <- ystar_chain %>%
    arrange(Mutation_index, s) %>%
    mutate(Parameter = factor(Parameter, levels = unique(Parameter)))
  
  chains <- list(w_chain = w_chain,
                 z_chain = z_chain,
                 ystar_chain = ystar_chain)
  return(chains)
}

relabel_w_chain <- function(w_chain, new_cluster_labels) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  new_w <- w_chain %>% 
    mutate(k = as.numeric(gsub("w\\[", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", 
                               sapply(w_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][2]))))
  if (length(new_cluster_labels) != length(unique(new_w$k))) {
    stop("number of supplied new cluster labels does not match the number of clusters in w_chain")
  }
  new_w <- new_w %>% 
    mutate(k_new = new_cluster_labels[new_w$k]) %>%
    mutate(Parameter = paste0("w[", k_new, ",", s, "]")) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_w)
}

relabel_z_chain <- function(z_chain, new_cluster_labels, mutation_indices) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  if (length(mutation_indices) != length(unique(z_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in z_chain")
  }
  ## would break when no mutation is assigned to a cluster
  ## poor choice of k, would prob lower the k
  if (length(new_cluster_labels) < length(unique(z_chain$value))) {
    stop("number of supplied new cluster labels does not match the number of clusters in z_chain")
  }
  new_z <- z_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("z\\[", "", 
                                    sapply(z_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_z <- new_z %>%
    mutate(new_i = mutation_indices[i],
           value = new_cluster_labels[new_z$value]) %>%
    mutate(Parameter = paste0("z[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_z)
}

relabel_ystar_chain <- function(ystar_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  i_s <- gsub("ystar\\[|]", "", ystar_chain$Parameter)
  i <- sapply(i_s, function(x) strsplit(x, ",")[[1]][1]) %>%
    as.numeric
  s <- sapply(i_s, function(x) strsplit(x, ",")[[1]][2]) %>%
    as.numeric
  new_ystar <- ystar_chain %>%
    mutate(i = i,
           s = s)
  new_ystar <- new_ystar %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("ystar[", new_i, ",", s, "]")) %>%
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_ystar)
}

relabel_z_chain_mut_only <- function(z_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  # cluster labels are left unchanged 
  if (length(mutation_indices) != length(unique(z_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in z_chain")
  }
  new_z <- z_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("z\\[", "", 
                                    sapply(z_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_z <- new_z %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("z[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_z)}