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
                                 y = input_data$y[type_indices[[types[t]]], ],
                                 n = input_data$n[type_indices[[types[t]]], ],
                                 tcn = input_data$tcn[type_indices[[types[t]]], ],
                                 m = input_data$m[type_indices[[types[t]]], ])
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

runMCMCForABox <- function(box, 
                           n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                           inits = list(".RNG.name" = "base::Wichmann-Hill",
                                        ".RNG.seed" = 123),
                           params = c("z", "w"),
                           max_K = 5) {
  # returns samps_list 
  box_input_data <- getBoxInputData(box)
  
  extdir <- system.file("extdata", package="clone.tools")
  if (box$I == 1) {
    jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1_I1.jags")
    box_input_data$I <- NULL
  } else {
    jags.file.K1 <- file.path(extdir, "spike_and_slab_purity_2_K1.jags")
  }
  
  jags.file <- file.path(extdir, "spike_and_slab_purity_2.jags")
  
  
  
  samps_K1 <- runMCMC(box_input_data, 1, jags.file.K1, inits, params, n.iter=n.iter, thin=thin, n.burn=n.burn)
  if(box$I == 1) {
    colnames(samps_K1[[1]])[which(colnames(samps_K1[[1]]) == "z")] <- "z[1]"
  }
  
  # Max number of clusters cannot be more than number of mutations
  max_K <- min(max_K, length(box$mutation_indices)) 
  if (max_K > 1) {
    samps_2 <- parallel::mclapply(2:max_K,
                                  function(k) runMCMC(box_input_data, k, jags.file, inits, params,
                                                      n.iter=n.iter, thin=thin,
                                                      n.burn=n.burn),
                                  mc.cores=mc.cores)
    return(c(list(samps_K1), samps_2))
  } else {
    return(samps_K1)
  }
  
}

#' @importFrom ggmcmc ggs
clusterSep <- function(input_data,
                       n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                       inits = list(".RNG.name" = "base::Wichmann-Hill",
                                    ".RNG.seed" = 123),
                       params = c("z", "w"),
                       max_K = 5) {
  # 1. separate mutations by sample presence
  sep_list <- separateMutationsBySamplePresence(input_data)
  
  # 2. For each presence set, run clustering MCMC and choose K by minimum BIC 
  best_K_vals <- c()
  sep_samps_list <- list()
  for (i in seq_len(length(sep_list))) {
    temp_box <- sep_list[[i]]
    # Max number of clusters cannot be more than number of mutations
    temp_max_K <- min(max_K, length(temp_box$mutation_indices))
    if (temp_max_K == 1) {
      temp_samps_list <- runMCMCForABox(temp_box,
                                        n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores,
                                        inits = inits,
                                        params = params,
                                        max_K = temp_max_K)
      sep_samps_list[[i]] <- temp_samps_list
      best_K_vals[i] <- 1
    } else {
      temp_samps_list <- runMCMCForABox(temp_box,
                                        n.iter = n.iter, n.burn = n.burn, thin = thin, mc.cores = mc.cores,
                                        inits = inits,
                                        params = params,
                                        max_K = temp_max_K)
      # calc BIC and choose K by minimum BIC 
      temp_BIC <- mapply(function(samps, k)
        calcBIC(temp_box$I*temp_box$S, k, calcChainLogLik(samps, getBoxInputData(temp_box), k)),
        samps = temp_samps_list, k = 1:temp_max_K)
      temp_min_BIC_K <- which.min(temp_BIC)
      best_K_vals[i] <- temp_min_BIC_K
      sep_samps_list[[i]] <- temp_samps_list[[temp_min_BIC_K]]
    }
  }
  
  # 3. Relabel z_chain and w_chain for alls and merge 
  # first set doens't need to change cluster labels
  w_chain <- get.parameter.chain("w", ggs(sep_samps_list[[1]])) %>%
    mutate(Parameter = as.character(Parameter))
  temp_z_chain <- get.parameter.chain("z", ggs(sep_samps_list[[1]])) %>%
    mutate(Parameter = as.character(Parameter))
  # still need to change mutation indices
  z_chain <- relabel_z_chain_mut_only(temp_z_chain, sep_list[[1]]$mutation_indices)
  
  if (length(sep_list) > 1) {
    for (i in 2:length(sep_samps_list)) {
      temp_w_chain <- get.parameter.chain("w", ggs(sep_samps_list[[i]]))
      temp_z_chain <- get.parameter.chain("z", ggs(sep_samps_list[[i]]))
      new_cluster_labels <- seq_len(best_K_vals[i]) + sum(best_K_vals[1:(i-1)])
      
      temp_relabeled_w_chain <- relabel_w_chain(temp_w_chain, new_cluster_labels)
      temp_relabeled_z_chain <- relabel_z_chain(temp_z_chain, new_cluster_labels, sep_list[[i]]$mutation_indices)
      
      w_chain <- rbind(w_chain, temp_relabeled_w_chain)
      z_chain <- rbind(z_chain, temp_relabeled_z_chain)
    }
  }
  
  # set levels for Parameter
  w_chain <- w_chain %>%
    mutate(Parameter = factor(Parameter, levels = unique(w_chain$Parameter)))
  z_chain_param_order <- tibble(Parameter = unique(z_chain$Parameter)) %>%
    mutate(Variant = as.numeric(gsub("z\\[", "", gsub("\\]", "", unique(z_chain$Parameter))))) %>%
    arrange(Variant)
  z_chain <- z_chain %>%
    mutate(Parameter = factor(Parameter, levels = z_chain_param_order$Parameter))
  
  return(list(w_chain, z_chain))
}

relabel_w_chain <- function(w_chain, new_cluster_labels) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  new_w <- w_chain %>% 
    mutate(k = as.numeric(gsub("w\\[", "", sapply(w_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", sapply(w_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][2]))))
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
  if (length(new_cluster_labels) != length(unique(z_chain$value))) {
    stop("number of supplied new cluster labels does not match the number of clusters in z_chain")
  }
  new_z <- z_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", gsub("z\\[", "", sapply(z_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_z <- new_z %>%
    mutate(new_i = mutation_indices[i],
           value = new_cluster_labels[new_z$value]) %>%
    mutate(Parameter = paste0("z[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_z)
}

relabel_z_chain_mut_only <- function(z_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  # cluster labels are left unchanged 
  if (length(mutation_indices) != length(unique(z_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in z_chain")
  }
  new_z <- z_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", gsub("z\\[", "", sapply(z_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_z <- new_z %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("z[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_z)}