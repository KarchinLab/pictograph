#' Determine most probable integer assignments by taking those with highest posterior probability
#' 
#' @export
#' @param icn_chain MCMC chain of integer copy number assignment values
estimateICN <- function(icn_chain) {
  it <- max(icn_chain$Iteration)
  mcmc_icn <- icn_chain %>%
    group_by(Parameter, value) %>%
    reframe(n=n(),
            maxiter=it) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_icn <- mcmc_icn %>%
    group_by(Parameter) %>%
    reframe(value=value[probability==max(probability)]) %>%
    ungroup()
  
  # choose first cluster if equal probability
  map_icn_count <- map_icn %>% 
    group_by(Parameter) %>%
    reframe(map_count = n()) %>%
    ungroup()
  if (any(map_icn_count$map_count > 1)) {
    mut_ind <- which(map_icn_count$map_count > 1)
    for (i in mut_ind) {
      dup_var <- as.numeric(gsub("icn\\[|]", "", map_icn_count$Parameter[i]))
      map_icn_dups <- which(gsub("icn\\[|]", "", map_icn$Parameter) == dup_var)
      dup_ind <- map_icn_dups[-1]
      map_icn <- map_icn[-dup_ind, ]
    }
  }
  return(map_icn)
}
#' Determine most probable integer copy number by taking those with highest posterior probability. 
#' 
#' @export
#' @param icn_chain MCMC chain of integer copy number
#' @param Mut_ID Vector of mutation IDs, same order as provided as input data (e.g. indata$Mut_ID)
#' @return A tibble listing mutation IDs and their cluster assignments
writeIcnTable <- function(icn_chain, Mut_ID = NULL) {
  map_icn <- estimateICN(icn_chain) 
  if (is.null(Mut_ID)) {
    Mut_ID <- paste0("Mut", 1:nrow(map_icn))
  }
  map_icn <- map_icn %>%
    mutate(Parameter_n = as.numeric(gsub("icn\\[(\\d+)\\]","\\1",Parameter)))%>%
    arrange(Parameter_n)%>%
    mutate(Mut_ID = Mut_ID, icn = value)%>%
    select(Mut_ID, icn)
  
  # map_icn <- map_icn %>%
  #   arrange(Cluster)
  return(map_icn)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability
#' 
#' @export
#' @param m_chain MCMC chain of multiplicity assignment values
estimateMultiplicity <- function(m_chain) {
  it <- max(m_chain$Iteration)
  mcmc_m <- m_chain %>%
    group_by(Parameter, value) %>%
    reframe(n=n(),
            maxiter=it) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_m <- mcmc_m %>%
    group_by(Parameter) %>%
    reframe(value=value[probability==max(probability)]) %>%
    ungroup()
  
  # choose first cluster if equal probability
  map_m_count <- map_m %>% 
    group_by(Parameter) %>%
    reframe(map_count = n()) %>%
    ungroup()
  if (any(map_m_count$map_count > 1)) {
    mut_ind <- which(map_m_count$map_count > 1)
    for (i in mut_ind) {
      dup_var <- as.numeric(gsub("m\\[|]", "", map_m_count$Parameter[i]))
      map_m_dups <- which(gsub("m\\[|]", "", map_m$Parameter) == dup_var)
      dup_ind <- map_m_dups[-1]
      map_m <- map_m[-dup_ind, ]
    }
  }
  return(map_m)
}

#' Determine most probable multiplicity assignments by taking those with highest posterior probability. 
#' 
#' @export
#' @param m_chain MCMC chain of mutation cluster assignment values
#' @param Mut_ID Vector of mutation IDs, same order as provided as input data (e.g. indata$Mut_ID)
#' @return A tibble listing mutation IDs and their cluster assignments
writeMultiplicityTable <- function(m_chain, Mut_ID = NULL) {
  map_m <- estimateMultiplicity(m_chain) 
  if (is.null(Mut_ID)) {
    Mut_ID <- paste0("Mut", 1:nrow(map_m))
  }
  map_m <- map_m %>%
    mutate(Parameter_n = as.numeric(gsub("m\\[(\\d+)\\]","\\1",Parameter)))%>%
    arrange(Parameter_n)%>%
    mutate(Mut_ID = Mut_ID, Multiplicity = value)%>%
    select(Mut_ID, Multiplicity)
  
  return(map_m)
}

#' Determine the most probable cluster MCF values by taking the mean of the posterior distributions
#' 
#' @export
#' @param mcf_chain MCMC chain of mCF values
#' @return matrix of estimated cluster MCFs
estimateMCFs <- function(mcf_chain) {
  S <- numberSamples(mcf_chain)
  K <- numberClusters(mcf_chain)
  temp <- mcf_chain %>% 
    mutate(
      I = as.numeric(gsub("mcf\\[([0-9]+),[0-9]+\\]", "\\1", Parameter)),
      J = as.numeric(gsub("mcf\\[[0-9]+,([0-9]+)\\]", "\\1", Parameter))) %>% 
    group_by(I,J) %>% 
    summarise(mean_value = mean(value), .groups = 'drop')
  mcf.map.matrix <- matrix(NA, nrow = K, ncol = S, dimnames = list(1:K, 1:S))
  for(row in 1:nrow(temp)) {
    mcf.map.matrix[temp$I[row], temp$J[row]] <- round(temp$mean_value[row],3)
  }
  return(mcf.map.matrix)
}

#' @importFrom stringr str_replace
numberSamples <- function(mcf_stats){
  params <- as.character(mcf_stats$Parameter)    
  nSamples <- strsplit(params, ",") %>%
    sapply("[", 2) %>%
    stringr::str_replace("\\]", "") %>%
    as.numeric() %>%
    max()
  nSamples
}

#' @importFrom stringr str_replace
numberClusters <- function(mcf_stats){
  params <- as.character(mcf_stats$Parameter)
  K <- strsplit(params, ",") %>%
    sapply("[", 1) %>%
    str_replace("mcf\\[", "") %>%
    as.numeric() %>%
    max()
  K
}

#' Determine the most probable cluster MCF values by taking the mode of the posterior distributions
#' 
#' @export
#' @param mcf_chain MCMC chain of MCF values
#' @param Sample_ID Vector of sample IDs, same order as provided as input data (e.g. indata$Sample_ID)
#' @return A tibble of estimated cluster MCFs in each sample 
writeClusterMCFsTable <- function(mcf_chain, Sample_ID = NULL) {
  map_mcf <- as.data.frame(estimateMCFs(mcf_chain))
  
  if (is.null(Sample_ID)) {
    Sample_ID <- paste0("Sample ", 1:ncol(map_mcf))
  }
  colnames(map_mcf) <- Sample_ID
  map_mcf <- map_mcf %>%
    as_tibble() %>%
    bind_cols(tibble(Cluster = 1:nrow(map_mcf)), .)
  return(map_mcf)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability
#' 
#' @export
#' @param z_chain MCMC chain of mutation cluster assignment values
estimateClusterAssignments <- function(z_chain) {
  it <- max(z_chain$Iteration)
  mcmc_z <- z_chain %>%
    group_by(Parameter, value) %>%
    reframe(n=n(),
              maxiter=it) %>%
    mutate(probability=n/maxiter) %>%
    ungroup()
  map_z <- mcmc_z %>%
    group_by(Parameter) %>%
    reframe(value=value[probability==max(probability)]) %>%
    ungroup()
  
  # choose first cluster if equal probability
  map_z_count <- map_z %>% 
    group_by(Parameter) %>%
    reframe(map_count = n()) %>%
    ungroup()
  if (any(map_z_count$map_count > 1)) {
    mut_ind <- which(map_z_count$map_count > 1)
    for (i in mut_ind) {
      dup_var <- as.numeric(gsub("z\\[|]", "", map_z_count$Parameter[i]))
      map_z_dups <- which(gsub("z\\[|]", "", map_z$Parameter) == dup_var)
      dup_ind <- map_z_dups[-1]
      map_z <- map_z[-dup_ind, ]
    }
  }
  return(map_z)
}

#' Determine most probable mutation cluster assignments by taking those with highest posterior probability. 
#' 
#' @export
#' @param z_chain MCMC chain of mutation cluster assignment values
#' @param Mut_ID Vector of mutation IDs, same order as provided as input data (e.g. indata$Mut_ID)
#' @return A tibble listing mutation IDs and their cluster assignments
writeClusterAssignmentsTable <- function(z_chain, mcf_chain=NULL, cncf=NULL, Mut_ID = NULL) {
  map_z <- estimateClusterAssignments(z_chain) 
  if (is.null(Mut_ID)) {
    Mut_ID <- paste0("Mut", 1:nrow(map_z))
  }
  map_z <- map_z %>%
    mutate(Parameter_n = as.numeric(gsub("z\\[(\\d+)\\]","\\1",Parameter)))%>%
    arrange(Parameter_n)%>%
    mutate(Mut_ID = Mut_ID, Cluster = value)%>%
    select(Mut_ID, Cluster)

  if (!is.null(cncf)) {
    if (is.null(mcf_chain)) {
      warning("mcf_chain information is required to add CNA to cluster assignment table")
    } else {
      w_mat <- estimateMCFs(mcf_chain)
      for (i in seq_len(nrow(cncf))) {
        cls = which(apply(w_mat, 1, function(x) return(all(x == cncf_update[i,]))))
        map_z <- map_z %>% add_row(Mut_ID=rownames(cncf)[i], Cluster=cls)
      }
    }
  }

  return(map_z)
}

#' Collect chains for best K of each mutation set 
#' 
#' @export
#' @param all_set_results List of MCMC results for each mutation set
#' @param chosen_K (Optional) Vector of K to choose for each mutation set, in the same order as all_set_results. If left blank, function will select best K automatically selected by \code{clusterSep}
collectBestKChains <- function(all_set_results, chosen_K = NULL) {
  # best_set_chains <- lapply(all_set_results, function(x) x$all_chains[[length(x$all_chains)]])
  if (is.null(chosen_K)) {
    best_set_chains <- lapply(all_set_results, function(x) x$silhouette_best_chains)
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
#' @param best_set_chains List of lists of MCMC chains (mcf_chain, z_chain, ystar_chain) for each mutation set
#' @param indata List of input data objects (same as provided to clusterSep)
mergeSetChains <- function(best_set_chains, indata) {
  best_K_vals <- unname(sapply(best_set_chains, function(x) max(x$z_chain$value)))
  sep_list <- separateMutationsBySamplePresence(indata)
  
  # first set doesn't need to change cluster labels
  mcf_chain <- best_set_chains[[1]]$mcf_chain
  temp_m_chain <- best_set_chains[[1]]$m_chain
  temp_icn_chain <-  best_set_chains[[1]]$icn_chain
  temp_z_chain <- best_set_chains[[1]]$z_chain
  temp_ystar_chain <- best_set_chains[[1]]$ystar_chain
  
  if (length(best_set_chains) > 1) {
    # still need to change mutation indices if more than 1 box
    z_chain <- relabel_z_chain_mut_only(temp_z_chain, sep_list[[1]]$mutation_indices)
    m_chain <- relabel_m_chain_mut_only(temp_m_chain, sep_list[[1]]$mutation_indices)
    icn_chain <- relabel_icn_chain_mut_only(temp_icn_chain, sep_list[[1]]$mutation_indices)
    ystar_chain <- relabel_ystar_chain(temp_ystar_chain,
                                       sep_list[[1]]$mutation_indices)
    for (i in 2:length(best_set_chains)) {
      temp_mcf_chain <- best_set_chains[[i]]$mcf_chain
      temp_m_chain <- best_set_chains[[i]]$m_chain
      temp_icn_chain <-  best_set_chains[[i]]$icn_chain
      temp_z_chain <- best_set_chains[[i]]$z_chain
      temp_ystar_chain <- best_set_chains[[i]]$ystar_chain
      new_cluster_labels <- seq_len(best_K_vals[i]) + sum(best_K_vals[1:(i-1)])

      temp_relabeled_mcf_chain <- relabel_mcf_chain(temp_mcf_chain, new_cluster_labels)
      temp_relabeled_icn_chain <- relabel_icn_chain(temp_icn_chain, new_cluster_labels, 
                                                sep_list[[i]]$mutation_indices)
      temp_relabeled_m_chain <- relabel_m_chain(temp_m_chain, new_cluster_labels, 
                                                sep_list[[i]]$mutation_indices)
      temp_relabeled_z_chain <- relabel_z_chain(temp_z_chain, new_cluster_labels, 
                                                sep_list[[i]]$mutation_indices)
      temp_relabeled_ystar_chain <- relabel_ystar_chain(temp_ystar_chain,
                                                        sep_list[[i]]$mutation_indices)
      
      mcf_chain <- rbind(mcf_chain, temp_relabeled_mcf_chain)
      icn_chain <- rbind(icn_chain, temp_relabeled_icn_chain)
      m_chain <- rbind(m_chain, temp_relabeled_m_chain)
      z_chain <- rbind(z_chain, temp_relabeled_z_chain)
      ystar_chain <- rbind(ystar_chain, temp_relabeled_ystar_chain)
    }
  } else {
    z_chain <- temp_z_chain
    icn_chain <- temp_icn_chain
    m_chain <- temp_m_chain
    ystar_chain <- temp_ystar_chain
  }
  
  # set levels for Parameter
  mcf_chain <- mcf_chain %>% 
    mutate(k = as.numeric(gsub("mcf\\[", "", 
                               sapply(mcf_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", 
                               sapply(mcf_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][2])))) %>%
    arrange(k, s) %>%
    mutate(Parameter = factor(Parameter, levels = unique(mcf_chain$Parameter))) %>%
    select(Iteration, Chain, Parameter, value)
  
  z_chain_param_order <- tibble(Parameter = unique(z_chain$Parameter)) %>%
    mutate(Variant = as.numeric(gsub("z\\[", "", 
                                     gsub("\\]", "", 
                                          unique(z_chain$Parameter))))) %>%
    arrange(Variant)
  z_chain <- z_chain %>%
    mutate(Parameter = factor(Parameter, levels = z_chain_param_order$Parameter))
  
  icn_chain_param_order <- tibble(Parameter = unique(icn_chain$Parameter)) %>%
    mutate(Variant = as.numeric(gsub("icn\\[", "", 
                                     gsub("\\]", "", 
                                          unique(icn_chain$Parameter))))) %>%
    arrange(Variant)
  icn_chain <- icn_chain %>%
    mutate(Parameter = factor(Parameter, levels = icn_chain_param_order$Parameter))
  
  m_chain_param_order <- tibble(Parameter = unique(m_chain$Parameter)) %>%
    mutate(Variant = as.numeric(gsub("m\\[", "", 
                                     gsub("\\]", "", 
                                          unique(m_chain$Parameter))))) %>%
    arrange(Variant)
  m_chain <- m_chain %>%
    mutate(Parameter = factor(Parameter, levels = m_chain_param_order$Parameter))
  
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
  
  chains <- list(mcf_chain = mcf_chain,
                 z_chain = z_chain,
                 icn_chain = icn_chain,
                 m_chain = m_chain,
                 ystar_chain = ystar_chain)
  return(chains)
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

relabel_m_chain <- function(m_chain, new_cluster_labels, mutation_indices) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  if (length(mutation_indices) != length(unique(m_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in m_chain")
  }
  new_m <- m_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("m\\[", "", 
                                    sapply(m_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_m <- new_m %>%
    mutate(new_i = mutation_indices[i],
           value = new_m$value) %>%
    mutate(Parameter = paste0("m[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_m)
}

relabel_icn_chain <- function(icn_chain, new_cluster_labels, mutation_indices) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  if (length(mutation_indices) != length(unique(icn_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in icn_chain")
  }
  new_icn <- icn_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("icn\\[", "", 
                                    sapply(icn_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_icn <- new_icn %>%
    mutate(new_i = mutation_indices[i],
           value = new_icn$value) %>%
    mutate(Parameter = paste0("icn[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_icn)
}

relabel_m_chain_mut_only <- function(m_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  # cluster labels are left unchanged 
  if (length(mutation_indices) != length(unique(m_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in m_chain")
  }
  new_m <- m_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("m\\[", "", 
                                    sapply(m_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_m <- new_m %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("m[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_m)
}

relabel_icn_chain_mut_only <- function(icn_chain, mutation_indices) {
  # mutation_indices = numeric vector of original mutation indices prior to separating by sample presence
  # cluster labels are left unchanged 
  if (length(mutation_indices) != length(unique(icn_chain$Parameter))) {
    stop("number of supplied mutation indices does not match the number of mutations in icn_chain")
  }
  new_icn <- icn_chain %>%
    mutate(i = as.numeric(gsub("\\]", "", 
                               gsub("icn\\[", "", 
                                    sapply(icn_chain$Parameter, 
                                           function(x) strsplit(as.character(x), ",")[[1]][1])))))
  new_icn <- new_icn %>%
    mutate(new_i = mutation_indices[i]) %>%
    mutate(Parameter = paste0("icn[", new_i, "]")) %>% 
    arrange(new_i) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_icn)
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

relabel_mcf_chain <- function(mcf_chain, new_cluster_labels) {
  # new_cluster_labels = numeric vector of labels that map to 1:length(new_cluster_labels)
  new_mcf <- mcf_chain %>% 
    mutate(k = as.numeric(gsub("mcf\\[", "", 
                               sapply(mcf_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(s = as.numeric(gsub("\\]", "", 
                               sapply(mcf_chain$Parameter, 
                                      function(x) strsplit(as.character(x), ",")[[1]][2]))))
  if (length(new_cluster_labels) != length(unique(new_mcf$k))) {
    stop("number of supplied new cluster labels does not match the number of clusters in mcf_chain")
  }
  new_mcf <- new_mcf %>% 
    mutate(k_new = new_cluster_labels[new_mcf$k]) %>%
    mutate(Parameter = paste0("mcf[", k_new, ",", s, "]")) %>%
    select(Iteration, Chain, Parameter, value)
  return(new_mcf)
}
