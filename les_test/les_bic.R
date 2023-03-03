calcBIC_les <- function(samp_list, temp_box){
  
  BIC <- rep(0,length(samp_list))
  
  for (k in 1:length(samp_list)){
    chains <- samp_list[[k]]
    #get one chain from a box and its K
    # get w matrix
    options(dplyr.summarise.inform = FALSE)
    w <- chains$w_chain%>%
      mutate(value = round(value,5),
             Cluster = as.numeric(gsub("w\\[(.*),.*","\\1", Parameter)),
             Sample = as.numeric(gsub(".*,(.*)\\]","\\1", Parameter)))%>%
      group_by(Cluster,Sample)%>%
      summarize(w = mean(value))%>%
      ungroup()%>%
      spread(key = Sample, value = w)
    # 
    # samp_list[[3]]$w_chain%>%
    #   mutate(Cluster = as.numeric(gsub("w\\[(.*),.*","\\1", Parameter)),
    #          Sample = as.numeric(gsub(".*,(.*)\\]","\\1", Parameter)))%>%
    #   ggplot(aes(x = Iteration, y = value))+
    #   geom_line() +
    #   theme_light() +
    #   facet_wrap("Cluster"~"Sample") +
    #   ylab("Cancer Cell Fraction")
    
    # get mutation assignment
    ww <- writeClusterAssignmentsTable(chains$z_chain)%>%
      mutate(Mut_ID = as.numeric(gsub("Mut","",Mut_ID)))%>%
      arrange(Mut_ID)%>%
      left_join(w, by = "Cluster")%>%
      select(-c("Mut_ID","Cluster"))%>%
      as.matrix()
    
    y <- temp_box$y
    n <- temp_box$n
    m <- temp_box$m
    p <- temp_box$purity
    pp <- diag(p,length(p),length(p))
    c <- temp_box$tcn
    i <- temp_box$I
    s <- temp_box$S
    
    vaf <- m*ww%*%pp/sweep(c%*%pp, 2, -2*(1-p))
    
    # vaf is the probability
    # trails is n
    #assume predicted n is from the input
    
    ############################################################# do the log first then sum()
    lik <- sum(dbinom(y,n,vaf,log = T))
    BIC[k] <- log(i*s)*k-2*lik
  }
  return(BIC)
}


runMutSetMCMC <- function(temp_box, 
                          n.iter = 10000, n.burn = 1000, thin = 10, mc.cores = 1,
                          inits = list(".RNG.name" = "base::Wichmann-Hill",
                                       ".RNG.seed" = 123),
                          temp_max_K = 5,
                          model_type = "spike_and_slab",
                          params = c("z", "w", "ystar"),
                          beta.prior = FALSE) {
  
  # Run MCMC
  if (temp_max_K == 1) {
    # only 1 possible cluster
    temp_samps_list <- runMCMCForABox(temp_box,
                                      n.iter = n.iter, n.burn = n.burn, 
                                      thin = thin, mc.cores = mc.cores,
                                      inits = inits,
                                      params = params,
                                      max_K = temp_max_K)
  } else {
    # assess range of K: [1, temp_max_K]
    temp_samps_list <- runMCMCForABox(temp_box,
                                      n.iter = n.iter, n.burn = n.burn, 
                                      thin = thin, mc.cores = mc.cores,
                                      inits = inits,
                                      params = params,
                                      max_K = temp_max_K,
                                      model = model_type,
                                      beta.prior = beta.prior)
  }
  
  # Format chains
  samps_list <- parallel::mclapply(temp_samps_list, formatChains,
                                   mc.cores = mc.cores)
  
  # Calculate BIC
  K_tested <- seq_len(temp_max_K)
  if (temp_max_K > 1) {
    bic_vec <- calcBIC_les(samps_list,temp_box)
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