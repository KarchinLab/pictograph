runTreeMH <- function(w_chain, 
                      max.num.root.children=1,
                      num_iter=10000, thin=10, burn=1000,
                      first_am=NULL, seed=123,
                      post=NULL,
                      mc.cores=1) {
  set.seed(seed)
  mcf_stats <- summarizeWChain(w_chain)
  mcf_matrix <- get.map.w(w_chain)
  
  # initialize chains ---------------------------------------------------------
  #if (is.null(post)) {
    if (is.null(first_am)) {
      am_chain <- list(initializeGraph(mcf_matrix, max.num.root.children))
    } else {
      am_chain <- list(first_am)
    }
  #} else {
    # @TODO update restrictions using posterior admat
  #}
  
  
  cpov <- create.cpov(mcf_stats)
  am_prev <- am_chain[[1]]
  score_chain <- calcTreeFitness(am_prev, cpov, mcf_matrix)
  fit_prev <- score_chain[1]
  
  # Metropolis Hastings -------------------------------------------------------
  num_accept <- 0
  
  for (i in seq_len(num_iter)) {
    # propse new am.long
    am_star <- sampleNewEdge(am_prev, 
                             max.num.root.children = max.num.root.children, 
                             mc.cores = mc.cores)
    fit_star <- calcTreeFitness(am_star, cpov, mcf_matrix)
    
    # accept or reject proposal
    r <- log(fit_star) - log(fit_prev)
    u <- log(runif(1, 0, 1))
    if(u <= r) {
      am_prev <- am_star
      num_accept <- num_accept + 1
      fit_prev <- fit_star
    }
    
    # keep thinned obs (and don't keep burn-in)
    if (i%%thin == 0 & i > burn) {
      am_chain <- append(am_chain, list(am_prev))
      score_chain <- c(score_chain, fit_prev)
    }
  }
  results <- list(am_chain = am_chain,
                  score_chain = score_chain,
                  accept_rate = num_accept/num_iter)
  return(results)
}

