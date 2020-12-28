runTreeMH <- function(w_chain, 
                      max.num.root.children=1,
                      num_iter=10000, thin=10, burn=1000,
                      first_am=NULL, seed=123,
                      post=NULL,
                      mc.cores=1) {
  set.seed(seed)
  mcf_stats <- summarizeWChain(w_chain)
  mcf_matrix <- get.map.w(w_chain)
  cpov <- create.cpov(mcf_stats)
  
  # initialize chains ---------------------------------------------------------
  if (is.null(post)) {
    if (is.null(first_am)) {
      am_chain <- list(initializeGraph(mcf_matrix, max.num.root.children))
    } else {
      am_chain <- list(first_am)
    }
  } else {
    am_chain <- list(initializeGraphFromPost2(post, max.num.root.children))
  }
  
  if (length(filter(am_chain[[1]], possible_edge == TRUE)$child) == 
      length(unique(filter(am_chain[[1]], possible_edge == TRUE)$child))) {
    print("Possible edges:")
    print(am_chain[[1]] %>% filter(possible_edge == TRUE))
    print("Only one possible tree, so stopping MH and returning single tree")
    results <- list(am_chain = am_chain,
                    score_chain = calcTreeFitness(am_chain[[1]], cpov, mcf_matrix))
    return(results)
  }

  am_prev <- am_chain[[1]]
  score_chain <- calcTreeFitness(am_prev, cpov, mcf_matrix)
  fit_prev <- score_chain[1]
  
  # Metropolis Hastings -------------------------------------------------------
  num_accept <- 0
  
  u_vec <- log(runif(num_iter, min = 0, max = 1))
  for (i in seq_len(num_iter)) {
    # propse new am.long
    am_star <- sampleNewEdge(am_prev, 
                             max.num.root.children = max.num.root.children, 
                             mc.cores = mc.cores)
    if (all(am_star$connected == am_prev$connected)) {
      am_star <- initializeGraphFromPost2(post, max.num.root.children)
    }
    fit_star <- calcTreeFitness(am_star, cpov, mcf_matrix)
    
    # accept or reject proposal
    r <- log(fit_star) - log(fit_prev)
    #u <- log(runif(1, 0, 1))
    if(u_vec[i] <= r) {
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

runTreeMH2 <- function(filtered_trees, w_chain = NULL, mcf_matrix = NULL, 
                       num_iter=10000, thin=10, burn=1000,
                       first_am=NULL, seed=123,
                       mc.cores=1) {
  set.seed(seed)
  
  if (!is.null(w_chain)) {
    mcf_stats <- summarizeWChain(w_chain)
    mcf_matrix <- get.map.w(w_chain)
    cpov <- create.cpov(mcf_stats)
  } else if (!is.null(mcf_matrix)) {
    # make fake mcf_stats -- set sd to 0.05
    w_params <- paste0("w[", 
                       rep(1:nrow(mcf_matrix), each = ncol(mcf_matrix)), ",", 
                       rep(1:ncol(mcf_matrix), times = nrow(mcf_matrix)), "]")
    mcf_stats <- tibble(Parameter = factor(w_params, levels = w_params),
                        sd = 0.05,
                        mean = c(t(mcf_matrix)))
    cpov <- create.cpov(mcf_stats)
  } else {
    stop("Must supply either w_chain or mcf_matrix")
  }
  
  # don't run MH if only 1 possible tree --------------------------------------
  if (length(filtered_trees) == 1) {
    print("Only one possible tree, so stopping MH and returning single tree")
    results <- list(am_chain = filtered_trees,
                    score_chain = calcTreeFitness(filtered_trees[[1]], cpov, mcf_matrix,
                                                  am_format = "edges"))
    return(results)
  }
  # initalize chains -----------------------------------------------------------
  am_chain <- list()
  if (!is.null(first_am)) {
    am_chain[[1]] <- first_am
  } else{
    am_chain[[1]] <- edgesToAmLong(filtered_trees[[sample(seq_len(length(filtered_trees)), 1)]])
  }
  score_chain <- calcTreeFitness(am_chain[[1]], cpov, mcf_matrix)

  # if tree space is smaller than number of iterations, calculate score for each tree first
  if (length(filtered_trees) < num_iter) {
    all_scores <- parallel::mclapply(filtered_trees, 
                                     function(edges) calcTreeFitness(edges, cpov, mcf_matrix, am_format = "edges"),
                                     mc.cores = mc.cores)
    all_scores <- unlist(all_scores)
    
    results <- preCalcMH(filtered_trees, all_scores, num_iter, am_chain, score_chain,
                         thin, burn)
    
  } else {
    results <- calcEachIterMH(filtered_trees, cpov, mcf_matrix, num_iter, am_chain, score_chain,
                              thin, burn)
  }
  
  return(results)
}

preCalcMH <- function(filtered_trees, all_scores, num_iter, am_chain, score_chain,
                      thin, burn) {
  # MH for when tree scores are already calculated 
  num_accept <- 0
  am_prev <- am_chain[[1]]
  fit_prev <- score_chain[[1]]
  u_vec <- log(runif(num_iter, min = 0, max = 1))
  
  # sample tree indices
  tree_ind <- sample(seq_len(length(filtered_trees)), num_iter, replace = TRUE)
  
  for (i in seq_len(num_iter)) {
    # propse new am.long
    am_star <- filtered_trees[[tree_ind[i]]]
    fit_star <- all_scores[tree_ind[i]]
    
    # accept or reject proposal
    r <- log(fit_star) - log(fit_prev)
    if(u_vec[i] <= r) {
      am_prev <- edgesToAmLong(am_star)
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

calcEachIterMH <- function(filtered_trees, cpov, mcf_matrix, num_iter, am_chain, score_chain,
                           thin, burn) {

  # initial values 
  num_accept <- 0
  am_prev <- am_chain[[1]]
  fit_prev <- score_chain[[1]]
  
  u_vec <- log(runif(num_iter, min = 0, max = 1))
  # sample tree indices
  tree_ind <- sample(seq_len(length(filtered_trees)), num_iter, replace = TRUE)
  
  for (i in seq_len(num_iter)) {
    # propse new am.long
    am_star <- edgesToAmLong(filtered_trees[[tree_ind[1]]])
    fit_star <- calcTreeFitness(am_star, cpov, mcf_matrix)
    
    # accept or reject proposal
    r <- log(fit_star) - log(fit_prev)
    if(u_vec[i] <= r) {
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