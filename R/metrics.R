# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Tree metric: proportion of correct ancestral relationships 

getParent <- function(node, am.long) {
  am.long <- am.long %>%
    mutate(parent = as.character(parent))
  e <- am.long %>% 
    filter(child == node, connected == 1)
  if (nrow(e) > 1) warning(paste0("node ", node, " has ", nrow(e), " parents"))
  return(e$parent[1])
}

getPathFromRoot <- function(node, am.long) {
  path <- node
  curr_parent <- getParent(node, am.long)
  while(curr_parent != "root") {
    path <- c(path, curr_parent)
    curr_parent <- getParent(curr_parent, am.long)
  }
  as.numeric(path)
}

getRelationshipTypeOfClusters <- function(cluster1, cluster2, am.long) {
  # input: cluster1 and cluster2 can be either characters or numeric 
  # returns type of ancestral relationship of cluster1 and cluster2
  #     1 = cluster1 and cluster2 are the same
  #     2 = cluster1 is ancestral to cluster2
  #     3 = cluster2 is ancestral to cluster1
  #     4 = cluster1 and cluster2 are on distinct branches
  
  if (is.numeric(cluster1)) cluster1 <- as.character(cluster1)
  if (is.numeric(cluster2)) cluster2 <- as.character(cluster2)
  
  # clusters are the same
  if (cluster1 == cluster2) return(1)
  
  # cluster1 is ancestral to cluster2
  temp_path <- getPathFromRoot(cluster2, am.long)
  if (cluster1 %in% temp_path) return(2)
  
  # cluster2 is ancestral to cluster1
  temp_path <- getPathFromRoot(cluster1, am.long)
  if (cluster2 %in% temp_path) return(3)
  
  # cluster1 and cluster2 are on distinct branches
  return(4)
}

getMutRelTb <- function(z, mut_ind, am) {
  # input: 
  # z = vector of cluster assignments for each mutation
  # mut_ind = vector of mutation indices e.g. 1:I
  # am = am.long format of graph 
  # output: tibble with columns mut1, mut2, cluster1, cluster2, relationship_type
  # relationship_type is 1:4 -- see getRelationshipTypeOfCluster()
  
  # mutation pairs are all pairwise combinations, not including self
  all_mut_pairs <- t(combn(mut_ind, 2)) 
  
  # cluster pairs are all permutations of 2, including self
  all_cluster_pairs <- rbind(t(combn(unique(z), 2)), t(combn(rev(unique(z)), 2)),
                             matrix(rep(unique(z), 2), nrow=length(unique(z)), byrow = F))
  all_cluster_pairs_tb <- tibble(cluster1 = all_cluster_pairs[, 1],
                                 cluster2 = all_cluster_pairs[, 2])
  
  true_rel_types <- mapply(function(c1, c2) getRelationshipTypeOfClusters(c1, c2, am), 
                           c1=all_cluster_pairs[, 1], c2=all_cluster_pairs[, 2])
  
  rel_tb <- tibble(cluster1 = all_cluster_pairs[, 1],
                   cluster2 = all_cluster_pairs[, 2],
                   relationship_type = true_rel_types)
  
  true_mut_rel <- tibble(mut1 = all_mut_pairs[, 1],
                         mut2 = all_mut_pairs[, 2],
                         cluster1 = z[all_mut_pairs[, 1]],
                         cluster2 = z[all_mut_pairs[, 2]])
  
  left_join(true_mut_rel, rel_tb, by = c("cluster1", "cluster2"))
}

calcPropRelationshipsCorrect <- function(test_mut_rel, true_mut_rel) {
  # input: tibbles from getMutRelTb()
  # returns proportion of relationships in test_mut_rel that match true_mut_rel
  test_comp <- left_join(true_mut_rel, test_mut_rel, by = c("mut1", "mut2"))
  rel_same <- test_comp$relationship_type.x == test_comp$relationship_type.y
  return(sum(rel_same) / length(rel_same))
}

calcTreeMetricSingleIter <- function(z, am, sim_data) {
  true_mut_rel <- getMutRelTb(sim_data$z, 1:sim_data$I, sim_data$am.long)
  if (length(unique(z$value)) == 1) {
    this_rel <- tibble(mut1 = true_mut_rel$mut1,
                       mut2 = true_mut_rel$mut2,
                       cluster1 = 1,
                       cluster2 = 1,
                       relationship_type = 1)
  } else {
    this_rel <- getMutRelTb(z$value, 1:sim_data$I, am)
  }
  
  prop_true <- calcPropRelationshipsCorrect(this_rel, true_mut_rel)
  return(prop_true)
}

calcTreeMetricChain <- function(z_chain_list, am_chain, sim_data,
                                mc.cores=1) {

  tree_metric <- parallel::mcmapply(function(z, am) calcTreeMetricSingleIter(z, am, sim_data), 
                        z_chain_list, am_chain,
                        mc.cores = mc.cores)
  return(tree_metric)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# CCF metrics (as described in PASTRI paper)


calcMetric1 <- function(true_w, w_star) {
  # Metric 1 -- A measure of sensitivity, matches the true clones to the 
  #   nearest reported clusters
  # input: matrices with ncol = number of samples, nrow = number of clusters
  # returns: numeric value 
  
  # number of clusters can be different, but number of samples must be equal
  if (ncol(true_w) != ncol(w_star)) {
    stop("true_w and w_star have different number of samples")
  }
  
  min_diff <- rep(Inf, nrow(true_w))
  
  for (i in 1:nrow(true_w)) {
    for (j in 1:nrow(w_star)) {
      temp_diff <- mean(abs(true_w[i, ] - w_star[j, ]))
      if (temp_diff < min_diff[i]) min_diff[i] <- temp_diff
    }
  }
  
  return(sum(min_diff))
}

calcMetric2 <- function(true_w, w_star) {
  # Metric 2 -- A measure of specificity, matches the reported clusters to
  #   the nearest true clones 
  # input: matrices with ncol = number of samples, nrow = number of clusters
  # returns: numeric value 
  
  # number of clusters can be different, but number of samples must be equal
  if (ncol(true_w) != ncol(w_star)) {
    stop("true_w and w_star have different number of samples")
  }
  
  min_diff <- rep(Inf, nrow(w_star))
  
  for (i in 1:nrow(w_star)) {
    for (j in 1:nrow(true_w)) {
      temp_diff <- mean(abs(true_w[i, ] - w_star[j, ]))
      if (temp_diff < min_diff[i]) min_diff[i] <- temp_diff
    }
  }
  
  return(sum(min_diff))
}

wTibbleToMatrix <- function(w_tb) {
  # input: w_chain tibble with columns Iteration, Parameter, value for a single Iteration
  # output: w_matrix with rows = clusters, columns = samples
  S <- numberSamples(w_tb)
  K <- numberClusters(w_tb)
  w_mat <- matrix(w_tb$value, nrow = K, ncol = S, byrow = T)
  return(w_mat)
}

calcMetric1Chain <- function(true_w, w_chain, mc.cores) {
  w_chain_list <- w_chain %>%
    group_by(Iteration) %>%
    group_split()
  w_chain_list_matrix <- parallel::mclapply(w_chain_list, 
                                  wTibbleToMatrix, 
                                  mc.cores = mc.cores)
  m1_chain <- parallel::mclapply(w_chain_list_matrix, 
                                 function(w_star) calcMetric1(true_w, w_star),
                                 mc.cores = mc.cores)
  return(unlist(m1_chain))
}

calcMetric2Chain <- function(true_w, w_chain, mc.cores) {
  w_chain_list <- w_chain %>%
    group_by(Iteration) %>%
    group_split()
  w_chain_list_matrix <- parallel::mclapply(w_chain_list, 
                                  wTibbleToMatrix, 
                                  mc.cores = mc.cores)
  m2_chain <- parallel::mclapply(w_chain_list_matrix, 
                       function(w_star) calcMetric2(true_w, w_star),
                       mc.cores = mc.cores)
  return(unlist(m2_chain))
}