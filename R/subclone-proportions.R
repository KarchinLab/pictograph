#' Calculate proportions of subclones in each sample
#'
#' @param w_mat Matrix of CCF estimates (from \code{estimateCCFs})
#' @param tree_edges Tibble of tree edges with columns edge, parent, and child 
#' @export
#' @import dplyr
calcSubcloneProportions <- function(w_mat, tree_edges) {
  K <- nrow(w_mat)
  S <- ncol(w_mat)
  subclone_props <- matrix(NA, nrow = K, ncol = S)
  
  leaves <- seq_len(K)[sapply(as.character(seq_len(K)), function(x) isLeaf(x, tree_edges))]
  
  # subclone proportions of leaves is just CCFs
  for (leaf in leaves) {
    subclone_props[leaf, ] <- w_mat[leaf, ]
  }
  
  nodes_left <- as.character(seq_len(K))[sapply(as.character(seq_len(K)), function(x) !isLeaf(x, tree_edges))]
  new_children <- leaves
  
  # calculate proportions of other subclones (bottom up)
  while(length(nodes_left) > 0) {
    temp_nodes <- tree_edges %>% 
      filter(child %in% new_children, parent != "root") %>% 
      pull(parent) %>%
      unique()
    
    # check if nodes are branch points 
    for (node in temp_nodes) {
      if (isBranchPoint(node, tree_edges)) {
        
        branch_children <- tree_edges %>%
          filter(parent == node) %>%
          pull(child)
        
        # only add ccf if all children have been added
        if (any(branch_children %in% as.character(nodes_left))) {
          temp_nodes <- temp_nodes[-which(temp_nodes == node)]
          next
        }
        
        node_ccf <- w_mat[as.numeric(node), ]
        ######sum_child_ccfs <- colSums(w_mat[as.numeric(branch_children), ])
        # sum up proportions for all children (all generations)
        temp_child_prop <- lapply(branch_children, 
                                  function(x) calcDescendantProps(x, tree_edges, subclone_props) +
                                    subclone_props[as.numeric(x), ]) %>%
          unlist() %>% 
          matrix(., nrow = length(branch_children), byrow = T) %>%
          colSums()
        
        subclone_props[as.numeric(node), ] <- sapply(node_ccf - temp_child_prop, function(x) max(x, 0))
        nodes_left <- nodes_left[-which(nodes_left == node)]
        
      } else {
        # add ccf like normal 
        
        # sum up proportions for all children (all generations)
        temp_child_prop <- calcDescendantProps(node, tree_edges, subclone_props)
        
        # proportion is node ccf - sum(all children proportions)
        subclone_props[as.numeric(node), ] <- sapply(w_mat[as.numeric(node), ] - temp_child_prop, function(x) max(x, 0))
        nodes_left <- nodes_left[-which(nodes_left == node)]
      }
    }
    new_children <- temp_nodes
  }
  
 # @TODO normalize proportions -- columns should be no more than 1  
  subclone_props <- normalizeProps(subclone_props)
  
  return(subclone_props)
}

normalizeProps <- function(subclone_props) {
  props_sums <- colSums(subclone_props)
  if (all(props_sums == 1)) return(subclone_props)
  
  for (s in seq_len(ncol(subclone_props))) {
    if (props_sums[s] != 1) {
      new_col <- round(subclone_props[, s] / props_sums[s], 3)
      subclone_props[, s] <- new_col
    }
  }
  return(subclone_props)
}

calcDescendantProps <- function(node, tree_edges, subclone_props) {
  temp_child <- tree_edges %>%
    filter(parent == node) %>%
    pull(child)
  temp_child_prop <- subclone_props[as.numeric(temp_child), ]
  while(!isLeaf(temp_child, tree_edges)) {
    temp_child <- tree_edges %>%
      filter(parent == temp_child) %>%
      pull(child)
    temp_child_prop <- temp_child_prop + subclone_props[as.numeric(temp_child), ]
  }
  return(temp_child_prop)
}

isBranchPoint <- function(node, tree_edges) {
  node <- as.character(node) # make sure node is character
  children <- tree_edges %>%
    filter(parent == node) %>%
    pull(child)
  return(length(children) > 1)
}

isLeaf <- function(node, tree_edges) {
  node <- as.character(node) # make sure node is character
  children <- tree_edges %>%
    filter(parent == node) %>%
    pull(child)
  return(length(children) == 0)
}

#' Plot pie charts for subclone proportions in each sample
#'
#' @param subclone_props matrix of subclone proportions (returned from \code{calcSubcloneProportions})
#' @param sample_names (Optional) Vector of sample names. Should be in the order of columns of subclone_props
#' @export
#' @import ggplot2
#' @importFrom magrittr set_colnames
#' @importFrom viridis viridis
plotSubclonePie <- function(subclone_props, sample_names = NULL) {
  if (is.null(sample_names)) sample_names <- paste0("Sample ", 1:ncol(subclone_props))
  props_tb <- subclone_props %>%
    magrittr::set_colnames(sample_names) %>%
    as_tibble() %>%
    mutate(Subclone = factor(paste0("Clone ", 1:nrow(.)),
                             levels = paste0("Clone ", 1:nrow(.)))) %>%
    pivot_longer(cols = sample_names,
                 names_to = "Sample",
                 values_to = "Proportion")
  
  clone_colors <- viridis::viridis(nrow(subclone_props))
  ggplot(props_tb, aes(x="", y=Proportion, fill = Subclone)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = clone_colors, drop = F) +
    theme_void() +
    facet_wrap(~Sample)
  
}

#' Force CCFs to comply with lineage precedence and sum condition
#'
#' @param w_mat Matrix of CCF estimates (from \code{estimateCCFs})
#' @param tree_edges Tibble of tree edges with columns edge, parent, and child 
#' @export
#' @import dplyr
forceCCFs <- function(w_mat, tree_edges) {
  
  K <- nrow(w_mat)
  S <- ncol(w_mat)
  fixed_w_mat <- matrix(NA, nrow = K, ncol = S)
  
  curr_edges <- tree_edges %>%
    filter(parent == "root")
  
  while (nrow(curr_edges) > 0) {
    next_edges <- tibble()
    for (p in unique(curr_edges$parent)) {
      # grab parent node CCF
      if (p == "root") {
        parent_ccfs <- rep(1, S)
      } else {
        parent_ccfs <- fixed_w_mat[as.numeric(p), ]
      }
      
      temp_curr_edges <- curr_edges %>%
        filter(parent == p)
      
      if (isBranchPoint(p, tree_edges)) {
        child_ccfs <- w_mat[as.numeric(temp_curr_edges$child), ]
        sample_sums <- colSums(child_ccfs)
        for (s in 1:S) {
          if (sample_sums[s] > parent_ccfs[s]) {
            for (i in seq_len(length(temp_curr_edges$child))) {
              child_num <- as.numeric(temp_curr_edges$child)[i]
              fixed_w_mat[i, s] <- max(round(child_ccfs[i, s] / sample_sums[s] * parent_ccfs[s], 2),
                                      0.01)
            }
          } else {
            for (child in as.numeric(temp_curr_edges$child)) {
              child_ccfs <- w_mat[child, ]
              fixed_ccfs <- mapply(function(child_ccf, parent_ccf) ifelse(child_ccf > parent_ccf, parent_ccf, child_ccf),
                                   child_ccf = child_ccfs,
                                   parent_ccf = parent_ccfs)
              fixed_w_mat[child, ] <- fixed_ccfs
            }
          }
        }
      } else {
        for (child in as.numeric(temp_curr_edges$child)) {
          child_ccfs <- w_mat[child, ]
          fixed_ccfs <- mapply(function(child_ccf, parent_ccf) ifelse(child_ccf > parent_ccf, parent_ccf, child_ccf),
                               child_ccf = child_ccfs,
                               parent_ccf = parent_ccfs)
          fixed_w_mat[child, ] <- fixed_ccfs
        }
      }
      
      
      child_next_edges <- edges_left %>% 
        filter(parent %in% temp_curr_edges$child)
      next_edges <- bind_rows(next_edges,
                              child_next_edges)
    }
    curr_edges <- next_edges
  }
  return(fixed_w_mat)
}
