plot.chain.trace <- function(chain, K) {
  ggplot(chain, aes(x = Iteration, y = value)) +
    geom_line() +
    facet_wrap(~Parameter, nrow = K) +
    theme_light()
}

get.parameter.chain <- function(param, chains) {
  chains[grep(paste0(param, "\\["), chains$Parameter), ]
}

#' Plot probabilities of mutation cluster assignments
#' 
#' @export
#' @import ggplot2
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
plotClusterAssignmentProb <- function(z_chain) {
  I <- length(unique(z_chain$Parameter))
  K <- max(unique(z_chain$value))
  num_iter <- max(z_chain$Iteration)
  
  mcmc_z <- z_chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              num_iter=num_iter) %>%
    mutate(probability=n/num_iter) %>%
    ungroup()
  
  z.seg.tb <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(z1 = min(value), z2 = max(value))
  
  z.plot <- ggplot(mcmc_z, aes(x = Parameter, y = value, color = probability)) +
    theme_light() +
    scale_x_discrete(labels = 1:I, name = "Variant") +
    scale_y_continuous(breaks = 1:K, name = "Cluster") +
    geom_segment(data = z.seg.tb, 
                 aes(x=Parameter, xend=Parameter,
                     y=z1, yend=z2),
                 color="black", linetype=2) +
    geom_point() +
    theme(panel.grid.minor = element_blank())
  return(z.plot)
}

summarizeZPost <- function(z_chain) {
  I <- length(unique(z_chain$Parameter))
  K <- max(unique(z_chain$value))
  num_iter <- max(z_chain$Iteration)
  mcmc_z <- z_chain %>%
    group_by(Parameter, value) %>%
    summarize(n=n(),
              num_iter=num_iter) %>%
    mutate(probability=n/num_iter) %>%
    ungroup()
  return(mcmc_z)
}

#' Plot probabilities of mutation cluster assignments - vertical
#' 
#' @export
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @param filter_thresh Lowest posterior probability to include cluster assignment. Default value is 0.05 (inclusive)
#' @param MutID (Optional) Vector of mutation IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Mut_ID)
#' @param SampleID (Optional) Vector of sample IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Sample_ID)
plotClusterAssignmentProbVertical <- function(z_chain, 
                                              w_chain,
                                              filter_thresh = 0.05,
                                              MutID = NULL,
                                              SampleID = NULL) {
  
  mcmc_z <- generateZPostSummary(z_chain, w_chain, filter_thresh, MutID, SampleID)
  K <- max(mcmc_z$value)
  z.seg.tb <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(z1 = min(value), z2 = max(value)) %>%
    ungroup() %>%
    mutate(Variant = as.numeric(gsub("z\\[|]", "", Parameter)),
           Mut_ID =  mcmc_z$Mut_ID[match(Variant, mcmc_z$Variant)],
           Sample_presence =  mcmc_z$Sample_presence[match(Variant, mcmc_z$Variant)])
  
  z.plot <- ggplot(mcmc_z, aes(x = value, y = Mut_ID, color = probability)) +
    theme_light() +
    scale_y_discrete(drop = T, name = "Variant") +
    scale_x_continuous(breaks = 1:K, name = "Cluster", labels = 1:K) +
    geom_segment(data = z.seg.tb, 
                 aes(y=Mut_ID, yend=Mut_ID,
                     x=z1, xend=z2),
                 color="black", linetype=2) +
    geom_point() +
    theme(panel.grid.minor = element_blank(),
          strip.background=element_blank(),
          strip.text = element_text(colour = 'black'),
          strip.text.y = element_text(angle = 0)) +
    facet_grid(Sample_presence~., scales = "free", space = "free") +
    scale_color_gradient(limits = c(0,1))
  return(z.plot)
}

generateZPostSummary <- function(z_chain, 
                                 w_chain,
                                 filter_thresh = 0.05,
                                 MutID = NULL,
                                 SampleID = NULL) {
  
  map_z <- estimateClusterAssignments(z_chain)
  map_w <- estimateCCFs(w_chain)
  
  I <- length(unique(z_chain$Parameter))
  K <- max(unique(z_chain$value))
  num_iter <- max(z_chain$Iteration)
  S <- ncol(map_w)
  
  if (is.null(MutID)) {
    mut_labels <- 1:I
  } else {
    mut_labels <- MutID
  }
  if (is.null(SampleID)) {
    sample_labels <- paste0("Sample ", 1:S)
  } else {
    sample_labels <- SampleID
  }
  

  
  tiers <- generateTiers(map_w, sample_labels)
  
  mcmc_z <- summarizeZPost(z_chain) %>%
    filter(probability >= filter_thresh)
  
  # Variant sample presence 
  var_sample_pres <- map_z %>%
    ungroup() %>%
    mutate(cluster_num = value,
           Variant = as.numeric(gsub("z\\[|]", "", Parameter)),
           Mut_ID = mut_labels[Variant],
           Sample_presence = tiers$samples[value])
  # sample presence order 
  sample_pres_order <- tiers %>% 
    select(samples, tier) %>% 
    distinct() %>% 
    arrange(-tier) %>% 
    pull(samples)
  # variant order 
  var_order <- map_z %>%
    arrange(-value) %>%
    mutate(Variant = as.numeric(gsub("z\\[|]", "", Parameter)),
           Mut_ID = mut_labels[Variant]) %>%
    pull(Mut_ID)
  
  mcmc_z <- mcmc_z %>%
    mutate(Variant = as.numeric(gsub("z\\[|]", "", Parameter)),
           Mut_ID = factor(mut_labels[Variant], var_order),
           Sample_presence = factor(var_sample_pres$Sample_presence[Variant],
                                    sample_pres_order))
  return(mcmc_z)
}

grabBIC <- function(all_set_results) {
  bic_list <- lapply(all_set_results, function(x) x$BIC) 
  bic_list_filt <- bic_list[sapply(bic_list, function(x) typeof(x) != "logical")]
  bic_tb <- bic_list_filt %>%
    bind_rows(.id = 'Set') %>%
    mutate(Set = factor(Set, levels = names(all_set_results)))
  return(bic_tb)
}

grabBestK <- function(all_set_results) {
  best_k <- lapply(all_set_results, function(x) x$best_K) %>%
    stack() %>%
    rename(Set = ind,
           Best_K = values) %>%
    as_tibble()
  return(best_k)
}

#' Plot probabilities of mutation cluster assignments (vertical) for tested K across all mutation sets
#' 
#' @export
#' @import ggplot2
#' @import dplyr
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param outdir Path to directory for output of plots 
#' @param SampleID (Optional) Vector of sample IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Sample_ID)
#' @param filter_thresh Lowest posterior probability to include cluster assignment. Default value is 0.05 (inclusive)
#' @param compare Option to only plot cluster probabilities for K chosen by minimum BIC, elbow or knee of plot when different (default FALSE plots all K tested)
plotAllZProb <- function(all_set_results, outdir, SampleID = NULL, filter_thresh = 0.05, compare = FALSE) {
  if (is.null(SampleID)) {
    S <- estimateCCFs(all_set_results[[1]]$all_chains[[1]]$w_chain) %>% 
      ncol
    sample_names <- paste0("Sample ", 1:S)
  } else {
    sample_names <- SampleID
  }
  
  num_sets <- length(all_set_results)
  set_names_bin <- names(all_set_results)
  
  if (compare) k_tb <- writeSetKTable(all_set_results)
  
  for (set in set_names_bin) {
    set_name_full <- sample_names[as.logical(as.numeric(strsplit(set, "")[[1]]))] %>%
      paste0(., collapse = ",\n")
    
    all_set_chains <- all_set_results[[set]]$all_chains
    k_tested <- as.numeric(gsub("K", "", names(all_set_chains)))
    S <- all_set_chains[[1]]$w_chain %>%
      estimateCCFs %>%
      ncol
    I <- all_set_chains[[1]]$z_chain %>%
      estimateClusterAssignments %>%
      nrow
    num_iter <- all_set_chains[[1]]$z_chain %>%
      pull(Iteration) %>%
      max
    
    if (I == 1) next # plot is trivial if only one mutation in set; skip plotting
    
    if (compare) {
      # only plot for K = minimum BIC and K = elbow
      k_to_plot <- k_tb %>% 
        filter(set_name_bin == set) %>% 
        select(min_BIC, elbow, knee) %>% 
        unlist %>% 
        unname %>% 
        unique
      if (any(is.na(k_to_plot)))  k_to_plot <- k_to_plot[-which(is.na(k_to_plot))]
     
    } else {
      k_to_plot <- k_tested
    }

    for (k in k_to_plot) {
      z_plot_file <- file.path(outdir, paste0(set, "_", names(all_set_chains)[k], "_z_plot.pdf"))
      plot_title <- paste0(set_name_full, ": ", names(all_set_chains)[k])
      temp_chains <- all_set_chains[[k]]
      
      mcmc_z <- summarizeZPost(temp_chains$z_chain) %>%
        filter(probability >= filter_thresh) %>%
        mutate(Variant = as.numeric(gsub("z\\[|]", "", Parameter)))
      
      # order varints by highest probability cluster assignment
      var_order <-  mcmc_z %>%
        group_by(Variant) %>%
        summarize(map_z = value[which.max(probability)]) %>%
        arrange(-map_z) %>%
        pull(Variant)
      
      mcmc_z <- mcmc_z %>%
        mutate(Variant = factor(Variant, levels = var_order))
      
      # segments connecting multiple assignments for each variant
      z.seg.tb <- mcmc_z %>%
        group_by(Variant) %>%
        summarize(z1 = min(value), z2 = max(value)) %>%
        ungroup() 
      
      z.plot <- ggplot(mcmc_z, aes(x = value, y = Variant, color = probability)) +
        theme_light() +
        scale_y_discrete(drop = T, name = "Variant", labels = NULL) +
        scale_x_continuous(breaks = 1:k, name = "Cluster", labels = 1:k) +
        geom_segment(data = z.seg.tb, 
                     aes(y=Variant, yend=Variant,
                         x=z1, xend=z2),
                     color="black", linetype=2) +
        geom_point() +
        theme(panel.grid.minor = element_blank(),
              strip.background=element_blank(),
              strip.text = element_text(colour = 'black'),
              strip.text.y = element_text(angle = 0)) +
        scale_color_gradient(limits = c(0,1)) +
        ggtitle(plot_title)
      plot_height <- max(3, I/15) + S/2
      ggsave(z_plot_file, plot = z.plot, height = plot_height, width = max(3, k))
    }
  }
}

#' Plot BIC for all mutation sets 
#' 
#' @export
#' @import dplyr
#' @import ggplot2
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param SampleID (Optional) Vector of sample IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Sample_ID)
plotBIC <- function(all_set_results, Sample_ID = NULL) {
  bic_tb <- grabBIC(all_set_results) 
  k_tb <- writeSetKTable(all_set_results)
  min_BIC <- grabBestK(all_set_results) %>%
    rename(K_tested = Best_K) %>%
    left_join(., bic_tb,
              by = c("Set", "K_tested")) %>%
    mutate(Choice = "Minimum")
  elbow <- k_tb %>%
    select(K_tested = elbow, Set = set_name_bin) %>%
    mutate(Set = factor(Set, levels(bic_tb$Set))) %>%
    left_join(., bic_tb,
              by = c("Set", "K_tested")) %>%
    mutate(Choice = "Elbow")
  knee <- k_tb %>%
    select(K_tested = knee, Set = set_name_bin) %>%
    mutate(Set = factor(Set, levels(bic_tb$Set))) %>%
    left_join(., bic_tb,
              by = c("Set", "K_tested")) %>%
    mutate(Choice = "Knee")
  
  if (any(is.na(min_BIC$BIC))) min_BIC <- min_BIC[-which(is.na(min_BIC$BIC)), ]
  if (any(is.na(elbow$BIC))) elbow <- elbow[-which(is.na(elbow$BIC)), ]
  if (any(is.na(knee$BIC))) knee <- knee[-which(is.na(knee$BIC)), ]
  
  # rename sets
  if (is.null(Sample_ID)) {
    S <- all_set_results[[1]]$best_chains$w_chain %>% 
      estimateCCFs %>% 
      ncol
    Sample_ID <- paste0("Sample ", 1:S)
  }
  
  set_name_tb <- tibble(set_bin = unique(bic_tb$Set),
                        set_name = sapply(unique(bic_tb$Set), function(x) getSetName(x, Sample_ID))) %>%
    mutate(set_name = factor(set_name, set_name))
   
  bic_tb <- bic_tb %>%
    mutate(Set_name = sapply(bic_tb$Set, function(x) getSetName(x, Sample_ID))) %>%
    mutate(Set_name = factor(Set_name, set_name_tb$set_name))
  k_choices <- bind_rows(min_BIC, elbow, knee) %>%
    mutate(Set_name = sapply(Set, function(x) getSetName(x, Sample_ID))) %>%
    mutate(Set_name = factor(Set_name, set_name_tb$set_name)) %>%
    mutate(Choice = factor(Choice, c("Minimum", "Elbow", "Knee")))
  
  # plot
  bic_plot <- ggplot(bic_tb, aes(x = K_tested, y = BIC)) +
    theme_light() +
    facet_wrap(~Set_name, scales = "free") +
    geom_line() +
    geom_point(data = k_choices, aes(color = Choice, size = Choice, shape = Choice), stroke = 1) +
    scale_shape_manual(values=c(19, 1, 1))+
    scale_color_manual(values=c('#999999',"#56B1F7", "#E69F00")) +
    scale_size_manual(values=c(1, 3, 5)) +
    #geom_point(data = elbow, color = , size = 4, shape = 1)
    theme(strip.background=element_blank(),
          strip.text = element_text(colour = 'black'),
          panel.grid.minor = element_blank()) +
    xlab("K")
  
  return(bic_plot)
}

#' Convert binary set names to long form with Sample_ID
getSetName <- function(binary_name, Sample_ID, collapse_string = ", \n") {
  split_bin <- strsplit(as.character(binary_name), "") %>%
    .[[1]] %>%
    as.numeric() %>%
    as.logical()
  samples_present <- Sample_ID[split_bin]
  set_name <- paste0(samples_present, collapse = collapse_string)
  return(set_name)
}

#' Plot cluster CCF posterior distributions
#' 
#' @export
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
plotDensityCCF <- function(w_chain) {
  w_chain <- w_chain %>% 
    mutate(Cluster = as.numeric(gsub("w\\[", "", sapply(w_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(Sample = gsub("\\]", "", sapply(w_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][2])))
  K <- max(as.numeric(w_chain$Cluster))
  S <- max(as.numeric(w_chain$Sample))
  suppressWarnings(print(
    ggplot(w_chain, aes(x = value)) +
      geom_density() +
      facet_wrap(~Cluster + Sample, ncol = S, scales = "free_y") +
      theme_light() +
      ylab("Density") + xlab("CCF") + 
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            strip.background = element_blank(),
            #strip.text = element_text(colour = 'black'),
            strip.text.x = element_blank(),
            #axis.text.x=element_blank(),
            #axis.ticks.x=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank()) +
      scale_x_continuous(breaks=c(0, 0.5, 1))
  ))
}

#' Plot cluster CCF posterior distributions as violin plots
#' 
#' @export
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @importFrom stringr str_replace_all
#' @param w_chain MCMC chain of CCF values
#' @param z_chain (Optional) MCMC chain of mutation cluster assignment values. If provided, cluster names will show the number of mutations assigned in brackets
#' @param indata (Optional) List of input data items 
plotCCFViolin <- function(w_chain, z_chain = NULL, indata = NULL) {
  # process data
  vdat <- violinProcessData(w_chain, indata)
  
  if (!is.null(z_chain)) {
    num_muts_in_clusters <- estimateClusterAssignments(z_chain) %>%
      group_by(value) %>%
      summarize(num_muts = n()) %>%
      ungroup() %>%
      rename(cluster = value)
    num_muts <- num_muts_in_clusters$num_muts[match(vdat$cluster_num, num_muts_in_clusters$cluster)]
    new_cluster_labels <- paste0("Cluster ", 
                                 vdat$cluster_num, 
                                 " [", num_muts,"]")
    vdat <- vdat %>%
      mutate(cluster = factor(new_cluster_labels, unique(new_cluster_labels)))
  }

  # plot violins
  vplot <- plotViolin(vdat)
  return(vplot)
}

generateTiers <- function(w_mat, Sample_ID) {
  clusters <- paste0("Cluster ", seq_len(nrow(w_mat)))
  bin <- w_mat > 0
  samples <- apply(bin, 1, function(x) paste(Sample_ID[x], collapse = ",\n"))
  tier <- rowSums(bin)
  tiers <- tibble(cluster = clusters,
                  cluster_num = seq_len(nrow(w_mat)),
                  samples = samples,
                  tier = tier)
  return(tiers)
}

violinProcessData <- function(w_chain, indata = NULL) {
  w_mat <- estimateCCFs(w_chain)
  est_K <- nrow(w_mat)
  
  if (is.null(indata$Sample_ID)) {
    sample_names <- paste0("Sample ", 1:ncol(w_mat))
  } else {
    sample_names <- indata$Sample_ID
  }
  
  vdat <- w_chain %>%
    mutate(sample=stringr::str_replace_all(Parameter, "w\\[[:digit:]+,", ""),
           sample=stringr::str_replace_all(sample, "\\]", ""),
           cluster=stringr::str_replace_all(Parameter, "w\\[", ""),
           cluster=stringr::str_replace_all(cluster, ",[:digit:]\\]", "")) %>%
    mutate(sample=as.numeric(sample),
           sample=sample_names[sample],
           sample=factor(sample, sample_names),
           cluster=as.numeric(cluster),
           cluster=paste0("Cluster ", cluster),
           cluster=factor(cluster, level=paste("Cluster", 1:est_K)))
  
  tiers <- generateTiers(w_mat, sample_names)
  
  vdat <- vdat %>%
    mutate(cluster=as.character(cluster)) %>%
    left_join(tiers, by="cluster") %>%
    mutate(cluster=factor(cluster, tiers$cluster),
           tier=factor(tier, sort(unique(tiers$tier))))
  return(vdat)
}

plotViolin <- function(vdat) {
  vplot <- ggplot(vdat, aes(sample, value)) +
    geom_violin(aes(fill=tier),
                alpha=0.6,
                scale="width",
                draw_quantiles=c(0.25, 0.5, 0.75),
                color="white") +
    geom_violin(fill="transparent", color="black",
                scale="width", draw_quantiles=0.5) +
    theme_bw(base_size=12) +
    theme(strip.background=element_blank(),
          axis.text.x=element_text(size=12),
          panel.grid=element_blank(),
          legend.pos="bottom") +
    facet_wrap(~cluster, nrow=1) +
    ylab("Posterior CCF") + xlab("") +ylim(c(0, 1)) +
    guides(fill=guide_legend("Sample-presence")) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  return(vplot)
}

#' Plot single tree 
#' 
#' @export
#' @param edges tibble of edges with columns edge, parent, child
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
plotTree <- function(edges, palette=viridis::viridis) {
  plotGraph(edgesToAmLong(edges), colorScheme(edges, palette))
}

#' generate colors for each vertice
#' @export
colorScheme <- function(edges, palette=viridis::viridis) {
  v_sorted = sort(unique(c(edges$parent, edges$child)))
  v_sorted = c(sort(as.integer(v_sorted[!v_sorted=='root'])), "root")
  # root_idx <- which(v_sorted=="root")
  colors <- c(palette(length(v_sorted)-1), "white")
  v_color <- tibble(v_sorted, colors)
  return(v_color)
}

#' Plot ensemble tree
#' 
#' @export
#' @param trees list of tibbles of edges, each with columns edge, parent, child
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
plotEnsembleTree <- function(trees, palette=viridis::viridis) {
  am_chain <- lapply(trees, edgesToAmLong)
  post_am <- getPosteriorAmLong(am_chain)
  plotPosteriorAmLong(post_am, colorScheme(trees[[1]], palette))
}

#' Plot CCF chain trace 
#'
#' @export
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{mergeSetChains}
plotChainsCCF <- function(w_chain) {
  cluster <- strsplit(as.character(w_chain$Parameter), ",") %>%
    sapply(., function(x) gsub("w\\[", "", x[1])) %>%
    as.numeric
  sample <- strsplit(as.character(w_chain$Parameter), ",") %>%
    sapply(., function(x) gsub("\\]", "", x[2])) %>%
    as.numeric
  
  w_chain <- w_chain %>%
    mutate(Cluster = factor(paste0("Cluster ", cluster), 
                            levels = paste0("Cluster ", sort(unique(cluster)))), 
           Sample = factor(paste0("Sample ", sample),
                           levels = paste0("Sample ", sort(unique(sample)))))
  
  ggplot(w_chain, aes(x = Iteration, y = value)) +
    geom_line() +
    theme_light() +
    facet_grid(Cluster ~ Sample) +
    ylab("Cancer Cell Fraction")
}

#' Plot posterior predictive distribution for number of variant reads
#' 
#' @export
#' @import ggplot2
#' @import tibble
#' @import tidyr
#' @param ystar_chain MCMC chain of ystar values, which is the third item in the list returned by \code{clusterSep}
#' @param indata List of input data items 
#' @param SampleID (Optional) Vector of sample IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Sample_ID). If not provided, function will use the Sample_ID in indata
#' @param Mutation_ID (Optional) Vector of mutation IDs. If not provided, function will use the MutID in indata
plotPPD <- function(ystar_chain, indata, 
                    Sample_ID = NULL,
                    Mutation_ID = NULL) {
  I <- indata$I
  if (is.null(Sample_ID)) Sample_ID <- indata$Sample_ID
  if (is.null(Mutation_ID)) Mutation_ID <- indata$Mut_ID
  
  ppd.summaries <- ystar_chain %>%
    group_by(Parameter) %>%
    summarize(mean=mean(value),
              median=median(value),
              q1=quantile(value, 0.025),
              q3=quantile(value, 0.975))
  
  observed_y <- indata$y %>%
    magrittr::set_colnames(Sample_ID) %>%
    as_tibble() %>%
    mutate(Mutation_index = 1:I) %>%
    pivot_longer(cols = Sample_ID,
                 names_to = "Sample",
                 values_to = "observed_y") %>%
    mutate(s = match(Sample, Sample_ID),
           Parameter = paste0("ystar[", Mutation_index, ",", s, "]"))
  
  ppd.summaries2 <- ppd.summaries %>%
    left_join(., observed_y, by = "Parameter")
  points <- ppd.summaries2 %>%
    select(Parameter, Mutation_index, Sample, observed_y, median) %>%
    rename("Observed variant read count" = observed_y,
           "Posterior median" = median) %>%
    pivot_longer(cols = c("Observed variant read count", "Posterior median"),
                 names_to = "type",
                 values_to = "value") %>%
    left_join(., ppd.summaries2, by = c("Parameter",
                                        "Mutation_index",
                                        "Sample"))
  points_order <- ppd.summaries2 %>%
    filter(Sample == Sample_ID[1]) %>%
    arrange(observed_y) %>%
    pull(Mutation_index)
  
  variant_gene_names <- Mutation_ID %>%
    sapply(., function(x) strsplit(x, "_")[[1]][1]) %>%
    as.character()
  variant_gene_names_ordered <- variant_gene_names[points_order]
  points2 <- points %>%
    mutate(Mutation_index = factor(Mutation_index,
                                   levels = points_order))
  
  # plot
  points3 <- points2 %>%
    mutate(type=gsub("Observed variant read count",
                     "Observed variant\nread count",
                     type))
  colors <- setNames(c("gray40", "steelblue"), unique(points3$type))
  fill <- setNames(c("white", "steelblue"), unique(points3$type))
  points3 %>%
    ggplot(aes(x=mean, y=Mutation_index,
               xmin=q1,
               xmax=q3)) +
    geom_errorbar(color="gray") +
    geom_point(aes(x = value, y = Mutation_index, color=type, fill=type),
               size=2, pch=21) +
    scale_y_discrete(labels = variant_gene_names_ordered,
                     breaks = points_order) +
    theme_bw(base_size=15) +
    theme(#axis.text.y=element_blank(),
      axis.title.x=element_text(size=20),
      axis.text.x=element_text(size=17),
      strip.text=element_text(size=22),
      panel.grid=element_blank(),
      axis.ticks.y=element_blank(),
      legend.text=element_text(size=18),
      panel.background=element_rect(fill="white",
                                    color="black"),
      legend.pos="bottom",
      strip.background=element_blank()) +
    scale_color_manual(name="",
                       labels=names(colors),
                       values=colors) +
    scale_fill_manual(name="",
                      labels=names(colors),
                      values=fill) +
    xlab("Variant allele count") +
    ylab("") +
    facet_wrap(~Sample) +
    guides(color=guide_legend(override.aes=list(size=3)))
}
