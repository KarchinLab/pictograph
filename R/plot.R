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
           Variant = 1:I,
           Mut_ID = mut_labels,
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
    mutate(Variant = as.numeric(Parameter),
           Mut_ID = mut_labels[Variant]) %>%
    pull(Mut_ID)
  
  z.seg.tb <- mcmc_z %>%
    group_by(Parameter) %>%
    summarize(z1 = min(value), z2 = max(value)) %>%
    ungroup() %>%
    mutate(Variant = 1:I,
           Mut_ID = factor(mut_labels, var_order),
           Sample_presence = factor(var_sample_pres$Sample_presence, sample_pres_order))
  
  mcmc_z <- mcmc_z %>%
    mutate(Variant = as.numeric(Parameter),
           Mut_ID = factor(mut_labels[Variant], var_order),
           Sample_presence = factor(var_sample_pres$Sample_presence[Variant],
                                    sample_pres_order))
  
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

#' Plot BIC for all mutation sets 
#' 
#' @export
#' @import dplyr
#' @import ggplot2
#' @param all_set_results List of MCMC results for each mutation set; returned by \code{clusterSep}
#' @param SampleID (Optional) Vector of sample IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Sample_ID)
plotBIC <- function(all_set_results, Sample_ID = NULL) {
  bic_tb <- grabBIC(all_set_results) 
  min_BIC <- grabBestK(all_set_results) %>%
    rename(K_tested = Best_K) %>%
    right_join(., bic_tb,
              by = c("Set", "K_tested"))
  
  # rename sets
  if (is.null(Sample_ID)) {
    S <- all_set_results[[1]]$best_chains$w_chain %>% 
      estimateCCFs %>% 
      ncol
    Sample_ID <- paste0("Sample ", 1:S)
  }
  
  binary_set_names <- unique(bic_tb$Set)
  set_names <- sapply(binary_set_names, function(x) getSetName(x, Sample_ID))
  
  bic_tb <- bic_tb %>%
    mutate(Set_name = factor(set_names[Set], levels = set_names))
  min_BIC <- min_BIC %>%
    mutate(Set_name = factor(set_names[Set], levels = set_names))
  
  # plot
  bic_plot <- ggplot(bic_tb, aes(x = K_tested, y = BIC)) +
    theme_light() +
    facet_wrap(~Set_name, scales = "free") +
    geom_line() +
    geom_point(data = min_BIC, color = "#56B1F7", size = 2) +
    theme(strip.background=element_blank(),
          strip.text = element_text(colour = 'black'),
          panel.grid.minor = element_blank()) +
    xlab("K")
  
  return(bic_plot)
}

#' Convert binary set names to long form with Sample_ID
getSetName <- function(binary_name, Sample_ID) {
  split_bin <- strsplit(as.character(binary_name), "") %>%
    .[[1]] %>%
    as.numeric() %>%
    as.logical()
  samples_present <- Sample_ID[split_bin]
  set_name <- paste0(samples_present, collapse = ", \n")
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
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @param indata List of input data items 
plotCCFViolin <- function(w_chain, indata) {
  # process data
  vdat <- violinProcessData(w_chain, indata)

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

violinProcessData <- function(w_chain, indata) {
  w_mat <- estimateCCFs(w_chain)
  est_K <- nrow(w_mat)
  
  if (is.null(indata$Sample_ID)) {
    sample_names <- paste0("Sample ", 1:indata$S)
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
    theme_bw(base_size=20) +
    theme(strip.background=element_blank(),
          axis.text.x=element_text(size=14),
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
plotTree <- function(edges) {
  plotGraph(edgesToAmLong(edges))
}

#' Plot ensemble tree
#' 
#' @export
#' @param trees list of tibbles of edges, each with columns edge, parent, child
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
plotEnsembleTree <- function(trees) {
  am_chain <- lapply(trees, edgesToAmLong)
  post_am <- getPosteriorAmLong(am_chain)
  plotPosteriorAmLong(post_am)
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
