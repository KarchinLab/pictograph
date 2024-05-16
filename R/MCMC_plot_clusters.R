#' Plot probabilities of mutation cluster assignments - vertical
#' 
#' @export
#' @import ggplot2
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @param z_chain MCMC chain of mutation cluster assignment values, which is the second item in the list returned by \code{clusterSep}
#' @param mcf_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
#' @param filter_thresh Lowest posterior probability to include cluster assignment. Default value is 0.05 (inclusive)
#' @param MutID (Optional) Vector of mutation IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Mut_ID)
#' @param SampleID (Optional) Vector of sample IDs for labeling purposes. Same order as supplied as input data (e.g. indata$Sample_ID)
plotClusterAssignmentProbVertical <- function(z_chain, 
                                              mcf_chain,
                                              filter_thresh = 0.05,
                                              MutID = NULL,
                                              SampleID = NULL) {
  
  mcmc_z <- generateZPostSummary(z_chain, mcf_chain, filter_thresh, MutID, SampleID)
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
                                 mcf_chain,
                                 filter_thresh = 0.05,
                                 MutID = NULL,
                                 SampleID = NULL) {
  
  map_z <- estimateClusterAssignments(z_chain)
  map_mcf <- estimateMCFs(mcf_chain)
  
  I <- length(unique(z_chain$Parameter))
  K <- max(unique(z_chain$value))
  num_iter <- max(z_chain$Iteration)
  S <- ncol(map_mcf)
  
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
  
  
  
  tiers <- generateTiers(map_mcf, sample_labels)
  
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

#' Plot CCF chain trace 
#'
#' @export
#' @param mcf_chain MCMC chain of CCF values, which is the first item in the list returned by \code{mergeSetChains}
plotChainsMCF <- function(mcf_chain) {
  cluster <- strsplit(as.character(mcf_chain$Parameter), ",") %>%
    sapply(., function(x) gsub("mcf\\[", "", x[1])) %>%
    as.numeric
  sample <- strsplit(as.character(mcf_chain$Parameter), ",") %>%
    sapply(., function(x) gsub("\\]", "", x[2])) %>%
    as.numeric
  
  mcf_chain <- mcf_chain %>%
    mutate(Cluster = factor(paste0("Cluster ", cluster), 
                            levels = paste0("Cluster ", sort(unique(cluster)))), 
           Sample = factor(paste0("Sample ", sample),
                           levels = paste0("Sample ", sort(unique(sample)))))
  
  ggplot(mcf_chain, aes(x = Iteration, y = value)) +
    geom_line() +
    theme_light() +
    facet_grid(Cluster ~ Sample) +
    ylab("Cancer Cell Fraction")
}

#' Plot cluster CCF posterior distributions as violin plots
#' 
#' @export
#' @param mcf_chain MCMC chain of CCF values
#' @param z_chain (Optional) MCMC chain of mutation cluster assignment values. If provided, cluster names will show the number of mutations assigned in brackets
#' @param indata (Optional) List of input data items 
plotMCFViolin <- function(mcf_chain, z_chain = NULL, indata = NULL) {
  # process data
  vdat <- violinProcessData(mcf_chain, indata)
  
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

violinProcessData <- function(mcf_chain, indata = NULL) {
  mcf_mat <- estimateMCFs(mcf_chain)
  est_K <- nrow(mcf_mat)
  
  if (is.null(indata$Sample_ID)) {
    sample_names <- paste0("Sample ", 1:ncol(mcf_mat))
  } else {
    sample_names <- indata$Sample_ID
  }
  
  vdat <- mcf_chain %>%
    mutate(sample=stringr::str_replace_all(Parameter, "mcf\\[[:digit:]+,", ""),
           sample=stringr::str_replace_all(sample, "\\]", ""),
           cluster=stringr::str_replace_all(Parameter, "mcf\\[", ""),
           cluster=stringr::str_replace_all(cluster, ",[:digit:]\\]", "")) %>%
    mutate(sample=as.numeric(sample),
           sample=sample_names[sample],
           sample=factor(sample, sample_names),
           cluster=as.numeric(cluster),
           cluster=paste0("Cluster ", cluster),
           cluster=factor(cluster, level=paste("Cluster", 1:est_K)))
  
  tiers <- generateTiers(mcf_mat, sample_names)
  
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