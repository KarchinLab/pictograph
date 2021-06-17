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
#' @import tidyverse
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
  z.plot
}

#' Plot cluster CCF posterior distributions
#' 
#' @export
#' @import tidyverse
#' @param w_chain MCMC chain of CCF values, which is the first item in the list returned by \code{clusterSep}
plotDensityCCF <- function(w_chain) {
  w_chain <- w_chain %>% 
    mutate(Cluster = as.numeric(gsub("w\\[", "", sapply(w_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][1])))) %>%
    mutate(Sample = gsub("\\]", "", sapply(w_chain$Parameter, function(x) strsplit(as.character(x), ",")[[1]][2])))
  K <- max(as.numeric(w_chain$Cluster))
  S <- max(as.numeric(w_chain$Sample))
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
}

#' Plot single tree 
#' 
#' @export
#' @param edges tibble of edges with columns edge, parent, child
#' @import tidyverse
plotTree <- function(edges) {
  plotGraph(edgesToAmLong(edges))
}

#' Plot ensemble tree
#' 
#' @export
#' @param trees list of tibbles of edges, each with columns edge, parent, child
#' @import tidyverse
plotEnsembleTree <- function(trees) {
  am_chain <- lapply(trees, edgesToAmLong)
  post_am <- getPosteriorAmLong(am_chain)
  plotPosteriorAmLong(post_am)
}