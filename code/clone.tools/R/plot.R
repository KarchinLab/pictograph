plot.chain.trace <- function(chain, K) {
  ggplot(chain, aes(x = Iteration, y = value)) +
    geom_line() +
    facet_wrap(~Parameter, nrow = K) +
    theme_light()
}

get.parameter.chain <- function(param, chains) {
  chains[grep(paste0(param, "\\["), chains$Parameter), ]
}

plot.variant.z <- function(z_chain) {
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
    geom_point()
}