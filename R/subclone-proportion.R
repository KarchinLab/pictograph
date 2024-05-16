#' Calculate proportions of subclones in each sample (assumes CCFs comply with lineage precedence and sum condition)
#'
#' @param w_mat Matrix of CCF estimates (from \code{estimateCCFs})
#' @param tree_edges Tibble of tree edges with columns edge, parent, and child
#' @export
calcSubcloneProportions <- function(w_mat, tree_edges) {
  K <- nrow(w_mat)
  S <- ncol(w_mat)
  subclone_props <- matrix(NA, nrow = K, ncol = S)
  
  for (i in seq_len(nrow(w_mat))) {
    children <- tree_edges %>%
      filter(parent == as.character(i)) %>%
      pull(child) %>%
      as.numeric()
    
    if (length(children) == 1) {
      children_ccfs <- w_mat[children, ]
    } else if (length(children) > 1) {
      children_ccfs <- w_mat[children, ,drop=FALSE] %>%
        colSums
    } else {
      children_ccfs <- rep(0, ncol(w_mat))
    }
    
    subclone_props[i, ] <- w_mat[i, ] - children_ccfs
  }
  
  # normalize subclone_props matrix so the props add up to 1
  subclone_props[subclone_props < 0] = 0
  subclone_props = round(t(t(subclone_props) / colSums(subclone_props)),digit = 3)
  
  return(subclone_props)
}

#' Plot pie charts for subclone proportions in each sample
#'
#' @param subclone_props matrix of subclone proportions (returned from \code{calcSubcloneProportions})
#' @param sample_names (Optional) Vector of sample names. Should be in the order of columns of subclone_props
#' @export
plotSubclonePie <- function(subclone_props, palette=viridis::viridis, sample_names = NULL, title_size=16, legend_size=10) {
  if (is.null(sample_names)) sample_names <- paste0("Sample ", 1:ncol(subclone_props))
  props_tb <- subclone_props %>%
    magrittr::set_colnames(sample_names) %>%
    as_tibble() %>%
    mutate(Subclone = factor(paste0("Clone ", 1:nrow(.)),
                             levels = paste0("Clone ", 1:nrow(.)))) %>%
    pivot_longer(cols = sample_names,
                 names_to = "Sample",
                 values_to = "Proportion")
  
  clone_colors <- palette(nrow(subclone_props))
  ggplot(props_tb, aes(x="", y=Proportion, fill = Subclone)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = clone_colors, drop = F) +
    theme_void() +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size=legend_size), legend.title = element_text(size=legend_size)) + 
    facet_wrap(~Sample) +
    theme(strip.text.x = element_text(size=title_size))
  
}

#' Plot subclone proportions in each sample as stacked bar chart
#'
#' @param subclone_props matrix of subclone proportions (returned from \code{calcSubcloneProportions})
#' @param sample_names (Optional) Vector of sample names. Should be in the order of columns of subclone_props
#' @param label_cluster (Default FALSE) Whether to add cluster label to text on stacked bar
#' @export
plotSubcloneBar <- function(subclone_props, palette=viridis::viridis, sample_names = NULL, label_cluster = FALSE) {
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample ", 1:ncol(subclone_props))
  }
  
  clone_colors <- palette(nrow(subclone_props))
  color_half <- nrow(subclone_props) / 2
  color_half_vec <- factor(ifelse(1:nrow(subclone_props) < color_half, "white", "black"),
                           c("white", "black"))
  
  props_tb <- subclone_props %>%
    magrittr::set_colnames(sample_names) %>%
    as_tibble %>%
    mutate(Clone = as.factor(1:nrow(.)),
           text_color = color_half_vec) %>%
    pivot_longer(cols = all_of(sample_names),
                 names_to = "Sample",
                 values_to = "Proportion of Tumor")
  if (label_cluster) {
    text_label <- paste0("Clone ", props_tb$Clone, ": ", props_tb$`Proportion of Tumor`)
    props_tb <- props_tb %>%
      mutate(text_label = text_label)
  } else {
    props_tb <- props_tb %>%
      mutate(text_label = `Proportion of Tumor`)
  }
  
  stacked_bar <- ggplot(props_tb, aes(x = Sample, y = `Proportion of Tumor`, fill = Clone)) +
    theme_light() +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = clone_colors) +
    geom_text(data = subset(props_tb, `Proportion of Tumor` != 0),
              aes(label = text_label, colour = text_color),
              size = 4, position = position_stack(vjust = 0.5)) +
    xlab("") +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    scale_color_manual(values = c("white", "black"), guide = "none") +
    ylim(0,1.05) + 
    theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
    theme(legend.text = element_text(size=12), legend.title = element_text(size=12))
  
  return(stacked_bar)
}