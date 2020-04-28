#!/usr/bin/env Rscript
suppressMessages(library(rjags))
suppressMessages(library(coda))
suppressMessages(library(ggmcmc))
suppressMessages(library(optparse))
suppressMessages(library(gridExtra))
suppressMessages(library(stringr))
suppressMessages(library(fs))
suppressMessages(library(devtools))
suppressMessages(library(tidyverse))
suppressMessages(library(GGally))
suppressMessages(library(network))
suppressMessages(library(sna))

pictograph_path <- "~/GitHub/pictograph"

devtools::load_all(file.path(pictograph_path, "code", "clone.tools"))

data.dir <- file.path(pictograph_path, "data")

# ===================================================================================== #
# Command line arguments

option_list = list(
  make_option(c("-r", "--results"), type="character", default=NULL,
              help="tree MH results RDS", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL,
              help="output directory", metavar="character"),
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="run name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

if (is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Output directory must be supplied", call.=FALSE)
}
if (is.null(opt$results)){
  print_help(opt_parser)
  stop("Tree MH results must be supplied", call.=FALSE)
}
if (is.null(opt$name)){
  print_help(opt_parser)
  stop("Run name must be supplied", call.=FALSE)
}

dir_create(opt$outdir)

# ===================================================================================== #
# Data
results <- readRDS(opt$results)

# ===================================================================================== #
# Plot score chain

print("plotting score chain")

score.tb <- tibble(score = results$score.chain, 
                   iter = seq_len(length(results$score.chain)),
                   log.score = log(results$score.chain))
score.plot <- ggplot(score.tb, aes(x = iter, y = log.score)) +
  geom_line() +
  theme_light() +
  ggtitle(paste0("Acceptance rate = ", round(results$acceptRate, 2)))
ggsave(file.path(opt$outdir, paste0(opt$name, "_score.chain.pdf")), score.plot, 
       width = 7, height = 4)

# ===================================================================================== #
# Plot 10 most frequent trees

print("plotting most frequent trees")

trees <- sapply(results$admat.chain, numericRepresentation)
tab <- table(trees)
freq <- as.numeric(tab)
prob <- freq/sum(freq)

tb <- tibble(prob = prob,
             tree = names(tab))
tb <- tb %>%
  mutate(index = match(tb$tree, trees)) 
tb <- tb %>%
  mutate(score = results$score.chain[tb$index])
tb <- tb %>%
  arrange(desc(prob))

top.freq.plots <- list()
for (i in 1:10) {
  adm <- results$admat.chain[[tb$index[i]]]
  dag <- plotDAG(adm) + 
    ggtitle(paste0("Rank ", i, 
                   ", freq = ", round(tb$prob[i], 4),
                   ", log_score = ", round(log(tb$score[i]), 3)))
  top.freq.plots[[i]] <- dag
}

freq.plots <- do.call(grid.arrange, c(top.freq.plots, nrow=4))
ggsave(file.path(opt$outdir, paste0(opt$name, "_top.10.trees.pdf")), 
       freq.plots, width = 15, height = 20)

# ===================================================================================== #
# Plot posterior adjacency matrix

print("plotting posterior adjacency matrix")

adm.chain <- rapply(results$admat.chain, f=function(x) ifelse(is.na(x),0,x), how="replace")
post.admat <- Reduce("+", adm.chain)/length(adm.chain)
print("post.admat:")
print(round(post.admat,2))

post.admat.tb <- post.admat %>%
  as_tibble() %>%
  mutate(from = row.names(post.admat))
post.admat.tb <- post.admat.tb %>% 
  pivot_longer(-from, names_to = "to", values_to = "post")
# order clusters for axes
post.admat.tb <- post.admat.tb %>%
  mutate(from = factor(from, levels = c(paste0("cluster", (nrow(post.admat)-1):1), "root")),
         to = factor(to, levels = paste0("cluster", 1:(nrow(post.admat)-1))))

post.admat.plot <- ggplot(post.admat.tb, aes(x = to, y = from)) +
  geom_tile(aes(fill = post), colour = "grey50")

ggsave(file.path(opt$outdir, paste0(opt$name, "_post.admat.pdf")), 
       post.admat.plot, width = 14, height = 13)

