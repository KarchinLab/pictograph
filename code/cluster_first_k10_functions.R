rand.admat <- function(admat) {
  for(col in 1:ncol(admat)) {
    ind.0 <- which(admat[,col] == 0) # possible positions (0's)
    rand.ind <- sample(ind.0, size=1)
    admat[rand.ind,col] <- 1
  }
  
  while (sum(admat[1, ]) == 0) {
    admat <- mutate.admat(admat)
  }
  
  admat
}
base.admat <- function(w, zero.thresh=0.01) {
  cluster.sample.presence <- apply(w, 1, function(x) which(x>zero.thresh))
  K <- nrow(w)
  S <- ncol(w)
  all.samples <- 1:S
  admat <- matrix(data=0, nrow=(1+K), ncol=K) # rows=from is root + 1:K, cols=to is 1:K
  
  # fill in restraints
  # can go from root to anyone, skip and start at nrow=2 (cluster 1)
  for (from in 2:(K+1)) {
    for (to in 1:K) {
  
      # can't go to self 
      if ((from-1) == to) {
        admat[from, to] <- NA
        #print(c(from, to, "self"))
        next
      }
      # hierarchy restraints
      from.samples <- cluster.sample.presence[[from-1]]
      to.samples <- cluster.sample.presence[[to]]
      
      ## no restraints if same sample presence
      if (setequal(from.samples, to.samples)) {
        #print(c(from, to, "same")) 
        next
      }
      
      ## restraint if # from.samples < # to.samples
      if(length(from.samples) < length(to.samples)) {
        #print(c(from, to, "from set is smaller than to set")) 
        admat[from, to] <- NA 
        next
      }
      
      ## no restraints if to.samples is subset of from.samples
      if (all(to.samples %in% from.samples)) {
        #print(c(from, to, "subset")) 
        next
      } else {
        #print(c(from, to, "not subset"))
        admat[from, to] <- NA 
      }
    }
  }
  admat
}

init.admat <- function(w, zero.thresh) {
  base <- base.admat(w, zero.thresh)
  rand.admat(base)
}


mutate.admat <- function(admat) {
  # choose a column to mutate
  K <- ncol(admat)
  rand.k <- sample(1:K, size=1)
  
  # mutate
  ## possible positions (0's)
  possiblePos <- which(!is.na(admat[, rand.k]) & admat[, rand.k] != 1)
  ## current position with 1
  ind.1 <- which(admat[, rand.k] == 1)
  ## select new position
  if (length(possiblePos) == 1) {
    new.1 <- possiblePos
  } else {
    new.1 <- sample(possiblePos, size=1)
  }
  
  new.admat <- admat 
  new.admat[ind.1, rand.k] <- 0
  new.admat[new.1, rand.k] <- 1
  
  while (sum(new.admat[1, ]) == 0) {
    new.admat <- mutate.admat(admat)
  }
  new.admat
}

getChildrenWSum <- function(admat, w, row) {
  curr.row <- admat[row,]
  children <- which(curr.row == 1)
  children.w <- w[children,]
  if(length(children) > 1) {
    return(colSums(children.w))
  } else {
    return(children.w)
  }
}

score.admat <- function(admat, w) {
  
  S <- ncol(w)
  # root, MCF=1.0
  children.w.sum <- getChildrenWSum(admat, w, 1)
  score <- log(sum(1 >= children.w.sum) / S)
  
  for(i in 2:nrow(admat)) {
    if(all(is.na(admat[i,]))) next #leaf
    
    curr.child.sum <- getChildrenWSum(admat, w, i)
    curr.parent.mcf <- w[i-1, ]
    curr.score <- log(sum(curr.parent.mcf >= curr.child.sum) / S)
    score <- score + curr.score
  }
  score
}

count.from.nodes <- function(admat) {
    sum(rowSums(admat, na.rm = T) > 0)
}

subset.w.chain <- function(w.chain.tb, k, s) {
  # returns list -- each item of list is a sample
  if(length(s) == 1) return(list(w.chain.tb[paste0("w[", k,",", s, "]")]))
  
  lapply(s, function(sample) w.chain.tb[paste0("w[", k,",", sample, "]")])
}

get.w.chain.sum <- function(w.chain.list) {
  if(dim(w.chain.list[[1]])[2] > 1) return(lapply(w.chain.list, rowSums))
  w.chain.list
}

get.children <- function(admat, parent.node) {
  # parent.node: root=1, clusters start at 2
  # returns children clusters
  which(admat[parent.node, ] == 1)
}

get.children.w.sum.chain <- function(admat, w.chain.tb, parent.node, S) {
  # parent.node: root=1, clusters start at 2
  children <- get.children(admat, parent.node)
  children.w <- subset.w.chain(w.chain.tb, k=children, s=1:S)
  get.w.chain.sum(children.w)
}

rule.chain.freq <- function(parent.chain, child.chain) {
  rule <- mapply(function(parent, child) parent > child, parent=parent.chain, child=child.chain)
  sum(rule) / prod(dim(rule))
}

score.admat.chain <- function(admat, w.chain) {
    w.ind <- colnames(w.chain)
    k.vec <- gsub("^w\\[(.*)\\,.*", "\\1", w.ind)
    K <- length(unique(k.vec))
    S <- length(k.vec) / K
    total.obs <- nrow(w.chain) * S

    #numFromNodes <- count.from.nodes(admat)
    w.chain.tb <- as_tibble(w.chain)
    # root, MCF=1.0
    children.w.sum <- get.children.w.sum.chain(admat, w.chain.tb, 1, S)
    #parent.w <- subset.w.chain(w.chain.tb, 1, 1:S)
    rule.freq <- rule.chain.freq(parent.chain=1, child.chain=children.w.sum)

    #score <- log(rule.freq)
    score <- rule.freq

    for(i in 2:nrow(admat)) {
      if(sum(admat[i,], na.rm=T) == 0) next #leaf

      curr.child.sum <- get.children.w.sum.chain(admat, w.chain.tb, i, S)
      parent.w <- subset.w.chain(w.chain.tb, i-1, 1:S)
      curr.freq <- rule.chain.freq(parent.chain=parent.w, child.chain=curr.child.sum)
      #score <- score + log(curr.freq)
      score <- score * curr.freq
    }
    score
}

scoreA <- function(A, w){
    w.ind <- colnames(w.chain)
    k.vec <- gsub("^w\\[(.*)\\,.*", "\\1", w.ind)
    K <- length(unique(k.vec))
    S <- length(k.vec) / K
    total.obs <- nrow(w.chain) * S

    #numFromNodes <- count.from.nodes(A)
    w.chain.tb <- as_tibble(w.chain)
    # root, MCF=1.0
    children.w.sum <- get.children.w.sum.chain(A, w.chain.tb, 1, S)
    #parent.w <- subset.w.chain(w.chain.tb, 1, 1:S)
    rule.freq <- rule.chain.freq(parent.chain=1, child.chain=children.w.sum)

    #score <- log(rule.freq)
    score <- rule.freq

    for(i in 2:nrow(A)) {
      if(sum(A[i,], na.rm=T) == 0) next #leaf

      curr.child.sum <- get.children.w.sum.chain(A, w.chain.tb, i, S)
      parent.w <- subset.w.chain(w.chain.tb, i-1, 1:S)
      curr.freq <- rule.chain.freq(parent.chain=parent.w, child.chain=curr.child.sum)
      #score <- score + log(curr.freq)
      score <- score * curr.freq
    }
    score
}

#score.admat.chain(admat, w.chain)

check.base.admat <- function(admat, base) {
 all(which(is.na(admat)) == which(is.na(base)))
}

clusterIndex <- function(x){
    as.character(x) %>%
        strsplit("[", fixed=TRUE) %>%
        sapply("[", 2) %>%
        strsplit(",") %>%
        sapply("[", 1) %>%
        as.integer()
}
sampleIndex <- function(x){
    as.character(x) %>%
        strsplit(",", fixed=TRUE) %>%
        sapply("[", 2) %>%
        str_replace("]", "") %>%
        as.integer()
}

mutationIndex <- function(x){
    as.character(x) %>%
        strsplit("[", fixed=TRUE) %>%
        sapply("[", 2) %>%
        str_replace("]", "") %>%
        as.integer()
}

plot.z <- function(samps, z) {
  mcmc_vals <- summary(samps)$statistics
  mcmc_z <- as.vector(mcmc_vals[substr(rownames(mcmc_vals), 1, 1) == "z", "Mean"])
  plot(z, mcmc_z, type = "p")
  z_comp <- data.frame(z, mcmc_z)
}

z.chain.to.tb <- function(z.chain) {
  z.chain.tb <- z.chain %>%
    as_tibble() %>%
    mutate(iter=1:nrow(z.chain)) %>%
      gather(variant, mcmc_z, -c(iter))
  L <- length(unique(z.chain.tb$variant))
  variantIndex <- function(x){
      xx <- strsplit(x, "[", fixed=TRUE)
      x2 <- sapply(xx, "[", 2) %>%
          gsub("]", "", .) %>%
          as.integer()
  }
  K <- 3
  z2 <- z.chain.tb %>%
      mutate(variant_index = variantIndex(variant)) %>%
      mutate(true_z = rep(1:10, each=nrow(z.chain)*K)) %>%
      group_by(variant, mcmc_z) %>% 
      mutate(count = n())
  z.chain.tb_simp <- distinct(select(z2, -c(iter)))
  z.chain.tb_simp %>% 
    group_by(variant) %>%
    mutate(prop = count/sum(count))
}
