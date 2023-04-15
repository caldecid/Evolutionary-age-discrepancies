weightedQuantile <- function(x, w, prob = 0.5) {
  # Forked from modi
  # sort the data and the weights in ascending order
  ord <- order(x)
  w <- w[ord]
  x <- x[ord]
  
  # compute weighted cdf
  w.ord <- cumsum(w) / sum(w)
  
  # get index with length of data vector x
  index <- seq_along(x)
  
  # if the min of the cdf is larger than prob,
  # then the quantile is based on only one observation (Warning)
  # thus, set index of lower quant to 1 (first observation)
  # otherwise, take max index where w.ord smaller or equal to prob
  lower.k.quant <- 1
  if (min(w.ord) <= prob) {
    lower.k.quant <- max(index[w.ord <= prob])
  } 
  
  # set index of upper quant to min index where w.ord greater than prob
  upper.k.quant <- min(index[w.ord > prob])
  
  # if inverse is non-unique, return weighted linear interpolation
  if (w.ord[lower.k.quant] < prob) {
    # returns quantile according to upper.k.quant
    return(x[upper.k.quant])
  }
  else {
    # returns weighted linear interpolation
    return((w[lower.k.quant] * x[lower.k.quant] + w[upper.k.quant] * x[upper.k.quant]) /
             (w[lower.k.quant] + w[upper.k.quant]))
  }
}


# Forked from BAMM tools to get correct branching times for non-ultrametric trees
nuBranchingTimes <- function(phy, return.type = 'bt'){
  if (!is.binary.phylo(phy)){
    stop("error. Need fully bifurcating (resolved) tree\n")
  }
  phy$begin <- rep(0, nrow(phy$edge)) 
  phy$end <-  rep(0, nrow(phy$edge))
  # Do it recursively
  fx <- function(phy, node){
    cur.time <- 0
    root <- length(phy$tip.label) + 1
    if (node > root){
      cur.time <- phy$end[which(phy$edge[,2] == node)]	
    }
    dset <- phy$edge[,2][phy$edge[,1] == node]
    i1 <- which(phy$edge[,2] == dset[1])
    i2 <- which(phy$edge[,2] == dset[2])
    phy$end[i1] <- cur.time + phy$edge.length[i1]
    phy$end[i2] <- cur.time + phy$edge.length[i2]
    if (dset[1] > length(phy$tip.label)){
      phy$begin[phy$edge[,1] == dset[1]] <- phy$end[i1]
      phy <- fx(phy, node = dset[1])	
    }
    if (dset[2] > length(phy$tip.label)){
      phy$begin[phy$edge[,1] == dset[2]] <- phy$end[i2]
      phy <- fx(phy, node = dset[2])
    }
    return(phy)
  }
  phy <- fx(phy, node = length(phy$tip.label) + 1)
  if (return.type == 'bt'){
    maxbt <- max(phy$end)
    nodes <- (length(phy$tip.label) + 1):(2*length(phy$tip.label) - 1)
    bt <- numeric(length(nodes))
    names(bt) <- nodes
    for (i in 1:length(bt)){
      tt <- phy$begin[phy$edge[,1] == nodes[i]][1]
      bt[i] <- maxbt - tt
    }
    return(bt)
  }else if (return.type == 'begin.end'){
    return(phy)		
  }
}


getAncestralNodes <- function(tree) {
  anc <- Ancestors(tree, 1:Ntip(tree))
  names(anc) <- tree$tip
  anc <- anc[order(names(anc))]
  return(anc)
}


getPotentialAges <- function(tree, anc, bt, tl, rootedge = FALSE) {
  nt <- Ntip(tree)
  a <- lapply(anc, function(x) bt[x - nt])
  names_a <- names(a)
  # add root
  if (rootedge && exists("root.edge", tree)) {
    rootlength <- tree$root.edge
    a <- lapply(a, function(x) c(x, x[length(x)] + rootlength))
  }
  # remove time between extinction and present
  diff_to_present <- sapply(1:length(a), function(x) a[[x]][1] - tl[x])
  a <- sapply(1:length(a), function(x) a[[x]] - diff_to_present[x])
  names(a) <- names_a
  return(a)
}


getTerminalEdgeLengths <- function(tree) {
  # There must be an easier way!
  edge <- data.frame(tree$edge, length = tree$edge.length, name = NA_character_)
  edge_idx <- edge[, 2]
  is_tip <- edge_idx <= Ntip(tree)
  edge_idx <- edge_idx[is_tip]
  edge[is_tip , 4] <- tree$tip.label[edge_idx]
  edge <- edge[is_tip , ]
  edge <- edge[order(edge$name), ]
  term_edge <- edge[, 3]
  names(term_edge) <- edge[, 4]
  return(term_edge)
}


getProbNodes <- function(ages) {
  p <- lapply(ages, function(x) dgeom(0:(length(x) - 1), 0.5))
  return(p)
}


getModalAge <- function (ages) {
  m <- lapply(ages, function(x) x[1])
  m <- unlist(m)
  return(m)
}


getMeanAge <- function(ages, prob) {
  m <- sapply(1:length(ages), function(x)
    sum(ages[[x]] * prob[[x]]/sum(prob[[x]])))
  return(m)
}


getCredintAges <- function(ages, prob) {
  lwr <- sapply(1:length(ages), function(x)
    weightedQuantile(ages[[x]], prob[[x]], 0.025))
  upr <- sapply(1:length(ages), function(x)
    weightedQuantile(ages[[x]], prob[[x]], 0.975))
  return(cbind(lwr, upr))
}




getGeometricAges <- function(tree, root.edge = TRUE) {
  anc_nodes <- getAncestralNodes(tree)
  if (!is.binary(tree)) {
    tree <- multi2di(tree)
  }
  branch_times <- nuBranchingTimes(tree)
  terminal_edge_lengths <- getTerminalEdgeLengths(tree) 
  pot_ages <- getPotentialAges(tree,
                               anc_nodes,
                               branch_times,
                               terminal_edge_lengths,
                               root.edge)
  prob_anc_nodes <- getProbNodes(pot_ages)
  ages_summary <- as.data.frame(matrix(NA_real_, nrow = Ntip(tree), ncol = 4))
  colnames(ages_summary) <- c('modal', 'mean', 'lwr_ci', 'upr_ci')
  rownames(ages_summary) <- names(pot_ages)
  ages_summary[, 1] <- getModalAge(pot_ages)
  ages_summary[, 2] <- getMeanAge(pot_ages, prob_anc_nodes)
  ages_summary[, 3:4] <- getCredintAges(pot_ages, prob_anc_nodes)
  return(ages_summary)
}
