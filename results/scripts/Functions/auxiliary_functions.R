###############################################################################
############                                                    ###############
###########             Auxiliary functions                     ###############
###############################################################################

# calculate tip ages--------------------------------------

##function for calculating the branch lenght (myrs) for each tip of a calibrated
## phylogeny

calculate_tip_ages <- function(tree) {
  phy.age <- node.age(tree)
  BL.position <- cbind(phy.age$edge, phy.age$age, tree$edge.length)
  dist.tip <- max(phy.age$age) - BL.position[, 3]
  BL.positions <- cbind(BL.position, dist.tip)
  ages <- BL.positions[, 5] + BL.positions[, 4]
  BL.positions <- cbind(BL.positions,ages)
  node.ages <- as.data.frame(BL.positions)
  colnames(node.ages) <- c("parental.node","daughter.node",
                           "dist.root","BL",
                           "dist.tip","mrca.age")
  ## node.ages is a data frame listing as variables the identity of parental and
  #daughter nodes, the distance from the root and from the present of each node,
  #the branch length and the age of the most recent common ancestor
  species.ages<- node.ages[node.ages[,2] < length(tree$tip) + 1, ]
  rownames(species.ages) <- tree$tip[species.ages$daughter.node]
  ## species ages is node.ages data frame reduced to the tips (species)
  species.ages <- species.ages[order(row.names(species.ages)), ]
  output.table <- as.data.frame(cbind(row.names(species.ages),
                                      species.ages$mrca.age))
  colnames(output.table) <- c('tip','tip.age')
  output.table$tip.age <- as.numeric(output.table$tip.age)
  return(output.table)
}



# get tip label order -----------------------------------------------------


getTipLabelOrder <- function (Phylo) {
  png_null_device(4, 4)
  plot(Phylo, plot = FALSE)
  dev.off()
  LastPhyloPlot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  Y <- LastPhyloPlot$yy[1:Ntip(Phylo)]
  TipLabelsOrdered <- Phylo$tip.label[order(Y, decreasing = TRUE)]
  return(TipLabelsOrdered)
}


# share legends -----------------------------------------------------------


#https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }


# Estimating lambda and mu from an extinction fr and sampling fr  --------

#' @title Speciation and extinction rates from phylogeny
#'
#' @description This function estimates speciation and extinction rate from a phylogeny under a given extinction and sampling fraction
#' 
#' @param phy An ultrametric bifurcating phylogenetic tree, in ape "phylo" format.
#' @param epsilon Extinction fraction
#' @param rho Sampling fraction
#' @param ml_optim Method to use for optimisation. May be one of "optim", "subplex", "nlminb", "nlm" (partial unambigious string is allowed).
#' 
#' @return A named vector of two parameters, lambda and mu.


estimate_bd <- function(phy, epsilon = 0, rho = 1, ml_optim = 'subplex') {
  if(epsilon >= 1) {
    stop("Extinction fraction should not be greater than 1")
  }
  # Birth-death likelihood
  bd_lik <- make.bd(tree = phy, sampling.f = rho) 
  # Constrain extinction fraction
  con <- paste0("mu ~ ", epsilon, " * lambda")
  bd_lik <- constrain(bd_lik, con)
  # Initial birth rate for ML search with yule rate
  b_init <- phy$Nnode / sum(phy$edge.length)
  if (epsilon > 0) {
    b_init <- b_init / epsilon
  }
  # Finde maximum likelihood rates
  bd_fit <- find.mle(bd_lik, x.init = b_init, method = ml_optim)
  # Format birth and death rate for output
  lambda <- coef(bd_fit)
  mu <- epsilon * lambda
  bd_rates <- c(lambda, mu)
  names(bd_rates) <- c("lambda", "mu")
  return(bd_rates)
}

