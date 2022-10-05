getTipLabelOrder <- function (Phylo) {
  png_null_device(4, 4)
  plot(Phylo, plot = FALSE)
  dev.off()
  LastPhyloPlot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  Y <- LastPhyloPlot$yy[1:Ntip(Phylo)]
  TipLabelsOrdered <- Phylo$tip.label[order(Y, decreasing = TRUE)]
  return(TipLabelsOrdered)
}


getAgesExtantSpecies <- function (Taxonomy, Phylo, Tol = NULL) {
  # Get edge lengths indepentendly of species assignments
  # E.g. when we have anagenesis
  EdgeMax <- aggregate(start ~ edge, data = Taxa, FUN = max)
  EdgeMin <- aggregate(end ~ edge, data = Taxa, FUN = min)
  EdgeLength <- EdgeMax[, 2] - EdgeMin[, 2]
  Taxonomy$edge_start <- EdgeLength[Taxonomy$edge]
  # What looks like end=0 is not exactly 0, so we need a tolerance 
  # to identify extant species
  if (is.null(Tol)) {
    Tol <- max(branching.times(Phylo)) / 1e6
  }
  # Keep only extant species
  ExtantSpecies <- unique(Taxonomy$sp[Taxonomy$end < Tol])
  TaxaExtantTmp <- Taxonomy[Taxonomy$sp %in% ExtantSpecies, ]
  # Order so that for each species the edge ending at the present is on top of older branches
  TaxaExtantTmp <- TaxaExtantTmp[order(TaxaExtantTmp$sp, TaxaExtantTmp$end), ]
  # A species may appear multiple times and we look for the max of start,
  # look-up the maximum of start
  TaxaExtant <- aggregate(start ~ sp, data = TaxaExtantTmp, FUN = max)
  colnames(TaxaExtant)[2] <- 'Age'
  TaxaExtant <- data.frame(TaxaExtant, TaxaExtantTmp[TaxaExtantTmp$end < Tol, -1])
  TipAgesCompl <- calculate_tip_ages(Phylo)
  colnames(TipAgesCompl) <- c('Tip_Compl', 'Tip_Length_Compl')
  # Tip ages for extant-only phylogeny
  PhyloExtant <- prune.fossil.tips(Phylo)
  TipAgesExtant <- as.data.frame(as_tibble(PhyloExtant))
  TipAgesExtant <- TipAgesExtant[!is.na(TipAgesExtant$label), ]
  colnames(TipAgesExtant)[3:4] <- c('Tip_Length_Reconstr', 'Tip_Reconstr')
  # Get number of sister tips for each tip
  TipAgesExtant$Sisters <- sapply(TipAgesExtant$parent, function(x)
    length(Descendants(PhyloExtant, x, type = 'tips')[[1]]) - 1)
  #
  TipAgesExtant <- TipAgesExtant[match(TipAgesCompl$Tip_Compl,
                                       TipAgesExtant$Tip_Reconstr), ]
  TipAges <- data.frame(TipAgesCompl, 
                        Tip_Length_Reconstr = TipAgesExtant$Tip_Length_Reconstr,
                        N_Sisters_Reconstr = TipAgesExtant$Sisters)
  # Order according to plot
  # (older species form sim.taxonomy always on the right branch)
  TipLabelsOrdered <- getTipLabelOrder(Phylo)
  TipAges <- TipAges[match(TipLabelsOrdered, TipAges$Tip_Compl), ]
  TipAges <- TipAges[!is.na(TipAges$Tip_Length_Reconstr), ]
  # Keep branch lengths of extant tips only
  # Match tip lengths and labels with true species ages
  RoundDigits <- 4
  M <- pmatch(round(TaxaExtant$edge_start, RoundDigits),
              round(TipAges$Tip_Length_Compl, RoundDigits))
  TipAges <- TipAges[M, ]
  TaxaExtant <- data.frame(TaxaExtant, TipAges)
  return(TaxaExtant)
}
