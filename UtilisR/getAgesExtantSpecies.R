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
  # What looks like end=0 is not exactly 0, so we need a tolerance
  if (is.null(Tol)) {
    Tol <- max(branching.times(Phylo)) / 1e6
  }
  ExtantSpecies <- unique(Taxonomy$sp[Taxonomy$end < Tol])
  TaxaExtantTmp <- Taxonomy[Taxonomy$sp %in% ExtantSpecies, ]
  TaxaExtantTmp <- TaxaExtantTmp[order(TaxaExtantTmp$sp, TaxaExtantTmp$end), ]
  # A species may appear multiple times and we look for the max of start,
  # look-up the maximum of start
  TaxaExtant <- aggregate(start ~ sp, data = TaxaExtantTmp, FUN = max)
  colnames(TaxaExtant)[2] <- 'Age'
  TaxaExtant <- data.frame(TaxaExtant, TaxaExtantTmp[TaxaExtantTmp$end < Tol, -1])
  TipAgesCompl <- calculate_tip_ages(Phylo)
  colnames(TipAgesCompl) <- c('Tip_Compl', 'Tip_Length_Compl')
  PhyloExtant <- prune.fossil.tips(Phylo)
  TipAgesExtant <- calculate_tip_ages(PhyloExtant)
  colnames(TipAgesExtant) <- c('Tip_Reconstr', 'Tip_Length_Reconstr')
  TipAgesExtant <- TipAgesExtant[match(TipAgesCompl$Tip_Compl, TipAgesExtant$Tip_Reconstr), ]
  TipAges <- data.frame(TipAgesCompl, Tip_Length_Reconstr = TipAgesExtant$Tip_Length_Reconstr)
  TipLabelsOrdered <- getTipLabelOrder(Phylo)
  TipAges <- TipAges[match(TipLabelsOrdered, TipAges$Tip_Compl), ]
  TipAges <- TipAges[!is.na(TipAges$Tip_Length_Reconstr), ]
  # Keep branch lengths of extant tips only
  # Match tip lengths and labels with species ages
  RoundDigits <- 4
  M <- pmatch(round(TaxaExtant$start, RoundDigits),
              round(TipAges$Tip_Length_Compl, RoundDigits))
  TipAges <- TipAges[M, ]
  TaxaExtant <- data.frame(TaxaExtant, TipAges)
  return(TaxaExtant)
}
