getTipLabelOrder <- function (Phylo) {
  png_null_device(4, 4)
  plot(Phylo, plot = FALSE)
  dev.off()
  LastPhyloPlot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  Y <- LastPhyloPlot$yy[1:Ntip(Phylo)]
  TipLabelsOrdered <- Phylo$tip.label[order(Y, decreasing = TRUE)]
  return(TipLabelsOrdered)
}


getAgesExtantSpecies <- function(Taxonomy, Phylo, Tol = NULL) {
  # Assign names of tips and nodes from the phylogeny to the species
  Nt <- Ntip(Phylo)
  Names <- c(Phylo$tip, as.character((Nt + 1):(Nt + Phylo$Nnode)))
  Taxonomy$label <- Names[Taxa$edge]
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
  # Tip ages for phylogeny
  TipAgesCompl <- calculate_tip_ages(Phylo)
  colnames(TipAgesCompl) <- c('tip_compl', 'tip_length_compl')
  # Add those tip ages that are for extant tips
  M <- match(TaxaExtant$label, TipAgesCompl$tip_compl)
  TaxaExtant <- data.frame(TaxaExtant, tip_length_compl = TipAgesCompl[M, 2])
  # Get number of sister tips for each tip
  PhyloExtant <- prune.fossil.tips(Phylo)
  PhyloExtantDF <- as.data.frame(as_tibble(PhyloExtant))
  PhyloExtantDF <- PhyloExtantDF[!is.na(PhyloExtantDF$label), ]
  M <- match(TaxaExtant$label, PhyloExtantDF$label)
  PhyloExtantDF <- PhyloExtantDF[M, ]
  TaxaExtant$n_sisters_reconstr <- sapply(PhyloExtantDF$parent, function(x)
    length(Descendants(PhyloExtant, x, type = 'tips')[[1]]) - 1)
  # Add tip length of reconstructed tree
  TaxaExtant$tip_length_reconstr <- PhyloExtantDF$branch.length
  return(TaxaExtant)
}
