# Evolutionary age discrepancies ------------------------------------------

library(TreeSim)
library(FossilSim)
library(ape)
library(picante)
library(phytools)
library(tidyverse)


##### Simulate phylogeny 

Ntips <- 200 # Number of extant sampled tips
Lambda <- 0.2 # Speciation rate
Mu <- 0.02 # Extinction rate
tree.1 <- sim.bd.taxa(n = Ntips, lambda = Lambda, mu = Mu,
                     complete = TRUE, numbsim = 10)[[1]]


###simulating species taxonomy (forms of speciation)
Taxa <- sim.taxonomy(tree.1, beta = 0, lambda.a = 0) 
#Bifurcation = beta
#Anagenetic = lambda.a
#Budding = 1-beta

Taxa$sp <- sub("^", "t", Taxa$sp)##for matching posterior datasets

##Prunning the extinct species
tree.extant <- prune.fossil.tips(tree.1)

##species might be associated with multiple edges;
##for calculating the evolutionary age we should pick the oldest edge

extant.species <- Taxa[Taxa$sp %in% tree.extant$tip.label, ] %>% 
                       group_by(sp) %>%
                       filter(start == max(start)) 


#function by Tanentzap 2019 to calculate species age
calculate_tip_ages <- function(tree){
  phy.age <- node.age(tree)
  BL.position <- cbind(phy.age$edge,phy.age$age, tree$edge.length)
  dist.tip <- max(phy.age$age)-BL.position[,3]
  BL.positions <- cbind(BL.position,dist.tip)
  ages<- BL.positions[,5] + BL.positions[,4]
  BL.positions <- cbind(BL.positions,ages)
  node.ages<- as.data.frame(BL.positions)
  names(node.ages) <- c("parental.node","daughter.node",
                        "dist.root","BL",
                        "dist.tip","mrca.age")
  ## node.ages is a data frame listing as variables the identity of parental and
  #daughter nodes, the distance from the root and from the present of each node,
  #the branch length and the age of the most recent common ancestor
  species.ages<- node.ages[node.ages[,2] < length(tree$tip)+1,]
  row.names(species.ages) <- tree$tip[species.ages$daughter.node]
  ## species ages is node.ages data frame reduced to the tips (species)
  species.ages <- species.ages[order(row.names(species.ages)),]
  output.table <- as.data.frame(cbind(row.names(species.ages),
                                      species.ages$mrca.age))
  colnames(output.table) <- c('sp','phylo.age')
  output.table$phylo.age <- as.numeric(output.table$phylo.age)
  return(output.table)
}

###phylo age calculation from the Tanentzap's function
neo.tip.age <- calculate_tip_ages(tree.extant)


###joining datasets and creating a measure of discrepancy
age.discrepancies <- left_join(extant.species, neo.tip.age,
                               by = "sp") %>%
                      dplyr::select(sp, start, phylo.age) %>%
                      rename(real.age = start) 
