library(TreeSim)
library(FossilSim)
library(geiger)
library(picante)
library(cowplot)


source('UtilisR/calculate_tip_ages.R')
source('UtilisR/getAgesExtantSpecies.R')


# Simulate phylogeny including extant and extinct species
#########################################################
set.seed(42)

Ntips <- 25 # Number of extant sampled tips
Lambda <- 0.2 # Speciation rate
Mu <- 0.1 # Extinction rate
Tree <- sim.bd.taxa(n = Ntips, lambda = Lambda, mu = Mu,
                    complete = TRUE, numbsim = 1)[[1]]


# Map species on given phylogeny
################################
# If beta = 1, branch lengths of tips are equal to the age of the extant species
Taxa <- sim.taxonomy(Tree, beta = 0.0, lambda.a = 0) 


# Get ages of extant species and the corresponding tip lengths
##############################################################
TaxaExtant <- getAgesExtantSpecies(Taxonomy = Taxa, Phylo = Tree)