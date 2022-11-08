
# budding test dataset ----------------------------------------------------

library(TreeSim)
library(FossilSim)
library(ape)
library(picante)
library(phytools)
library(tidyverse)
library(purrr)
library(ggplot2)
library(geiger)
library(readr)
library(writexl)
library(tidytree)
library(phangorn)

##setting seed for replication
set.seed(17)
##diversification rates random vector
r <- runif(25, 0, 1)

##turnover rates random vector
a <- runif(25, 0, 1)

##lambda vector
lam.vector <- r/(1-a)

##mu vector
mu.vector <- -a*r/(a-1)

##Ntip random vector
Ntip <- round(runif(25, min = 20, max = 100))

##simulating 100 trees
Tree.1 <- vector("list", 25)

##for loop due to different extinction rates
for(i in seq_along(mu.vector)){
  Tree.1[i] <- sim.bd.taxa(n = Ntip[i], lambda = lam.vector[i], mu = mu.vector[i],
                           complete = TRUE, numbsim = 1)
}

###simulating species taxonomy; only budding speciation

Taxa.1.list <- lapply(Tree.1, sim.taxonomy, beta = 0, lambda.a = 0)

##saving Trees and taxonomies
save(Tree.1, Taxa.1.list, file = "bud.test.tree.taxa.RData")

##Calculating estimated and true ages
TaxaExtant <- Map(getAgesExtantSpecies, Taxa.1.list, Tree.1)

budding.df <- do.call("rbind", TaxaExtant)

##diversification
budding.df$div <- rep(r, Ntip)

##turnover rates
budding.df$turnover <- rep(a, Ntip)

##number of species
budding.df$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.1, tree.max)
budding.df$root.age <- rep(root.ages, Ntip)

##anagenetic speciation =1 , not anagenetic = 0
budding.df$anag <- rep(0, nrow(budding.df))

##mode: budding = 0; cladogenetic = 1
budding.df$mode <- rep(0, nrow(budding.df))

##estimating turnover and diversification rates from the sim trees
##prunning fossils
Tree.extant <- lapply(Tree.1, prune.fossil.tips)

##estimating turnover and div
bd <- lapply(Tree.extant, birthdeath)

bd.turnover <- rep(NA, length(bd))
bd.div <- rep(NA, length(bd))

for(i in seq_along(bd)){
  bd.turnover[i] <- bd[[i]]$para[1]
  bd.div[i] <- bd[[i]]$para[2]
}

##repeating the rate to match the species
Ntip <- sapply(Tree.extant, Ntip)
turnover.est <- rep(bd.turnover, Ntip)
div.est <- rep(bd.div, Ntip)

##joining to the dataframe
budding.testing.df$turnover.est <- turnover.est
budding.testing.df$div.est <- div.est

##arranging dataset 
budding.testing.df <- budding.df %>% rename(True.age = Age,
                                            Estimated.age = tip_length_reconstr,
                                            N_sisters = n_sisters_reconstr)  

##selecting predictors
budding.predictors <- budding.testing.df %>% 
                               select(Estimated.age, root.age, div, turnover,
                                             number.sp, mode, N_sisters, anag,
                                      turnover.est, div.est) 


##budding response
budding.response <- budding.testing.df %>% select(True.age)


# budding testing dataset ----------------------------------

write_xlsx(budding.testing.df, "budding.testing.dataset.xlsx")
