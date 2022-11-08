# Anagenetic-cladogenetic speciation with diversification and turnover rates -------------------
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
library(phangorn)
library(tidytree)

##setting seed for replication
set.seed(20)
##diversification rates random vector
r <- runif(100, 0, 1)

##turnover rates random vector
a <- runif(100, 0, 1)

##lambda vector
lam.vector <- r/(1-a)

##mu vector
mu.vector <- -a*r/(a-1)

##Ntip random vector
Ntip <- round(runif(100, min = 20, max = 100))

##simulating 100 trees
Tree.1 <- vector("list", 100)

##for loop due to different extinction rates
for(i in seq_along(mu.vector)){
  Tree.1[i] <- sim.bd.taxa(n = Ntip[i], lambda = lam.vector[i], mu = mu.vector[i],
                           complete = TRUE, numbsim = 1)
  print(i)
}

###simulating species taxonomy; with cladogenetic and anagenetic speciation
##the anagenetic rate is proportional (0.5) to the lambda value

Taxa.1.list <- Map(function(x, a) sim.taxonomy(x, beta = 1, lambda.a = 0.5*a),
                   Tree.1, lam.vector)


save(Tree.1, Taxa.1.list, file = "ana.cladogenetic.train.Tree.Taxa.RData")

##Calculating estimated and true ages
TaxaExtant <- Map(getAgesExtantSpecies, Taxa.1.list, Tree.1)

ana.clado.train.df <- do.call("rbind", TaxaExtant)

##diversification
ana.clado.train.df$div <- rep(r, Ntip)

##turnover rates
ana.clado.train.df$turnover <- rep(a, Ntip)

##number of species
ana.clado.train.df$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.1, tree.max)
ana.clado.train.df$root.age <- rep(root.ages, Ntip)
##budding speciation = 0 ; cladogenetic = 1
ana.clado.train.df$mode <- rep(1, nrow(ana.clado.train.df))
##anagenetic = 1; not anagenetic = 0
ana.clado.train.df$anag <- rep(1, nrow(ana.clado.train.df))
##renaming columns
ana.clado.train.df <- ana.clado.train.df %>% rename(True.age = Age,
                                Estimated.age = tip_length_reconstr,
                                N_sisters = n_sisters_reconstr) 

##selecting predictors
ana.clado.train.predictors <- ana.clado.train.df %>% select(Estimated.age, root.age,
                                        div, turnover, number.sp, mode,
                                        N_sisters, anag) 

#anagenetic-clado.response
ana.clado.train.response <- ana.clado.train.df %>% select(True.age)


# anagenetic-clado training dataset ----------------------------------

write_xlsx(ana.clado.train.df, "ana.clado.training.dataset.xlsx")



