
# Budding speciation diversification and turnover rates -------------------
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
set.seed(16)
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
}

###simulating species taxonomy; only budding speciation

Taxa.1.list <- lapply(Tree.1, sim.taxonomy, beta = 0, lambda.a = 0)


##saving taxonomies and Trees
save(Tree.1, Taxa.1.list, file = "budding.train.Tree.Taxa.RData")

##get extant species


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
##budding speciation = 0 ; cladogenetic = 1
budding.df$mode <- rep(0, nrow(budding.df))
##anagenetic = 1; not anagenetic = 0
budding.df$anag <- rep(0, nrow(budding.df))

##Renaming dataframe
budding.df <- budding.df %>% rename(True.age = Age,
                                      Estimated.age = Tip_Length_Reconstr,
                                      N_sisters = N_Sisters_Reconstr)  

##predictor
budding.predictors <- budding.df %>% 
                     select(Estimated.age, root.age, div, turnover, number.sp,
                            mode, N_sisters, anag) 



##budding response
budding.reponse <- budding.df %>% select(True.age)



# budding training dataset ----------------------------------

write_xlsx(budding.df, "budding.training.dataset.xlsx")

