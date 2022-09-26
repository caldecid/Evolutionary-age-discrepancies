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

##selecting predictors
budding.predictors <- budding.df %>% rename(True.age = Age,
                                      Estimated.age = Tip_Length_Reconstr)  %>% 
                     select(Estimated.age, root.age, div, turnover, number.sp, mode) 

##budding speciation = 0 ; cladogenetic = 1
budding.predictors$mode <- rep(0, nrow(budding.predictors))

##budding response
budding.true.age <- as_tibble(budding.df$Age)
colnames(budding.true.age) <- "true.age"

budding.true.age$mode <- rep(0, nrow(budding.predictors))


# budding training dataset ----------------------------------

write_xlsx(budding.df, "budding.training.dataset.xlsx")
