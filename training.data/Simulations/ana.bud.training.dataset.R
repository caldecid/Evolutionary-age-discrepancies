# Anagenetic-budding speciation with diversification and turnover rates -----
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
set.seed(21)
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

###simulating species taxonomy; with budding and anagenetic speciation
##the anagenetic rate is proportional (0.5) to the lambda value

Taxa.1.list <- Map(function(x, a) sim.taxonomy(x, beta = 0, lambda.a = 0.5*a),
                   Tree.1, lam.vector)


save(Tree.1, Taxa.1.list, file = "ana.budding.train.Tree.Taxa.RData")

##Calculating estimated and true ages
TaxaExtant <- Map(getAgesExtantSpecies, Taxa.1.list, Tree.1)

##binding the lists
ana.bud.train.df <- do.call("rbind", TaxaExtant)

##diversification
ana.bud.train.df$div <- rep(r, Ntip)

##turnover rates
ana.bud.train.df$turnover <- rep(a, Ntip)

##number of species
ana.bud.train.df$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.1, tree.max)
ana.bud.train.df$root.age <- rep(root.ages, Ntip)
##budding speciation = 0 ; cladogenetic = 1
ana.bud.train.df$mode <- rep(0, nrow(ana.bud.train.df))
##anagenetic = 1; not anagenetic = 0
ana.bud.train.df$anag <- rep(1, nrow(ana.bud.train.df))
##renaming columns
ana.bud.train.df <- ana.bud.train.df %>% rename(True.age = Age,
                                          Estimated.age = tip_length_reconstr,
                                          N_sisters = n_sisters_reconstr) 

##selecting predictors
ana.bud.train.predictors <- ana.bud.train.df %>% select(Estimated.age,
                                                        root.age, div, turnover,
                                                        number.sp, mode,
                                                        N_sisters, anag) 

#anagenetic-budding response
ana.bud.train.response <- ana.bud.train.df %>% select(True.age)

# anagenetic-budding training dataset ----------------------------------

write_xlsx(ana.bud.train.df, "ana.bud.training.dataset.xlsx")


# ages predictors -----------------------------------------------------------
ages.predictors <- rbind(budding.predictors, clado.predictors,
                         ana.bud.train.predictors, ana.clado.train.predictors)

rownames(ages.predictors)<- NULL

write_xlsx(ages.predictors, "ages.training.predictors.xlsx")
write_csv(ages.predictors, file = "ages.training.predictors.txt")

# true.age ----------------------------------------------------------------

ages.true.df <- rbind(budding.reponse, clado.response, ana.bud.train.response,
                      ana.clado.train.response)
rownames(ages.true.df) <- NULL
write_xlsx(ages.true.df, "ages.training.response.xlsx")
write_csv(ages.true.df, file = "ages.training.response.txt")
