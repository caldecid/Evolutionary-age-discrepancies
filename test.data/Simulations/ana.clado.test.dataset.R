
# ana.clado.test.dataset.R ------------------------------------------------

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
set.seed(33)
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
  print(i)
}

###simulating species taxonomy; with cladogenetic and anagenetic speciation
##the anagenetic rate is proportional (0.5) to the lambda value

Taxa.1.list <- Map(function(x, a) sim.taxonomy(x, beta = 1, lambda.a = 0.5*a),
                   Tree.1, lam.vector)

##saving data
save(Tree.1, Taxa.1.list, file = "ana.clado.test.Tree.taxa.RData")

##Calculating estimated and true ages
TaxaExtant <- Map(getAgesExtantSpecies, Taxa.1.list, Tree.1)

ana.clado.test.df <- do.call("rbind", TaxaExtant)

##diversification
ana.clado.test.df$div <- rep(r, Ntip)

##turnover rates
ana.clado.test.df$turnover <- rep(a, Ntip)

##number of species
ana.clado.test.df$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.1, tree.max)
ana.clado.test.df$root.age <- rep(root.ages, Ntip)

##budding speciation = 0 ; cladogenetic = 1
ana.clado.test.df$mode <- rep(1, nrow(ana.clado.test.df))
##anagenetic = 1; not anagenetic = 0
ana.clado.test.df$anag <- rep(1, nrow(ana.clado.test.df))

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
ana.clado.test.df$turnover.est <- turnover.est

ana.clado.test.df$div.est <- div.est

##renaming columns
ana.clado.test.df <- ana.clado.test.df %>% rename(True.age = Age,
                                            Estimated.age = tip_length_reconstr,
                                            N_sisters = n_sisters_reconstr) 

##selecting predictors
ana.clado.test.predictors <- ana.clado.test.df %>% select(Estimated.age, 
                                                  root.age, div, turnover, 
                                                  number.sp, mode,  
                                                  N_sisters, anag,
                                                  turnover.est, 
                                                  div.est) 

#anagenetic-clado.response
ana.clado.test.response <- ana.clado.test.df %>% select(True.age)


# anagenetic-clado training dataset ----------------------------------

write_xlsx(ana.clado.test.df, "ana.clado.test.dataset.xlsx")
