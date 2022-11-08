# Cladogenetic speciation diversification and turnover rates -------------------
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
set.seed(18)
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

###simulating species taxonomy; only cladogenetic speciation

Taxa.1.list <- lapply(Tree.1, sim.taxonomy, beta = 1, lambda.a = 0)

save(Tree.1, Taxa.1.list, file = "cladogenetic.train.Tree.Taxa.RData")

##Calculating estimated and true ages
TaxaExtant <- Map(getAgesExtantSpecies, Taxa.1.list, Tree.1)

clado.df <- do.call("rbind", TaxaExtant)

##diversification
clado.df$div <- rep(r, Ntip)

##turnover rates
clado.df$turnover <- rep(a, Ntip)

##number of species
clado.df$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.1, tree.max)
clado.df$root.age <- rep(root.ages, Ntip)
##budding speciation = 0 ; cladogenetic = 1
clado.df$mode <- rep(1, nrow(clado.df))
##anagenetic = 1; not anagenetic = 0
clado.df$anag <- rep(0, nrow(clado.df))
##renaming columns
clado.df <- clado.df %>% rename(True.age = Age,
                               Estimated.age = Tip_Length_Reconstr,
                               N_sisters = N_Sisters_Reconstr) 

##selecting predictors
clado.predictors <- clado.df %>% select(Estimated.age, root.age,
                                        div, turnover, number.sp, mode,
                                        N_sisters, anag) 

#clado.response
clado.response <- clado.df %>% select(True.age)


# clado training dataset ----------------------------------

write_xlsx(clado.df, "clado.training.dataset.xlsx")



# ages predictors -----------------------------------------------------------
ages.predictors <- rbind(budding.predictors, clado.predictors)

write_xlsx(ages.predictors, "ages.training.predictors.xlsx")

# true.age ----------------------------------------------------------------

ages.true.df <- rbind(budding.true.age, clado.true.age)

write_xlsx(ages.true.df, "ages.training.response.xlsx")
