
# cladogenetic test dataset -----------------------------------------------

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
set.seed(19)
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

###simulating species taxonomy; only cladogenetic speciation

Taxa.1.list <- lapply(Tree.1, sim.taxonomy, beta = 1, lambda.a = 0)

##saving data
save(Tree.1, Taxa.1.list, file = "clado.test.Tree.taxa.RData")

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

##mode: budding = 0; cladogenetic = 1
clado.df$mode <- rep(1, nrow(clado.df))

##anagenetic = 1; not anagenetic = 0
clado.df$anag <- rep(0, nrow(clado.df))

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

##joining to the dataframe (clado.df was not anymore, so directly to predictors)
clado.predictors$turnover.est <- turnover.est

clado.predictors$div.est <- div.est


##test dataset
clado.test.df <- clado.df%>% rename(True.age = Age,
                                    Estimated.age = tip_length_reconstr,
                                    N_sisters = n_sisters_reconstr)  
##selecting predictors
clado.predictors <- clado.test.df %>% 
                      select(Estimated.age, root.age, div, turnover,
                             number.sp, mode, N_sisters, anag, turnover.est,
                             div.est) 


##clado response
clado.response <- clado.test.df %>% select(True.age)



# clado test dataset ----------------------------------

write_xlsx(cbind(clado.response, clado.predictors), "clado.test.dataset.xlsx")



# ages predictors -----------------------------------------------------------
ages.predictors <- rbind(budding.predictors, clado.predictors,
                         ana.bud.test.predictors, ana.clado.test.predictors)



write_xlsx(ages.predictors, "ages.test.predictors.xlsx")
write_csv(ages.predictors, file = "ages.test.predictors.txt")

# true.age ----------------------------------------------------------------

ages.true.df <- rbind(budding.response, clado.response, ana.bud.test.response,
                      ana.clado.test.response)

write_xlsx(ages.true.df, "ages.training.response.xlsx")
write_csv(ages.true.df, file = "ages.test.response.txt")

