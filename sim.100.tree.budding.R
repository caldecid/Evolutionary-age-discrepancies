# First simulation, 100 trees, only budding speciation--------------------------
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

###############Tree simulation##################
set.seed(13)

Ntip = 100 ##number of extant species
lam = 0.4
mu.vector <- runif(100, 0, 0.4)


##simulating 100 trees
Tree.1 <- vector("list", 100)

##for loop due to different extinction rates
for(i in seq_along(mu.vector)){
  Tree.1[i] <- sim.bd.taxa(n = Ntip, lambda = lam, mu = mu.vector[i],
                             complete = TRUE, numbsim = 1)
}

###simulating species taxonomy; only budding speciation

Taxa.1.list <- lapply(Tree.1, sim.taxonomy, beta = 0, lambda.a = 0)

##for matching posterior datasets
for(i in 1:length(Taxa.1.list)){
  Taxa.1.list[[i]]$sp <- sub("^", "t", Taxa.1.list[[i]]$sp)
}
                    
##Prunning the extinct species
Tree.1.extant <- lapply(Tree.1, prune.fossil.tips) 

##species might be associated with multiple edges;
##for calculating the evolutionary age we should pick the oldest edge

extant.sp.list <- vector("list", 100)
for(i in seq_along(extant.sp.list)){
  extant.sp.list[[i]] <- Taxa.1.list[[i]][Taxa.1.list[[i]]$sp 
                                 %in% Tree.1.extant[[i]]$tip.label, ] %>% 
                                 group_by(sp) %>%
                                 filter(start == max(start))%>%
                                 rename(true.age = start)
}

#####Calculating the phylo ages from the "calculate_tip_ages" function

phylo.list.age <- lapply(Tree.1.extant, calculate_tip_ages)

##joining each dataframe of both lists

discrepancies.df <- do.call("rbind",map2(extant.sp.list, phylo.list.age,
                                           left_join, by = "sp"))

##lambda
discrepancies.df$lambda <- rep(lam, nrow(discrepancies.df))

##mu
discrepancies.df$mu <- rep(mu.vector, each = 100)

##new variables
discrepancies.df.1 <- discrepancies.df %>% 
                      mutate(div = lambda - mu,
                             absolute_error = abs(phylo.age - true.age),
                             relative_error = (phylo.age-true.age)/true.age)

discrepancies.df.1$mu.int <- cut_interval(discrepancies.df.1$mu, 4)

discrepancies.df.1$div.int <- cut_interval(discrepancies.df.1$div, 4)

##root ages for rescaling the ages (using "tree.max" from FossilSim)

root.ages <- sapply(Tree.1, tree.max)
discrepancies.df.1$root.age <- rep(root.ages, each = 100) 

##escaling the ages with root age, and consequently the errors

discrepancies.df.1 <- discrepancies.df.1 %>% 
                    mutate(phylo.age.esc = phylo.age/root.age,
                    true.age.esc = true.age/root.age,
                    absolute_error_esc = abs(phylo.age.esc - true.age.esc),
                    relative_error_esc = (phylo.age.esc-
                                            true.age.esc)/true.age.esc) 
                    
##absolute error dataset
discrepancies.df.2 <- discrepancies.df.1 %>% group_by(mu) %>% 
                      summarize(mean_abs_error = mean(absolute_error),
                                mean_abs_esc_error = mean(absolute_error_esc)) 

  
##saving dataframes

write_xlsx(discrepancies.df.1, "dis.1.xlsx")
write_xlsx(discrepancies.df.2, "dis.2.xlsx")
