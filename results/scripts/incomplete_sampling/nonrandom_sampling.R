########nonrandom sampling and age-dependent extinction

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

##simulating trees with intermediate extinction for the the missing scenarios
## a constant 100 extant tips

## 0% missing species
trees.0.missing <- sim.bd.taxa(100, numbsim = 1000, lambda = 0.3, mu = 0.15)

## 25% missing species
trees.25.missing <- sim.bd.taxa(134, numbsim = 1000, lambda = 0.3, mu = 0.15)

## 50% missing species
trees.50.missing <- sim.bd.taxa(200, numbsim = 1000, lambda = 0.3, mu = 0.15)


##prunning trees

## 0% missing species
trees.0.extant <- lapply(trees.0.missing, prune.fossil.tips) 

##root ages
root.0 <- sapply(trees.0.missing, tree.max)

## 25% missing species
trees.25.extant <- lapply(trees.25.missing, prune.fossil.tips)

##root ages
root.25 <- sapply(trees.25.missing, tree.max)

### 50% missing species
trees.50.extant <- lapply(trees.50.missing, prune.fossil.tips)

## root ages
root.50 <- sapply(trees.50.missing, tree.max)


# Bifurcating speciation --------------------------------------------------

## 0 % missing species

##simulating speciation via bifurcation
tax.0.bif <- lapply(trees.0.missing, sim.taxonomy, beta = 1)

##get ages extant species
ages.0.bif <- Map(getAgesExtantSpecies, tax.0.bif, trees.0.missing,
                  Tol = 1e-6) 

## 25% missing species
tax.25.bif <- lapply(trees.25.missing, sim.taxonomy, beta = 1)

##get ages extant species
ages.25.bif <- Map(getAgesExtantSpecies, tax.25.bif, trees.25.missing,
                   Tol = 1e-6)

## 50% missing species
tax.50.bif <- lapply(trees.50.missing, sim.taxonomy, beta = 1)

## get ages extant species
ages.50.bif <- Map(getAgesExtantSpecies, tax.50.bif, trees.50.missing,
                   Tol = 1e-6)

###########nonrandom incomplethe sampling of extant species#######################

##the probability of not getting sample increases with age

##get ages extant species with a nonrandom incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bif, trees.25.missing,
                          SamplingFrac = 0.25,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


##get ages extant species with a nonrandom incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bif, trees.50.missing,
                          SamplingFrac = 0.50,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))


##0% missing species

#collapsing in a DF
ages.0.bif <- do.call("rbind", ages.0.bif) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.0.bif$root.age <- rep(root.0, each = 100)

##adding rTrue.age
ages.0.bif <- ages.0.bif %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.0.bif$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species
ages.0.bif$species <- paste0(ages.0.bif$label,".",
                             ages.0.bif$tree)

##assigning the positive ADE through the IUCN status

####intermediate extinction
ages.0.bif$status <- discretize(ages.0.bif$rTrue.age, 
                                method = "frequency", 
                                breaks = 5,
                                labels = c("LC", "NT", "VU", "EN", "CR"))

ages.0.bif$status <- factor(ages.0.bif$status,
                            levels = c("LC", "NT", "VU", "EN", "CR"),
                            ordered = TRUE)


##25% missing species

ages.25.bif <- do.call("rbind", ages.25.bif) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.25.bif$root.age <- rep(root.25, each = 134)

##adding rTrue.age
ages.25.bif <- ages.25.bif %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.25.bif$tree <- rep(paste0("tree.", 1:1000), each = 134)

##species
ages.25.bif$species <- paste0(ages.25.bif$label,".",
                              ages.25.bif$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.25.bif$status <- discretize(ages.25.bif$rTrue.age, 
                                 method = "frequency", 
                                 breaks = 5,
                                 labels = c("LC", "NT", "VU", "EN", "CR"))

ages.25.bif$status <- factor(ages.25.bif$status,
                             levels = c("LC", "NT", "VU", "EN", "CR"))


##50% missing species
ages.50.bif <- do.call("rbind", ages.50.bif) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.50.bif$root.age <- rep(root.0, each = 200)

##adding rTrue.age
ages.50.bif <- ages.50.bif %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.50.bif$tree <- rep(paste0("tree.", 1:1000), each = 200)

##species
ages.50.bif$species <- paste0(ages.50.bif$label,".",
                              ages.50.bif$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.50.bif$status <- discretize(ages.50.bif$rTrue.age, 
                                 method = "frequency", 
                                 breaks = 5,
                                 labels = c("LC", "NT", "VU", "EN", "CR"))

ages.50.bif$status <- factor(ages.50.bif$status,
                             levels = c("LC", "NT", "VU", "EN", "CR"))


##selecting ages.bif.0.25 columns
ages.bif.incomp.25 <- ages.bif.incomp.25 %>% select(species,  
                                                    Incomp.age) %>% 
  rename(Incomp.age.25 = Incomp.age)

##merging ages.bif with ages.bif.0.25
ages.total.0.25 <- left_join(ages.bif.incomp.25, ages.25.bif, 
                             by = "species")

##selecting ages.bif.0.5 columns
ages.bif.incomp.50 <- ages.bif.incomp.50 %>% select(species,
                                                    Incomp.age) %>% 
  rename(Incomp.age.50 = Incomp.age)

##merging the three dataframes
ages.total.0.50.bif <- left_join(ages.bif.incomp.50, ages.50.bif,
                                 by = "species")



# Correcting ages ---------------------------------------------------------

###Previous probabilistic approach

##fully sampled trees
full.list <- list()

for(i in 1:nrow(ages.0.bif)){
  full.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                        mu = 0.15,
                                        node_age = ages.0.bif$Estimated.age[i],
                                        rho = 1)
}

##mean age
full.mean.age <- sapply(full.list, function(x) sum(x$time * x$prob))

##median age
full.median.age <- sapply(full.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))


###testing new funciton provided by the reviewer

#mean

full.mean.new <- vector(mode = "numeric",
                           length = nrow(ages.0.bif))

for(i in 1:nrow(ages.0.bif)){
  full.mean.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0.15,
                                           rho = 1,
                                           v = ages.0.bif$Estimated.age[i])
}

##median
full.median.new <- vector(mode = "numeric",
                        length = nrow(ages.0.bif))


for(i in 1:nrow(ages.0.bif)){
  full.median.new[[i]] <- medianAgeFunction(lam = 0.3,
                                               mu = 0.15,
                                               rho = 1,
                                               v = ages.0.bif$Estimated.age[i])
}



##binding
ages.total.full <- cbind(ages.0.bif, full.mean.age, 
                         full.median.age, full.mean.new,
                         full.median.new) %>% 
  rename(mean.age = full.mean.age,
         median.age = full.median.age,
         mean.new = full.mean.new,
         median.new = full.median.new)

##saving
write_csv(ages.total.full,
          file = "results/data/processed/conservation_status/ages.full.sampled.csv")




#######25% missing species

#previous function

f25.list.non <- list()

for(i in 1:nrow(ages.total.0.25)){
  f25.list.non[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                           mu = 0.15,
                                           node_age = ages.total.0.25$Incomp.age.25[i],
                                           rho = 0.75)
}

##mean age
f25.mean.age <- sapply(f25.list.non, function(x) sum(x$time * x$prob))

##median age
f25.median.age <- sapply(f25.list.non, function(x)
  weighted.median(x = x$time,
                  w = x$prob))

###testing new funciton provided by the reviewer

#mean

f25.mean.new <- vector(mode = "numeric",
                           length = nrow(ages.total.0.25))

for(i in 1:nrow(ages.total.0.25)){
  f25.mean.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0.15,
                                           rho = 0.75,
                                           v = ages.total.0.25$Incomp.age.25[i])
}

##median
f25.median.new <- vector(mode = "numeric",
                             length = nrow(ages.total.0.25))

for(i in 1:nrow(ages.total.0.25)){
  f25.median.new[[i]] <- medianAgeFunction(lam = 0.3,
                                               mu = 0.15,
                                               rho = 0.75,
                                               v = ages.total.0.25$Incomp.age.25[i])
}



##binding
ages.total.f25 <- cbind(ages.total.0.25, f25.mean.age, 
                        f25.median.age,
                        f25.mean.new, f25.median.new) %>% 
  rename(mean.age = f25.mean.age,
         median.age = f25.median.age,
         mean.new = f25.mean.new,
         median.new = f25.median.new)

##saving
write_csv(ages.total.f25,
          file = "results/data/processed/conservation_status/ages.f25.sampled.csv")

ages.total.f25 <- read_csv("results/data/processed/conservation_status/ages.f25.sampled.csv")

###50% missing species

#previous function

f50.list.non <- list()

for(i in 1:nrow(ages.total.0.50.bif)){
  f50.list.non[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                           mu = 0.15,
                                           node_age = ages.total.0.50.bif$Incomp.age.50[i],
                                           rho = 0.5)
}

##mean age
f50.mean.age <- sapply(f50.list.non, function(x) sum(x$time * x$prob))

##median age
f50.median.age <- sapply(f50.list.non, function(x)
  weighted.median(x = x$time,
                  w = x$prob))


###testing new funciton provided by the reviewer

#mean

f50.mean.new <- vector(mode = "numeric",
                       length = nrow(ages.total.0.50.bif))

for(i in 1:nrow(ages.total.0.50.bif)){
  f50.mean.new[[i]] <- meanAgeFunction(lam = 0.3,
                                       mu = 0.15,
                                       rho = 0.50,
                                       v = ages.total.0.50.bif$Incomp.age.50[i])
}

##median
f50.median.new <- vector(mode = "numeric",
                         length = nrow(ages.total.0.50.bif))

for(i in 1:nrow(ages.total.0.50.bif)){
  f50.median.new[[i]] <- medianAgeFunction(lam = 0.3,
                                           mu = 0.15,
                                           rho = 0.50,
                                           v = ages.total.0.50.bif$Incomp.age.50[i])
}


##binding
ages.total.f50 <- cbind(ages.total.0.50.bif, f50.mean.age, 
                        f50.median.age,
                        f50.mean.new, f50.median.new) %>% 
  rename(mean.age = f50.mean.age,
         median.age = f50.median.age,
         mean.new = f50.mean.new,
         median.new = f50.median.new)

##saving
write_csv(ages.total.f50,
          file = "results/data/processed/conservation_status/ages.f50.sampled.csv")

ages.total.f50 <- read_csv("results/data/processed/conservation_status/ages.f50.sampled.csv")


##Calculating the mean of each age for each category within each tree

#fully sampled
ages.mean.full <- ages.total.full %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.full = mean(Estimated.age, na.rm = TRUE),
            mean.mean = mean(mean.age, na.rm = TRUE),
            mean.median = mean(median.age, na.rm = TRUE),
            mean.mean.new = mean(mean.new, na.rm = TRUE),
            mean.median.new = mean(median.new, na.rm = TRUE))

##25% unsampled
ages.mean.f25 <- ages.total.f25 %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.f25 = mean(Incomp.age.25, na.rm = TRUE),
            mean.mean = mean(mean.age, na.rm = TRUE),
            mean.median = mean(median.age, na.rm = TRUE),
            mean.mean.new = mean(mean.new, na.rm = TRUE),
            mean.median.new = mean(median.new, na.rm = TRUE))

#50% unsampled
ages.mean.f50 <- ages.total.f50 %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.f50 = mean(Incomp.age.50, na.rm = TRUE),
            mean.mean = mean(mean.age, na.rm = TRUE),
            mean.median = mean(median.age, na.rm = TRUE),
            mean.mean.new = mean(mean.new, na.rm = TRUE),
            mean.median.new = mean(median.new, na.rm = TRUE))


####ages ranking and evaluating wrong estimations

#fully sampled
ages.rank.full <- ages.mean.full %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         full.rank = dense_rank(mean.full),
         mean.rank = dense_rank(mean.mean),
         median.rank = dense_rank(mean.median),
         mean.new.rank = dense_rank(mean.mean.new),
         median.new.rank = dense_rank(mean.median.new)) %>% 
  mutate(resp_full = if_else(true.rank == full.rank, 1, 0),
         resp_mean = if_else(true.rank == mean.rank, 1, 0),
         resp_median = if_else(true.rank == median.rank, 1,0),
         resp_mean_new = if_else(true.rank == mean.new.rank, 1, 0),
         resp_median_new = if_else(true.rank == median.new.rank, 1, 0))


#25% unsampled
ages.rank.f25 <- ages.mean.f25 %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         f25.rank = dense_rank(mean.f25),
         mean.rank = dense_rank(mean.mean),
         median.rank = dense_rank(mean.median),
         mean.new.rank = dense_rank(mean.mean.new),
         median.new.rank = dense_rank(mean.median.new)) %>% 
  mutate(resp_f25 = if_else(true.rank == f25.rank, 1, 0),
         resp_mean = if_else(true.rank == mean.rank, 1, 0),
         resp_median = if_else(true.rank == median.rank, 1,0),
         resp_mean_new = if_else(true.rank == mean.new.rank, 1, 0),
         resp_median_new = if_else(true.rank == median.new.rank, 1, 0)) 

#50% unsampled
ages.rank.f50 <- ages.mean.f50 %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         f50.rank = dense_rank(mean.f50),
         mean.rank = dense_rank(mean.mean),
         median.rank = dense_rank(mean.median),
         mean.new.rank = dense_rank(mean.mean.new),
         median.new.rank = dense_rank(mean.median.new)) %>% 
  mutate(resp_f50 = if_else(true.rank == f50.rank, 1, 0),
         resp_mean = if_else(true.rank == mean.rank, 1, 0),
         resp_median = if_else(true.rank == median.rank, 1,0),
         resp_mean_new = if_else(true.rank == mean.new.rank, 1, 0),
         resp_median_new = if_else(true.rank == median.new.rank, 1, 0)) 

##comparing categories and defining error rates

####fully sampled

#phylogenetic
ages.comparison.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_full = sum(resp_full, na.rm = TRUE)) %>% 
  filter(sum_full < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)


#Mean age previous function
ages.comparison.mean.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#median age previous function
ages.comparison.median.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

###Mean age NEW function
ages.comparison.mean.new.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_mean_new = sum(resp_mean_new, na.rm = TRUE)) %>% 
  filter(sum_mean_new < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#median age NEW function
ages.comparison.median.new.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_median_new = sum(resp_median_new, na.rm = TRUE)) %>% 
  filter(sum_median_new < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)


##25% unsampled

#phylogenetic
ages.comparison.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25, na.rm = TRUE)) %>% 
  filter(sum_f25 < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)


#mean age previous function
ages.comparison.mean.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#median age previous function
ages.comparison.median.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#Mean age NEW function
ages.comparison.mean.new.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_mean_new = sum(resp_mean_new, na.rm = TRUE)) %>% 
  filter(sum_mean_new < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#Median age NEW function
ages.comparison.median.new.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_median_new = sum(resp_median_new, na.rm = TRUE)) %>% 
  filter(sum_median_new < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

###50% unsampled

#phylogenetic
ages.comparison.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)


#Mean age previous function
ages.comparison.mean.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#Median age previous function
ages.comparison.median.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#Mean age NEW function
ages.comparison.mean.new.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_mean_new = sum(resp_mean_new, na.rm = TRUE)) %>% 
  filter(sum_mean_new < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)

#Median age NEW function
ages.comparison.median.new.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_median_new = sum(resp_median_new, na.rm = TRUE)) %>% 
  filter(sum_median_new < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)


###Selecting trees

##full sampled

##correct estimations 
tree.correct.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_full = sum(resp_full)) %>% 
  filter(sum_full == 5) %>% 
  pull(tree)

tree.correct.full <- sample(tree.correct.full,
                            size = round(length(tree.correct.full)*0.1))

##incorrect estimations
tree.incorrect.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_full = sum(resp_full, na.rm = TRUE)) %>% 
  filter(sum_full < 5) %>% 
  pull(tree)

tree.incorrect.full <- sample(tree.incorrect.full,
                              size = round(length(tree.incorrect.full)*0.1))

##mean prob
##correct estimations 
tree.correct.mean.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean)) %>% 
  filter(sum_mean == 5) %>% 
  pull(tree)

tree.correct.mean.full <- sample(tree.correct.mean.full,
                                 size = round(length(tree.correct.mean.full)*0.1))

##incorrect estimations
tree.incorrect.mean.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.mean.full <- sample(tree.incorrect.mean.full,
                                   size = round(length(tree.incorrect.mean.full)*0.1))


##median prob
##correct estimations 
tree.correct.median.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median)) %>% 
  filter(sum_median == 5) %>% 
  pull(tree)

tree.correct.median.full <- sample(tree.correct.median.full,
                                   size = round(length(tree.correct.median.full)*0.1))

##incorrect estimations
tree.incorrect.median.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.median.full <- sample(tree.incorrect.median.full,
                                     size = round(length(tree.incorrect.median.full)*0.1))


####new function
##mean prob
##correct estimations 
tree.correct.mean.new.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_mean_new = sum(resp_mean_new)) %>% 
  filter(sum_mean_new == 5) %>% 
  pull(tree)

tree.correct.mean.new.full <- sample(tree.correct.mean.new.full,
                                 size = round(length(tree.correct.mean.new.full)*0.1))

##incorrect estimations
tree.incorrect.mean.new.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_mean_new = sum(resp_mean_new, na.rm = TRUE)) %>% 
  filter(sum_mean_new < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.mean.new.full <- sample(tree.incorrect.mean.new.full,
                                   size = round(length(tree.incorrect.mean.new.full)*0.1))


##median prob
##correct estimations 
tree.correct.median.new.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median_new)) %>% 
  filter(sum_median == 5) %>% 
  pull(tree)

tree.correct.median.new.full <- sample(tree.correct.median.new.full,
                                   size = round(length(tree.correct.median.new.full)*0.1))

##incorrect estimations
tree.incorrect.median.new.full <- ages.rank.full %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median_new, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.median.new.full <- sample(tree.incorrect.median.new.full,
                                     size = round(length(tree.incorrect.median.new.full)*0.1))


###25% missing species

##correct estimations 
tree.correct.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25)) %>% 
  filter(sum_f25 == 5) %>% 
  pull(tree)

tree.correct.f25 <- sample(tree.correct.f25,
                           size = round(length(tree.correct.f25)*0.1))

##incorrect estimations
tree.incorrect.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25, na.rm = TRUE)) %>% 
  filter(sum_f25 < 5) %>% 
  pull(tree)

tree.incorrect.f25 <- sample(tree.incorrect.f25,
                             size = round(length(tree.incorrect.f25)*0.1))

##mean prob
##correct estimations 
tree.correct.mean.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean)) %>% 
  filter(sum_mean == 5) %>% 
  pull(tree)

tree.correct.mean.f25 <- sample(tree.correct.mean.f25,
                                size = round(length(tree.correct.mean.f25)*0.1))

##incorrect estimations
tree.incorrect.mean.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.mean.f25 <- sample(tree.incorrect.mean.f25,
                                  size = round(length(tree.incorrect.mean.f25)*0.1))


##median prob
##correct estimations 
tree.correct.median.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median)) %>% 
  filter(sum_median == 5) %>% 
  pull(tree)

tree.correct.median.f25 <- sample(tree.correct.median.f25,
                                  size = round(length(tree.correct.median.f25)*0.1))

##incorrect estimations
tree.incorrect.median.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.median.f25 <- sample(tree.incorrect.median.f25,
                                    size = round(length(tree.incorrect.median.f25)*0.1))



##NEW mean prob
##correct estimations 
tree.correct.mean.new.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean_new)) %>% 
  filter(sum_mean == 5) %>% 
  pull(tree)

tree.correct.mean.new.f25 <- sample(tree.correct.mean.new.f25,
                                size = round(length(tree.correct.mean.new.f25)*0.1))

##incorrect estimations
tree.incorrect.mean.new.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean_new, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.mean.new.f25 <- sample(tree.incorrect.mean.new.f25,
                                  size = round(length(tree.incorrect.mean.new.f25)*0.1))


##median prob
##correct estimations 
tree.correct.median.new.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median_new)) %>% 
  filter(sum_median == 5) %>% 
  pull(tree)

tree.correct.median.new.f25 <- sample(tree.correct.median.new.f25,
                                  size = round(length(tree.correct.median.new.f25)*0.1))

##incorrect estimations
tree.incorrect.median.new.f25 <- ages.rank.f25 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median_new, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.median.new.f25 <- sample(tree.incorrect.median.new.f25,
                                    size = round(length(tree.incorrect.median.new.f25)*0.1))


## 50% missing species

##correct estimations 
tree.correct.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50)) %>% 
  filter(sum_f50 == 5) %>% 
  pull(tree)

tree.correct.f50 <- sample(tree.correct.f50,
                           size = round(length(tree.correct.f50)*0.1))

##incorrect estimations
tree.incorrect.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 < 5) %>% 
  pull(tree)

tree.incorrect.f50 <- sample(tree.incorrect.f50,
                             size = round(length(tree.incorrect.f50)*0.1))

##mean prob
##correct estimations 
tree.correct.mean.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean)) %>% 
  filter(sum_mean == 5) %>% 
  pull(tree)

tree.correct.mean.f50 <- sample(tree.correct.mean.f50,
                                size = round(length(tree.correct.mean.f50)*0.1))

##incorrect estimations
tree.incorrect.mean.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.mean.f50 <- sample(tree.incorrect.mean.f50,
                                  size = round(length(tree.incorrect.mean.f50)*0.1))


##median prob
##correct estimations 
tree.correct.median.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median)) %>% 
  filter(sum_median == 5) %>% 
  pull(tree)

tree.correct.median.f50 <- sample(tree.correct.median.f50,
                                  size = round(length(tree.correct.median.f50)*0.1))

##incorrect estimations
tree.incorrect.median.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.median.f50 <- sample(tree.incorrect.median.f50,
                                    size = round(length(tree.incorrect.median.f50)*0.1))

##mean NEW prob
##correct estimations 
tree.correct.mean.new.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean_new)) %>% 
  filter(sum_mean == 5) %>% 
  pull(tree)

tree.correct.mean.new.f50 <- sample(tree.correct.mean.new.f50,
                                size = round(length(tree.correct.mean.new.f50)*0.1))

##incorrect estimations
tree.incorrect.mean.new.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean_new, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.mean.new.f50 <- sample(tree.incorrect.mean.new.f50,
                                  size = round(length(tree.incorrect.mean.new.f50)*0.1))


##median prob
##correct estimations 
tree.correct.median.new.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median_new)) %>% 
  filter(sum_median == 5) %>% 
  pull(tree)

tree.correct.median.new.f50 <- sample(tree.correct.median.new.f50,
                                  size = round(length(tree.correct.median.new.f50)*0.1))

##incorrect estimations
tree.incorrect.median.new.f50 <- ages.rank.f50 %>% group_by(tree) %>% 
  summarise(sum_median = sum(resp_median_new, na.rm = TRUE)) %>% 
  filter(sum_median < 5) %>% 
  drop_na() %>% 
  pull(tree)

tree.incorrect.median.new.f50 <- sample(tree.incorrect.median.new.f50,
                                    size = round(length(tree.incorrect.median.new.f50)*0.1))


##extinction categories
ages.mean.full$status <- factor(ages.mean.full$status,
                                levels = c("LC", "NT", "VU", "EN", "CR"),
                                ordered= TRUE)

ages.mean.f25$status <- factor(ages.mean.f25$status,
                               levels = c("LC", "NT", "VU", "EN", "CR"),
                               ordered= TRUE)

ages.mean.f50$status <- factor(ages.mean.f50$status,
                               levels = c("LC", "NT", "VU", "EN", "CR"),
                               ordered= TRUE)



##false facets
ages.mean.full$full <- "Fully sampled"
ages.mean.f25$f25 <- "25% missing species"
ages.mean.f50$f50 <- "50% missing species"



# Plotting ----------------------------------------------------------------




##########plots#########################################

## themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (10),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (11),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (10)),
                      legend.background = element_rect(fill = "gray90",
                                                       size = 0.5, linetype = "dotted"),
                      legend.position = "bottom")


##Fully sampled

##correct
ages.correct.full <- ages.mean.full %>% filter(tree %in% tree.correct.full)

##incorrect
ages.incorrect.full <- ages.mean.full %>% 
  filter(tree %in% tree.incorrect.full)

full.plot <- ggplot(ages.correct.full,
                    aes(x = status, y = log(mean.full +1), group = tree,
                        color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.3, color = "#d95f02")+
  geom_point(data = ages.incorrect.full,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.incorrect.full,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.full, aes(x = 2.0, y = 2.0,
                                             label = paste0("Error rate:", " ",error,
                                                            "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab("Phylogenetic age (log)")+
  xlab(NULL)+
  facet_wrap(~full)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x =element_blank())

##mean datasets
ages.correct.mean.full <- ages.mean.full %>% filter(tree %in%
                                                      tree.correct.mean.full)

##incorrect
ages.incorrect.mean.full <- ages.mean.full %>% 
  filter(tree %in% tree.incorrect.mean.full)

##mean
full.mean.plot<- ggplot(ages.correct.mean.full,
                        aes(x = status, y = log(mean.mean +1), group = tree,
                            color = "#7fc97f"))+
  geom_point(size=3, shape = 21, fill="white", color = "#7fc97f")+
  geom_line(alpha = 0.5, color = "#7fc97f")+
  geom_point(data = ages.incorrect.mean.full,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.mean.full,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.mean.full, aes(x = 1.7, y = 1.3,
                                                  label = paste0("Error rate:", " ",error,
                                                                 "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab("Mean age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##median dataset
ages.correct.median.full <- ages.mean.full %>% filter(tree %in%
                                                        tree.correct.median.full)

##incorrect
ages.incorrect.median.full <- ages.mean.full %>% 
  filter(tree %in% tree.incorrect.median.full)

##median
full.median.plot<- ggplot(ages.correct.median.full,
                          aes(x = status, y = log(mean.median +1), group = tree,
                              color = "#beaed4"))+
  geom_point(size=3, shape = 21, fill="white", color = "#beaed4")+
  geom_line(alpha = 0.5, color = "#beaed4")+
  geom_point(data = ages.incorrect.median.full,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.median.full,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.median.full, aes(x = 1.7, y = 1.25,
                                                    label = paste0("Error rate:", " ",error,
                                                                   "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab("Median age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")


##mean NEW datasets
ages.correct.mean.new.full <- ages.mean.full %>% filter(tree %in%
                                                      tree.correct.mean.new.full)

##incorrect
ages.incorrect.mean.new.full <- ages.mean.full %>% 
                               filter(tree %in% tree.incorrect.mean.new.full)

##mean
full.mean.new.plot<- ggplot(ages.correct.mean.new.full,
                        aes(x = status, y = log(mean.mean.new +1), group = tree,
                            color = "blue"))+
  geom_point(size=3, shape = 21, fill="white", color = "blue")+
  geom_line(alpha = 0.5, color = "blue")+
  geom_point(data = ages.incorrect.mean.new.full,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.mean.new.full,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.mean.new.full, aes(x = 1.7, y = 1.3,
                                   label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab("Mean New age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##median dataset
ages.correct.median.new.full <- ages.mean.full %>% filter(tree %in%
                                                        tree.correct.median.new.full)

##incorrect
ages.incorrect.median.new.full <- ages.mean.full %>% 
  filter(tree %in% tree.incorrect.median.new.full)

##median
full.median.new.plot<- ggplot(ages.correct.median.new.full,
                          aes(x = status, y = log(mean.median.new +1), group = tree,
                              color = "skyblue"))+
  geom_point(size=3, shape = 21, fill="white", color = "skyblue")+
  geom_line(alpha = 0.5, color = "skyblue")+
  geom_point(data = ages.incorrect.median.new.full,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.median.new.full,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.median.new.full, aes(x = 1.7, y = 1.25,
                                                    label = paste0("Error rate:", " ",error,
                                                                   "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab("Median New age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")



#25% missing species

##phylogenetic age datasets

##correct
ages.correct.f25 <- ages.mean.f25 %>% filter(tree %in% tree.correct.f25)

##incorrect
ages.incorrect.f25 <- ages.mean.f25 %>% 
  filter(tree %in% tree.incorrect.f25)

##phylogenetic
f25.plot <- ggplot(ages.correct.f25,
                   aes(x = status, y = log(mean.f25 +1), group = tree,
                       color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.3, color = "#d95f02")+
  geom_point(data = ages.incorrect.f25,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.incorrect.f25,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f25, aes(x = 2.0, y = 2.1,
                                            label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f25)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##mean datasets
ages.correct.mean.f25 <- ages.mean.f25 %>% filter(tree %in%
                                                    tree.correct.mean.f25)

##incorrect
ages.incorrect.mean.f25 <- ages.mean.f25 %>% 
  filter(tree %in% tree.incorrect.mean.f25)

##mean
f25.mean.plot<- ggplot(ages.correct.mean.f25,
                       aes(x = status, y = log(mean.mean +1), group = tree,
                           color = "#7fc97f"))+
  geom_point(size=3, shape = 21, fill="white", color = "#7fc97f")+
  geom_line(alpha = 0.5, color = "#7fc97f")+
  geom_point(data = ages.incorrect.mean.f25,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.mean.f25,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.mean.f25, aes(x = 1.8, y = 1.2,
                                                 label = paste0("Error rate:", " ",error,
                                                                "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##median dataset
ages.correct.median.f25 <- ages.mean.f25 %>% filter(tree %in%
                                                      tree.correct.median.f25)

##incorrect
ages.incorrect.median.f25 <- ages.mean.f25 %>% 
  filter(tree %in% tree.incorrect.median.f25)

##Median
f25.median.plot<- ggplot(ages.correct.median.f25,
                         aes(x = status, y = log(mean.median +1), group = tree,
                             color = "#beaed4"))+
  geom_point(size=3, shape = 21, fill="white", color = "#beaed4")+
  geom_line(alpha = 0.5, color = "#beaed4")+
  geom_point(data = ages.incorrect.median.f25,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.median.f25,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.median.f25, aes(x = 1.72, y = 1.1,
                                                   label = paste0("Error rate:", " ",error,
                                                                  "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")

##NEW mean datasets
ages.correct.mean.new.f25 <- ages.mean.f25 %>% filter(tree %in%
                                                    tree.correct.mean.new.f25)

##incorrect
ages.incorrect.mean.new.f25 <- ages.mean.f25 %>% 
  filter(tree %in% tree.incorrect.mean.new.f25)

##mean
f25.mean.new.plot<- ggplot(ages.correct.mean.new.f25,
                       aes(x = status, y = log(mean.mean.new +1), group = tree,
                           color = "blue"))+
  geom_point(size=3, shape = 21, fill="white", color = "blue")+
  geom_line(alpha = 0.5, color = "blue")+
  geom_point(data = ages.incorrect.mean.new.f25,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "gray")+
  geom_line(data = ages.incorrect.mean.new.f25,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.mean.new.f25, aes(x = 1.8, y = 1.2,
                                                 label = paste0("Error rate:", " ",error,
                                                                "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##median dataset
ages.correct.median.new.f25 <- ages.mean.f25 %>% filter(tree %in%
                                                      tree.correct.median.new.f25)

##incorrect
ages.incorrect.median.new.f25 <- ages.mean.f25 %>%  filter(tree %in% tree.incorrect.median.new.f25)

##Median
f25.median.new.plot<- ggplot(ages.correct.median.new.f25,
                         aes(x = status, y = log(mean.median.new +1), group = tree,
                             color = "skyblue"))+
  geom_point(size=3, shape = 21, fill="white", color = "skyblue")+
  geom_line(alpha = 0.6, color = "skyblue")+
  geom_point(data = ages.incorrect.median.new.f25,
             size=3, shape = 21, fill="white",
             alpha = 0.3, color = "#969696")+
  geom_line(data = ages.incorrect.median.new.f25,
            alpha = 0.3, color = "#969696")+
  geom_text(data = ages.comparison.median.new.f25, aes(x = 1.72, y = 1.1,
                                                   label = paste0("Error rate:", " ",error,
                                                                  "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")



##50% missing species

##phylogenetic age datasets

##correct
ages.correct.f50 <- ages.mean.f50 %>% filter(tree %in% tree.correct.f50)

##incorrect
ages.incorrect.f50 <- ages.mean.f50 %>% 
  filter(tree %in% tree.incorrect.f50)

#Phylogenetic
f50.plot <- ggplot(ages.correct.f50,
                   aes(x = status, y = log(mean.f50 +1), group = tree,
                       color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.5, color = "#d95f02")+
  geom_point(data = ages.incorrect.f50,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.incorrect.f50,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f50, aes(x = 1.62, y = 2.1,
                                            label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f50)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##mean datasets
ages.correct.mean.f50 <- ages.mean.f50 %>% filter(tree %in%
                                                    tree.correct.mean.f50)

##incorrect
ages.incorrect.mean.f50 <- ages.mean.f50 %>% 
  filter(tree %in% tree.incorrect.mean.f50)

##Mean
f50.mean.plot<- ggplot(ages.correct.mean.f50,
                       aes(x = status, y = log(mean.mean +1), group = tree,
                           color = "#7fc97f"))+
  geom_point(size=3, shape = 21, fill="white", color = "#7fc97f")+
  geom_line(alpha = 0.5, color = "#7fc97f")+
  geom_point(data = ages.incorrect.mean.f50,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.mean.f50,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.mean.f50, aes(x = 1.65, y = 1.1,
                                                 label = paste0("Error rate:", " ",error,
                                                                "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##median dataset
ages.correct.median.f50 <- ages.mean.f50 %>% filter(tree %in%
                                                      tree.correct.median.f50)

##incorrect
ages.incorrect.median.f50 <- ages.mean.f50 %>% 
  filter(tree %in% tree.incorrect.median.f50)

##median
f50.median.plot<- ggplot(ages.correct.median.f50,
                         aes(x = status, y = log(mean.median +1), group = tree,
                             color = "#beaed4"))+
  geom_point(size=3, shape = 21, fill="white", color = "#beaed4")+
  geom_line(alpha = 0.3, color = "#beaed4")+
  geom_point(data = ages.incorrect.median.f50,
             size=3, shape = 21, fill="white",
             alpha = 0.6, color = "#969696")+
  geom_line(data = ages.incorrect.median.f50,
            alpha = 0.6, color = "#969696")+
  geom_text(data = ages.comparison.median.f50, aes(x = 1.7, y = 0.95,
                                                   label = paste0("Error rate:", " ",error,
                                                                  "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")


##NEW mean datasets
ages.correct.mean.new.f50 <- ages.mean.f50 %>% filter(tree %in%
                                                    tree.correct.mean.new.f50)

##incorrect
ages.incorrect.mean.new.f50 <- ages.mean.f50 %>% 
  filter(tree %in% tree.incorrect.mean.new.f50)

##Mean
f50.mean.new.plot<- ggplot(ages.correct.mean.new.f50,
                       aes(x = status, y = log(mean.mean.new +1), group = tree,
                           color = "blue"))+
  geom_point(size=3, shape = 21, fill="white", color = "blue")+
  geom_line(alpha = 0.5, color = "blue")+
  geom_point(data = ages.incorrect.mean.new.f50,
             size=3, shape = 21, fill="white",
             alpha = 0.3, color = "#969696")+
  geom_line(data = ages.incorrect.mean.new.f50,
            alpha = 0.3, color = "#969696")+
  geom_text(data = ages.comparison.mean.new.f50, aes(x = 1.65, y = 1.1,
                                                 label = paste0("Error rate:", " ",error,
                                                                "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none",
        axis.text.x = element_blank())

##median dataset
ages.correct.median.new.f50 <- ages.mean.f50 %>% filter(tree %in%
                                                      tree.correct.median.new.f50)

##incorrect
ages.incorrect.median.new.f50 <- ages.mean.f50 %>% 
  filter(tree %in% tree.incorrect.median.new.f50)

##median
f50.median.new.plot<- ggplot(ages.correct.median.new.f50,
                         aes(x = status, y = log(mean.median.new +1), group = tree,
                             color = "skyblue"))+
  geom_point(size=3, shape = 21, fill="white", color = "skyblue")+
  geom_line(alpha = 0.6, color = "skyblue")+
  geom_point(data = ages.incorrect.median.new.f50,
             size=3, shape = 21, fill="white",
             alpha = 0.3, color = "#969696")+
  geom_line(data = ages.incorrect.median.new.f50,
            alpha = 0.3, color = "#969696")+
  geom_text(data = ages.comparison.median.new.f50, aes(x = 1.7, y = 0.95,
                                   label = paste0("Error rate:", " ",error,
                                                                  "%")),
            inherit.aes = FALSE,
            size = 3.2, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")



###unifying plots

###png object
png("text/figures/age_dependent_samp.png", 
    width = 25, height = 22,
    units = "cm", 
    pointsize = 8, res = 300)

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
x.full = grid.arrange(full.plot, full.mean.plot, full.median.plot,
                      full.mean.new.plot, full.median.new.plot, nrow = 5)
x.f25 = grid.arrange(f25.plot, f25.mean.plot, f25.median.plot,
                     f25.mean.new.plot, f25.median.new.plot, nrow = 5)
x.f50 = grid.arrange(f50.plot, f50.mean.plot, f50.median.plot,
                     f50.mean.new.plot, f50.median.new.plot, nrow = 5)

grid.arrange(x.full, x.f25, x.f50, ncol = 3, bottom = bottom)
dev.off()
