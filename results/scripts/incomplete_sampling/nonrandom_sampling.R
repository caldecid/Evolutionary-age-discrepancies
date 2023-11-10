########nonrandom sampling and extinction

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

###########nonrandom incomple sampling of extant species#######################

##the nonrandom probability increases with age

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
                            levels = c("LC", "NT", "VU", "EN", "CR"))


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

##Calculating the mean of each age for each tree for 0% missing sp
ages.mean.0.bif <-ages.0.bif %>% 
                  group_by(tree, status) %>% 
                  summarise(mean.true = mean(True.age),
                            mean.phy = mean(Estimated.age))
write_csv(ages.mean.0.bif,
      file = "results/data/processed/incomplete_sampling/ages.mean.0.bif.csv")

ages.mean.0.bif <- read_csv(file = "results/data/processed/incomplete_sampling/ages.mean.0.bif.csv")
##Calculating the mean of each age for each tree for 25% missing sp                        
ages.mean.25.bif <- ages.total.0.25 %>% 
                        group_by(tree, status) %>% 
                        summarise(mean.true = mean(True.age),
                        mean.0.25 = mean(Incomp.age.25, na.rm = TRUE))

write_csv(ages.mean.25.bif,
          file = "results/data/processed/incomplete_sampling/ages.mean.25.bif.csv")

ages.mean.25.bif <- read_csv(file = "results/data/processed/incomplete_sampling/ages.mean.25.bif.csv")

####calculating the mean of each age for each tree for 50% missing sp
ages.mean.50.bif <- ages.total.0.50.bif %>% 
                    group_by(tree, status) %>% 
                    summarise(mean.true = mean(True.age),
                  mean.0.50 = mean(Incomp.age.50, na.rm = TRUE)) 

write_csv(ages.mean.50.bif,
        file = "results/data/processed/incomplete_sampling/ages.mean.50.bif.csv")

ages.mean.50.bif <- read_csv(file = "results/data/processed/incomplete_sampling/ages.mean.50.bif.csv")

#######ages ranking 0%
ages.rank.0.bif <- ages.mean.0.bif %>% group_by(tree) %>%
              mutate(true.rank = dense_rank(mean.true),
                     phy.rank = dense_rank(mean.phy)) %>% 
                     #f25.rank = dense_rank(mean.0.25),
                     #f50.rank = dense_rank(mean.0.50)) %>% 
              mutate(resp_phy = if_else(true.rank == phy.rank, 1, 0))
                     #resp_f25 = if_else(true.rank == f25.rank, 1, 0),
                     #resp_f50 = if_else(true.rank == f50.rank, 1, 0))

###comparing phylogenetic age
ages.comparison.f0.bif <- ages.rank.0.bif %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

## ages rank 25%
ages.rank.25.bif <- ages.mean.25.bif %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         f25.rank = dense_rank(mean.0.25)) %>% 
  #f50.rank = dense_rank(mean.0.50)) %>% 
  mutate(resp_f25 = if_else(true.rank == f25.rank, 1, 0))
#resp_f50 = if_else(true.rank == f50.rank, 1, 0))

##comparing incomplete sampling 25%
ages.comparison.f25.bif <- ages.rank.25.bif %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25, na.rm = TRUE)) %>% 
  filter(sum_f25 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

## ages rank 50%
ages.rank.50.bif <- ages.mean.50.bif %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         f50.rank = dense_rank(mean.0.50)) %>% 
  mutate(resp_f50 = if_else(true.rank == f50.rank, 1, 0))

##comparing incomplete sampling 50%
ages.comparison.f50.bif <- ages.rank.50.bif %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)


# budding speciation ------------------------------------------------------

## 0 % missing species

##simulating speciation via budding
tax.0.bud <- lapply(trees.0.missing, sim.taxonomy, beta = 0)

##get ages extant species
ages.0.bud <- Map(getAgesExtantSpecies, tax.0.bud, trees.0.missing,
                  Tol = 1e-6) 

## 25% missing species
tax.25.bud <- lapply(trees.25.missing, sim.taxonomy, beta = 0)

##get ages extant species
ages.25.bud <- Map(getAgesExtantSpecies, tax.25.bud, trees.25.missing,
                   Tol = 1e-6)

## 50% missing species
tax.50.bud <- lapply(trees.50.missing, sim.taxonomy, beta = 0)

## get ages extant species
ages.50.bud <- Map(getAgesExtantSpecies, tax.50.bud, trees.50.missing,
                   Tol = 1e-6)

###########nonrandom incomple sampling of extant species#######################

##the nonrandom probability increases with age

##get ages extant species with a nonrandom incomple sampling (0.25)
ages.bud.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bud, trees.25.missing,
                          SamplingFrac = 0.25,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bud.incomp.25 <- do.call("rbind", ages.bud.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bud.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species name
ages.bud.incomp.25$species <- paste0(ages.bud.incomp.25$label,".",
                                     ages.bud.incomp.25$tree)

##sampling fraction
ages.bud.incomp.25$fraction <- rep(0.25, nrow(ages.bud.incomp.25))


##get ages extant species with a nonrandom incomple sampling (0.50)
ages.bud.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bud, trees.50.missing,
                          SamplingFrac = 0.50,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bud.incomp.50 <- do.call("rbind", ages.bud.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bud.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species name
ages.bud.incomp.50$species <- paste0(ages.bud.incomp.50$label,".",
                                     ages.bud.incomp.50$tree)

##sampling fraction
ages.bud.incomp.50$fraction <- rep(0.50, nrow(ages.bud.incomp.50))


##0% missing species

#collapsing in a DF
ages.0.bud <- do.call("rbind", ages.0.bud) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.0.bud$root.age <- rep(root.0, each = 100)

##adding rTrue.age
ages.0.bud <- ages.0.bud %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.0.bud$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species
ages.0.bud$species <- paste0(ages.0.bud$label,".",
                             ages.0.bud$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.0.bud$status <- discretize(ages.0.bud$rTrue.age, 
                                method = "frequency", 
                                breaks = 5,
                                labels = c("LC", "NT", "VU", "EN", "CR"))

ages.0.bud$status <- factor(ages.0.bud$status,
                            levels = c("LC", "NT", "VU", "EN", "CR"))


##25% missing species

ages.25.bud <- do.call("rbind", ages.25.bud) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.25.bud$root.age <- rep(root.25, each = 134)

##adding rTrue.age
ages.25.bud <- ages.25.bud %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.25.bud$tree <- rep(paste0("tree.", 1:1000), each = 134)

##species
ages.25.bud$species <- paste0(ages.25.bud$label,".",
                              ages.25.bud$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.25.bud$status <- discretize(ages.25.bud$rTrue.age, 
                                 method = "frequency", 
                                 breaks = 5,
                                 labels = c("LC", "NT", "VU", "EN", "CR"))

ages.25.bud$status <- factor(ages.25.bud$status,
                             levels = c("LC", "NT", "VU", "EN", "CR"))


##50% missing species
ages.50.bud <- do.call("rbind", ages.50.bud) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.50.bud$root.age <- rep(root.0, each = 200)

##adding rTrue.age
ages.50.bud <- ages.50.bud %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.50.bud$tree <- rep(paste0("tree.", 1:1000), each = 200)

##species
ages.50.bud$species <- paste0(ages.50.bud$label,".",
                              ages.50.bud$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.50.bud$status <- discretize(ages.50.bud$rTrue.age, 
                                 method = "frequency", 
                                 breaks = 5,
                                 labels = c("LC", "NT", "VU", "EN", "CR"))

ages.50.bud$status <- factor(ages.50.bud$status,
                             levels = c("LC", "NT", "VU", "EN", "CR"))


##selecting ages.bud.0.25 columns
ages.bud.incomp.25 <- ages.bud.incomp.25 %>% select(species,  
                                                    Incomp.age) %>% 
  rename(Incomp.age.25 = Incomp.age)

##merging ages.bud with ages.bud.0.25
ages.total.0.25.bud <- left_join(ages.bud.incomp.25, ages.25.bud, 
                             by = "species")

##selecting ages.bud.0.5 columns
ages.bud.incomp.50 <- ages.bud.incomp.50 %>% select(species,
                                                    Incomp.age) %>% 
  rename(Incomp.age.50 = Incomp.age)

##merging the three dataframes
ages.total.0.50.bud <- left_join(ages.bud.incomp.50, ages.50.bud,
                                 by = "species")

##Calculating the mean of each age for each tree for 0% missing sp
ages.mean.0.bud <-ages.0.bud %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.phy = mean(Estimated.age))
write_csv(ages.mean.0.bud,
          file = "results/data/processed/incomplete_sampling/ages.mean.0.bud.csv")
##Calculating the mean of each age for each tree for 25% missing sp                        
ages.mean.25.bud <- ages.total.0.25.bud %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.0.25 = mean(Incomp.age.25, na.rm = TRUE))

write_csv(ages.mean.25.bud,
          file = "results/data/processed/incomplete_sampling/ages.mean.25.bud.csv")
####calculating the mean of each age for each tree for 50% missing sp
ages.mean.50.bud <- ages.total.0.50.bud %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.0.50 = mean(Incomp.age.50, na.rm = TRUE)) 

write_csv(ages.mean.50.bud,
          file = "results/data/processed/incomplete_sampling/ages.mean.50.bud.csv")

#######ages ranking 0%
ages.rank.0.bud <- ages.mean.0.bud %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         phy.rank = dense_rank(mean.phy)) %>% 
  mutate(resp_phy = if_else(true.rank == phy.rank, 1, 0))


###comparing phylogenetic age
ages.comparison.0.bud <- ages.rank.0.bud %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

## ages rank 25%
ages.rank.25.bud <- ages.mean.25.bud %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         f25.rank = dense_rank(mean.0.25)) %>% 
  #f50.rank = dense_rank(mean.0.50)) %>% 
  mutate(resp_f25 = if_else(true.rank == f25.rank, 1, 0))
#resp_f50 = if_else(true.rank == f50.rank, 1, 0))

##comparing incomplete sampling 25%
ages.comparison.f25.bud <- ages.rank.25.bud %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25, na.rm = TRUE)) %>% 
  filter(sum_f25 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

## ages rank 50%
ages.rank.50.bud <- ages.mean.50.bud %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         f50.rank = dense_rank(mean.0.50)) %>% 
  mutate(resp_f50 = if_else(true.rank == f50.rank, 1, 0))

##comparing incomplete sampling 50%
ages.comparison.f50.bud <- ages.rank.50.bud %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)



##########plots#########################################

## themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (9),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (11),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (10)),
                      legend.background = element_rect(fill = "gray90",
                                                       size = 0.5, linetype = "dotted"),
                      legend.position = "bottom")




# Bifurcating -------------------------------------------------------------

###0% missing species

##correct estimations
tree.bif.0.correct <- ages.rank.0.bif %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy == 5) %>% 
  pull(tree)

tree.bif.0.correct <- sample(tree.bif.0.correct,
                             size = round(length(tree.bif.0.correct)*0.1))

##incorrect estimations 
tree.bif.0.incorrect <- ages.rank.0.bif %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy < 5) %>% 
  pull(tree)

tree.bif.0.incorrect <- sample(tree.bif.0.incorrect,
                             size = round(length(tree.bif.0.incorrect)*0.1))

###25% missing species

##correct estimations 
tree.bif.correct.f25 <- ages.rank.25.bif %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25)) %>% 
  filter(sum_f25 == 5) %>% 
  pull(tree)

tree.bif.correct.f25 <- sample(tree.bif.correct.f25,
                             size = round(length(tree.bif.correct.f25)*0.1))
##incorrect estimations (19%)
tree.bif.incorrect.f25 <- ages.rank.25.bif %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25, na.rm = TRUE)) %>% 
  filter(sum_f25 < 5) %>% 
  pull(tree)

tree.bif.incorrect.f25 <- sample(tree.bif.incorrect.f25,
                               size = round(length(tree.bif.incorrect.f25)*0.1))

###50% missing species

##correct estimation 
tree.bif.correct.f50 <- ages.rank.50.bif %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 == 5) %>%
  pull(tree)

tree.bif.correct.f50 <- sample(tree.bif.correct.f50,
                               size = round(length(tree.bif.correct.f50)*0.1))
##incorrect estimation 
tree.bif.incorrect.f50 <- ages.rank.50.bif %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 < 5) %>% 
  pull(tree)

tree.bif.incorrect.f50 <- sample(tree.bif.incorrect.f50,
                               size = round(length(tree.bif.incorrect.f50)*0.1))




###fake facets for bifurcating
ages.mean.0.bif$f0 <- "Fully sampled"
ages.mean.25.bif$f25 <- "25% missing species"
ages.mean.50.bif$f50 <- "50% missing species"



ages.mean.0.bif$status <- factor(ages.mean.0.bif$status, 
                                 levels = c("LC", "NT", "VU", "EN", "CR"),
                                 ordered= TRUE)
ages.mean.25.bif$status <- factor(ages.mean.25.bif$status,
                                  levels = c("LC", "NT", "VU", "EN", "CR"),
                                  ordered= TRUE)

ages.mean.50.bif$status <- factor(ages.mean.50.bif$status,
                                  levels = c("LC", "NT", "VU", "EN", "CR"),
                                  ordered= TRUE)

###0% missing species

##correct
ages.0.correct.bif <- ages.mean.0.bif %>%
                     filter(tree %in% tree.bif.0.correct)

##incorrect
ages.0.incorrect.bif <- ages.mean.0.bif %>% 
              filter(tree %in% tree.bif.0.incorrect)


f0.bif <- ggplot(ages.0.correct.bif,
                  aes(x = status, y = log(mean.phy +1), group = tree,
                           color = "#1b9e77"))+
  geom_point(size=3, shape = 21, fill="white", color = "#1b9e77")+
  geom_line(alpha = 0.3, color = "#1b9e77")+
  geom_point(data = ages.0.incorrect.bif,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.0.incorrect.bif,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f0.bif, aes(x = 2.0, y = 2.0,
                     label = paste0("Error rate:", " ",error,
                                    "%")),
            inherit.aes = FALSE,
            size = 4, 
            family = "serif")+
  theme_bw()+
  ylab("log(Phylogenetic age)")+
  xlab(NULL)+
  #ggtitle("Bifurcating speciation")+
  facet_wrap(~f0)+
  mynamestheme+
  theme(legend.position = "none")

##Incomplete sampling 25%

##correct
ages.25.correct.bif <- ages.mean.25.bif %>% filter(tree %in% tree.bif.correct.f25)

##incorrect
ages.25.incorrect.bif <- ages.mean.25.bif%>% 
                           filter(tree %in% tree.bif.incorrect.f25)


f25.bif <- ggplot(ages.25.correct.bif,
                  aes(x = status, y = log(mean.0.25 +1), group = tree,
             color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.3, color = "#d95f02")+
  geom_point(data = ages.25.incorrect.bif,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.25.incorrect.bif,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f25.bif, aes(x = 2.0, y = 2.0,
                              label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 4, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f25)+
  mynamestheme+
  theme(legend.position = "none")


##Incomplete sampling 50%

###correct
ages.50.correct.bif <- ages.mean.50.bif %>%
                          filter(tree %in% tree.bif.correct.f50)

##incorrect
ages.50.incorrect.bif <- ages.mean.50.bif %>%
                          filter(tree %in% tree.bif.incorrect.f50)


f50.bif <- ggplot(ages.50.correct.bif,
                  aes(x = status, y = log(mean.0.50 +1), group = tree,
             color = "red"))+
  geom_point(size=3, shape = 21, fill="white", color = "red")+
  geom_line(alpha = 0.4, color = "red")+
  geom_point(data = ages.50.incorrect.bif,
             size=3, shape = 21, fill="white",
             alpha = 0.3, color = "#969696")+
  geom_line(data = ages.50.incorrect.bif,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f50.bif, aes(x = 2.0, y = 2.35,
                               label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 4, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f50)+
  mynamestheme+
  theme(legend.position = "none")


# Budding -----------------------------------------------------------------

###0% missing species

##correct estimations 
tree.bud.0.correct <- ages.rank.0.bud %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy == 5) %>% 
  pull(tree)

tree.bud.0.correct <- sample(tree.bud.0.correct,
                             size = round(length(tree.bud.0.correct)*0.1))

##incorrect estimations 
tree.bud.0.incorrect <- ages.rank.0.bud %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy < 5) %>% 
  pull(tree)

tree.bud.0.incorrect <- sample(tree.bud.0.incorrect,
                             size = round(length(tree.bud.0.incorrect)*0.1))

###25% missing species

##correct estimations 
tree.bud.correct.f25 <- ages.rank.25.bud %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25)) %>% 
  filter(sum_f25 == 5) %>% 
  pull(tree)

tree.bud.correct.f25 <- sample(tree.bud.correct.f25,
                             size = round(length(tree.bud.correct.f25)*0.1))

##incorrect estimations 
tree.bud.incorrect.f25 <- ages.rank.25.bud %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25, na.rm = TRUE)) %>% 
  filter(sum_f25 < 5) %>% 
  pull(tree)

tree.bud.incorrect.f25 <- sample(tree.bud.incorrect.f25,
                               size = round(length(tree.bud.incorrect.f25)*0.1))

###50% missing species

##correct estimation 
tree.bud.correct.f50 <- ages.rank.50.bud %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 == 5) %>% 
  pull(tree)

tree.bud.correct.f50 <- sample(tree.bud.correct.f50,
                               size = round(length(tree.bud.correct.f50)*0.1))
##incorrect estimation 
tree.bud.incorrect.f50 <- ages.rank.50.bud %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50, na.rm = TRUE)) %>% 
  filter(sum_f50 < 5) %>% 
  pull(tree)

tree.bud.incorrect.f50 <- sample(tree.bud.incorrect.f50,
                               size = round(length(tree.bud.incorrect.f50)*0.1))

###fake facets for bifurcating
ages.mean.0.bud$f0 <- "0% missing species"
ages.mean.25.bud$f25 <- "25% missing species"
ages.mean.50.bud$f50 <- "50% missing species"


###0% missing species
age.bud.mean.0.correct <- ages.mean.0.bud %>%
                                 filter(tree %in% tree.bud.0.correct) 


age.bud.mean.0.incorrect <- ages.mean.0.bud %>% 
                               filter(tree %in% tree.bud.0.incorrect)

f0.bud <- ggplot(age.bud.mean.0.correct,
                 aes(x = status, y = log(mean.phy +1), group = tree,
                     color = "#1b9e77"))+
  geom_point(size=3, shape = 21, fill="white", color = "#1b9e77")+
  geom_line(alpha = 0.4, color = "#1b9e77")+
  geom_point(data = age.bud.mean.0.incorrect,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = age.bud.mean.0.incorrect,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.0.bud, aes(x = 2.0, y = 2.0,
                                             label = paste0("Error rate:", " ",error,
                                                               "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle("Budding speciation")+
  facet_wrap(~f0)+
  mynamestheme+
  theme(legend.position = "none")

#####25% missing species

##correct trees
age.bud.mean.25.correct <-ages.mean.25.bud %>%
                    filter(tree %in% tree.bud.correct.f25) 


age.bud.mean.25.incorrect <- ages.mean.25.bud %>% 
                            filter(tree %in% tree.bud.incorrect.f25)

f25.bud <- ggplot(age.bud.mean.25.correct,
                  aes(x = status, y = log(mean.0.25 +1), group = tree,
             color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.4, color = "#d95f02")+
  geom_point(data = age.bud.mean.25.incorrect,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = age.bud.mean.25.incorrect,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f25.bud, aes(x = 2.0, y = 1.7,
                                       label = paste0("Error rate:", " ",error,
                                                               "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f25)+
  mynamestheme+
  theme(legend.position = "none")


###50% missing species

##correct trees

age.bud.mean.50.correct <-ages.mean.50.bud %>%
                      filter(tree %in% tree.bud.correct.f50) 


age.bud.mean.50.incorrect <- ages.mean.50.bud %>% 
                              filter(tree %in% tree.bud.incorrect.f50)

##Incomplete sampling 50%
f50.bud <- ggplot(age.bud.mean.50.correct,
                  aes(x = status, y = log(mean.0.50 +1), group = tree,
             color = "red"))+
            geom_point(size=3, shape = 21, fill="white", color = "red")+
            geom_line(alpha = 0.9, color = "red")+
            geom_point(data = age.bud.mean.50.incorrect,
                       size=3, shape = 21, fill="white",
                       alpha = 0.4, color = "#969696")+
            geom_line(data = age.bud.mean.50.incorrect,
                      alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.f50.bud, aes(x = 1.7, y = 1.9,
                              label = paste0("Error rate:", " ",error,
                                                               "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f50)+
  mynamestheme+
  theme(legend.position = "none")

###unifying plots

###png object
png("text/figures/nonrandom_CS.png", 
    width = 14, height = 17,
    units = "cm", 
    pointsize = 8, res = 300)

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
x.bif = grid.arrange(f0.bif, f25.bif, f50.bif, nrow = 3)
x.bud = grid.arrange(f0.bud, f25.bud, f50.bud, nrow = 3)

grid.arrange(x.bif, x.bud, ncol = 2, bottom = bottom)
dev.off()



##Only bifurcating speciation
###png object
png("text/figures/bif_nonrandom_CS.png", 
    width = 25, height = 10,
    units = "cm", 
    pointsize = 8, res = 300)

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##top
top <- text_grob("Bifurcating speciation",
                 family = "serif", size = 13,  face = "bold")

grid.arrange(f0.bif, f25.bif, f50.bif, ncol = 3, top = top,
             bottom = bottom)
dev.off()
