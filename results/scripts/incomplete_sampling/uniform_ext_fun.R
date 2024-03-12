
# Testing the function wit different extinction and missing species -------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


# no extinction scenario -------------------------------------------------

## 25% missing species
trees.25.missing <- sim.bd.taxa(134, numbsim = 1000,
                                lambda = 0.3, mu = 0)

## 50% missing species
trees.50.missing <- sim.bd.taxa(200, numbsim = 1000,
                                lambda = 0.3, mu = 0)

## 25% missing species
trees.25.extant <- lapply(trees.25.missing, prune.fossil.tips)

##root ages
root.25 <- sapply(trees.25.missing, tree.max)

### 50% missing species
trees.50.extant <- lapply(trees.50.missing, prune.fossil.tips)

## root ages
root.50 <- sapply(trees.50.missing, tree.max)

# bifurcating  ------------------------------------------------------------

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

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bif, trees.25.missing,
                          SamplingFrac = 0.25)


ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.25$root.age <- rep(root.25, each = 100)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.25 <- ages.bif.incomp.25 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.25$extinction <- rep("no_ext", nrow(ages.bif.incomp.25))

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


###speciation mode
ages.bif.incomp.25$speciation <- rep("bif", nrow(ages.bif.incomp.25))

#####sampling fraction = 0.50##########

##get ages extant species with an uniform incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bif, trees.50.missing,
                          SamplingFrac = 0.50)


ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.50$root.age <- rep(root.50, each = 100)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.50 <- ages.bif.incomp.50 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.50$extinction <- rep("no_ext", nrow(ages.bif.incomp.50))

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))

###speciation mode
ages.bif.incomp.50$speciation <- rep("bif", nrow(ages.bif.incomp.50))

##binding dataframes
ages_incomplete_no <- rbind(ages.bif.incomp.25, ages.bif.incomp.50)

##saving
write_csv(ages_incomplete_no, file = 
            "results/data/processed/incomplete_sampling/ages.incomplete.no.ex.csv")

##read
ages_incomplete_no <- read_csv("results/data/processed/incomplete_sampling/ages.incomplete.no.ex.csv")

# Probabilistic function --------------------------------------------------

###for 25% incomplete sampling
ages.incomplete_no.f25 <- ages_incomplete_no %>% filter(fraction == "0.25")

##trying new function provided by the reviewer

###mean age
f25_mean_new <- vector(mode = "numeric", length = nrow(ages.incomplete_no.f25))


for(i in 1:nrow(ages.incomplete_no.f25)){
  f25_mean_new[[i]] <- meanAgeFunction(lambda = 0.3,
                                  mu = 0,
                                  rho = 0.75,
                                  v = ages.incomplete_no.f25$Incomp.age[i])
}

##median age
f25_median_new <- vector(mode = "numeric", length = nrow(ages.incomplete_no.f25))


for(i in 1:nrow(ages.incomplete_no.f25)){
  f25_median_new[[i]] <- medianAgeFunction(lambda = 0.3,
                                       mu = 0,
                                       rho = 0.75,
                                       v = ages.incomplete_no.f25$Incomp.age[i])
}

##binding columns
ages.f25_no <- cbind(ages.incomplete_no.f25,
                     f25_mean_new,
                     f25_median_new) %>% 
  rename(mean.age = f25_mean_new,
         median.age = f25_median_new)



# 50% missing species -----------------------------------------------------

ages.incomplete_no.f50 <- ages_incomplete_no %>% filter(fraction == "0.5")

## new function provided by the reviewer

###mean age
f50_mean_new <- vector(mode = "numeric", length = nrow(ages.incomplete_no.f50))


for(i in 1:nrow(ages.incomplete_no.f50)){
  f50_mean_new[[i]] <- meanAgeFunction(lambda = 0.3,
                                       mu = 0,
                                       rho = 0.50,
                                       v = ages.incomplete_no.f50$Incomp.age[i])
}

##median age
f50_median_new <- vector(mode = "numeric", length = nrow(ages.incomplete_no.f50))


for(i in 1:nrow(ages.incomplete_no.f50)){
  f50_median_new[[i]] <- medianAgeFunction(lambda = 0.3,
                                           mu = 0,
                                           rho = 0.50,
                                           v = ages.incomplete_no.f50$Incomp.age[i])
}

##binding columns
ages.f50_no <- cbind(ages.incomplete_no.f50,
                     f50_mean_new,
                     f50_median_new) %>% 
  rename(mean.age = f50_mean_new,
         median.age = f50_median_new)


ages.random_no_ext <- rbind(ages.f25_no, ages.f50_no)

###saving
write_csv(ages.random_no_ext,
      file = "results/data/processed/incomplete_sampling/ages.random_no_ext.csv")

#ages.random_no_ext <- read_csv("results/data/processed/incomplete_sampling/ages.random_no_ext.csv")

###calculating MAPE

ages.mape.no_ext <- ages.random_no_ext %>% group_by(fraction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Incomp.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)


ages.mape.no_ext$extinction <- "no_ext"



##binding mape dataframes
ages.mape.sampling_no <- ages.mape.no_ext %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.sampling_no$estimate <- factor(ages.mape.sampling_no$estimate,
                                          levels = c("mape.incomplete",
                                                     "mape.mean",
                                                     "mape.median"),
                                          ordered = TRUE)

# Intermediate extinction -------------------------------------------------

## 25% missing species
trees.25.missing <- sim.bd.taxa(134, numbsim = 1000,
                                lambda = 0.3, mu = 0.15)

## 50% missing species
trees.50.missing <- sim.bd.taxa(200, numbsim = 1000,
                                lambda = 0.3, mu = 0.15)

## 25% missing species
trees.25.extant <- lapply(trees.25.missing, prune.fossil.tips)

##root ages
root.25 <- sapply(trees.25.missing, tree.max)

### 50% missing species
trees.50.extant <- lapply(trees.50.missing, prune.fossil.tips)

## root ages
root.50 <- sapply(trees.50.missing, tree.max)

# bifurcating  ------------------------------------------------------------

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

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bif, trees.25.missing,
                          SamplingFrac = 0.25)


ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.25$root.age <- rep(root.25, each = 100)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.25 <- ages.bif.incomp.25 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.25$extinction <- rep("intermediate", nrow(ages.bif.incomp.25))

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


###speciation mode
ages.bif.incomp.25$speciation <- rep("bif", nrow(ages.bif.incomp.25))

#####sampling fraction = 0.50##########

##get ages extant species with an uniform incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bif, trees.50.missing,
                          SamplingFrac = 0.50)


ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.50$root.age <- rep(root.50, each = 100)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.50 <- ages.bif.incomp.50 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.50$extinction <- rep("intermediate", nrow(ages.bif.incomp.50))

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))

###speciation mode
ages.bif.incomp.50$speciation <- rep("bif", nrow(ages.bif.incomp.50))

##binding dataframes
ages_incomplete_int <- rbind(ages.bif.incomp.25, ages.bif.incomp.50)


### new funciton provided by the reviewer

#mean

f25.mean.age.new <- vector(mode = "numeric",
                           length = nrow(ages.bif.incomp.25))

for(i in 1:nrow(ages.bif.incomp.25)){
  f25.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0.15,
                                           rho = 0.75,
                                           v = ages.bif.incomp.25$Incomp.age[i])
}

##median
f25.median.age.new <- vector(mode = "numeric",
                           length = nrow(ages.bif.incomp.25))

for(i in 1:nrow(ages.bif.incomp.25)){
  f25.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                           mu = 0.15,
                                           rho = 0.75,
                                           v = ages.bif.incomp.25$Incomp.age[i])
}

##binding columns
ages.f25 <- cbind(ages.bif.incomp.25, f25.mean.age.new,
                  f25.median.age.new) %>% 
  rename(mean.age = f25.mean.age.new,
         median.age = f25.median.age.new)



## 50% sampling


##new function provided by the reviewer
#mean

f50.mean.age.new <- vector(mode = "numeric",
                           length = nrow(ages.bif.incomp.50))

for(i in 1:nrow(ages.bif.incomp.50)){
  f50.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0.15,
                                           rho = 0.50,
                                           v = ages.bif.incomp.50$Incomp.age[i])
}

##median
f50.median.age.new <- vector(mode = "numeric",
                             length = nrow(ages.bif.incomp.50))

for(i in 1:nrow(ages.bif.incomp.50)){
  f50.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                               mu = 0.15,
                                               rho = 0.50,
                                               v = ages.bif.incomp.50$Incomp.age[i])
}


##binding columns
ages.f50 <- cbind(ages.bif.incomp.50,
                  f50.mean.age.new, f50.median.age.new) %>% 
  rename(mean.age = f50.mean.age.new,
         median.age = f50.median.age.new)


##binding dataframes
ages.random_int_ext <- rbind(ages.f25, ages.f50) 

ages.random_int_ext$extinction <- "intermediate"

##saving
write_csv(ages.random_int_ext,
         file = "results/data/processed/incomplete_sampling/ages.random_int_ext.csv")

#ages.random_int_ext <- read_csv(file = "results/data/processed/incomplete_sampling/ages.random_int_ext.csv")

##MAPE
ages.mape.int <- ages.random_int_ext %>% group_by(fraction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Incomp.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)

ages.mape.int$extinction <- "intermediate"

##pivot longer
ages.mape.sampling_int <- ages.mape.int %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.sampling_int$estimate <- factor(ages.mape.sampling_int$estimate,
                                           levels = c("mape.incomplete",
                                                      "mape.mean",
                                                      "mape.median"),
                                           ordered = TRUE)



# High extinction ---------------------------------------------------------

## 25% missing species
trees.25.missing <- sim.bd.taxa(134, numbsim = 1000,
                                lambda = 0.3, mu = 0.25)

## 50% missing species
trees.50.missing <- sim.bd.taxa(200, numbsim = 1000,
                                lambda = 0.3, mu = 0.25)

## 25% missing species
trees.25.extant <- lapply(trees.25.missing, prune.fossil.tips)

##root ages
root.25 <- sapply(trees.25.missing, tree.max)

### 50% missing species
trees.50.extant <- lapply(trees.50.missing, prune.fossil.tips)

## root ages
root.50 <- sapply(trees.50.missing, tree.max)

# bifurcating  ------------------------------------------------------------

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

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bif, trees.25.missing,
                          SamplingFrac = 0.25)


ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.25$root.age <- rep(root.25, each = 100)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.25 <- ages.bif.incomp.25 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.25$extinction <- rep("high", nrow(ages.bif.incomp.25))

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


###speciation mode
ages.bif.incomp.25$speciation <- rep("bif", nrow(ages.bif.incomp.25))

#####sampling fraction = 0.50##########

##get ages extant species with an uniform incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bif, trees.50.missing,
                          SamplingFrac = 0.50)


ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.50$root.age <- rep(root.50, each = 100)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.50 <- ages.bif.incomp.50 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.50$extinction <- rep("high", nrow(ages.bif.incomp.50))

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))

###speciation mode
ages.bif.incomp.50$speciation <- rep("bif", nrow(ages.bif.incomp.50))

##binding dataframes
ages_incomplete_high <- rbind(ages.bif.incomp.25, ages.bif.incomp.50)

##saving
write_csv(ages_incomplete_high, file = 
            "results/data/processed/incomplete_sampling/ages.incomplete.high.ex.csv")

#ages.incomplete_high <- read_csv("results/data/processed/incomplete_sampling/ages.incomplete.high.ex.csv")

# Probabilistic function --------------------------------------------------


##new function provided by the reviewer
#mean

f25.mean.age.new <- vector(mode = "numeric",
                           length = nrow(ages.bif.incomp.25))

for(i in 1:nrow(ages.bif.incomp.25)){
  f25.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0.25,
                                           rho = 0.75,
                                           v = ages.bif.incomp.25$Incomp.age[i])
}

##median
f25.median.age.new <- vector(mode = "numeric",
                             length = nrow(ages.bif.incomp.25))

for(i in 1:nrow(ages.bif.incomp.25)){
  f25.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                               mu = 0.25,
                                               rho = 0.75,
                                               v = ages.bif.incomp.25$Incomp.age[i])
}


##binding columns
ages.f25_high <- cbind(ages.bif.incomp.25,
                     f25.mean.age.new, 
                     f25.median.age.new) %>% 
                  rename(mean.age = f25.mean.age.new,
                         median.age = f25.median.age.new)


# 50% missing species -----------------------------------------------------

##new function provided by the reviewer

#mean

f50.mean.age.new <- vector(mode = "numeric",
                           length = nrow(ages.bif.incomp.50))

for(i in 1:nrow(ages.bif.incomp.50)){
  f50.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0.25,
                                           rho = 0.50,
                                           v = ages.bif.incomp.50$Incomp.age[i])
}

##median
f50.median.age.new <- vector(mode = "numeric",
                             length = nrow(ages.bif.incomp.50))

for(i in 1:nrow(ages.bif.incomp.50)){
  f50.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                               mu = 0.25,
                                               rho = 0.50,
                                               v = ages.bif.incomp.50$Incomp.age[i])
}



##binding columns
ages.f50_high <- cbind(ages.bif.incomp.50,
                       f50.mean.age.new,
                     f50.median.age.new) %>% 
                    rename(mean.age = f50.mean.age.new,
                           median.age = f50.median.age.new)


ages.random_high_ext <- rbind(ages.f25_high, ages.f50_high)

write_csv(ages.random_high_ext,
          file = "results/data/processed/incomplete_sampling/ages.random_high_ext.csv")

#ages.random_high_ext <- read_csv(file = "results/data/processed/incomplete_sampling/ages.random_high_ext.csv")

###calculating MAPE

ages.mape.high <- ages.random_high_ext %>% group_by(fraction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Incomp.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)

ages.mape.high$extinction <- "high"


##binding mape dataframes
ages.mape.sampling_high <- ages.mape.high %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.sampling_high$estimate <- factor(ages.mape.sampling_high$estimate,
                                         levels = c("mape.incomplete",
                                                    "mape.mean",
                                                    "mape.median"),
                                         ordered = TRUE)

ages.mape.sampling_high$fraction <- as.factor(ages.mape.sampling_high$fraction)





# Fully sampled evaluation ------------------------------------------------

#using the generated dataframe because we estimated the species ages without considering
#the sampling shortfall

# no extinction -----------------------------------------------------------

ages.full_no <- ages.random_no_ext %>% filter(fraction == 0.25) %>% 
                select(-c(mean.age, median.age))


# probability function ----------------------------------------------------


#using the new function provided by the reviewer
#mean

full.no.mean.age.new <- vector(mode = "numeric",
                           length = nrow(ages.full_no))

for(i in 1:nrow(ages.full_no)){
  full.no.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                           mu = 0,
                                           rho = 1,
                                           v = ages.full_no$Estimated.age[i])
}

##median
full.no.median.age.new <- vector(mode = "numeric",
                             length = nrow(ages.full_no))

for(i in 1:nrow(ages.full_no)){
  full.no.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                               mu = 0,
                                               rho = 1,
                                               v = ages.full_no$Estimated.age[i])
}


##dataframe
ages.full_no <- cbind(ages.full_no,
                           full.no.mean.age.new, full.no.median.age.new) %>% 
                     rename(mean.age = full.no.mean.age.new,
                            median.age = full.no.median.age.new)


# intermediate extinction -------------------------------------------------

ages.full_int <- ages.random_int_ext %>% filter(fraction == 0.25) %>% 
               select(-c(mean.age, median.age))


# probability function ----------------------------------------------------

######new function 
#mean

full.int.mean.age.new <- vector(mode = "numeric",
                               length = nrow(ages.full_int))

for(i in 1:nrow(ages.full_no)){
  full.int.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                               mu = 0.15,
                                               rho = 1,
                                               v = ages.full_int$Estimated.age[i])
}

##median
full.int.median.age.new <- vector(mode = "numeric",
                                 length = nrow(ages.full_int))

for(i in 1:nrow(ages.full_int)){
  full.int.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                                   mu = 0.15,
                                                   rho = 1,
                                                   v = ages.full_int$Estimated.age[i])
}

##dataframe
ages.full_int_prob <- cbind(ages.full_int,
                            full.int.mean.age.new, full.int.median.age.new) %>% 
                  rename(mean.age = full.int.mean.age.new,
                         median.age = full.int.median.age.new)


# High extinction ---------------------------------------------------------

ages.full_high  <- ages.random_high_ext %>% filter(fraction == 0.25) %>% 
  select(-c(mean.age, median.age))


# probability function ----------------------------------------------------


##new function
#mean

full.high.mean.age.new <- vector(mode = "numeric",
                                length = nrow(ages.full_high))

for(i in 1:nrow(ages.full_high)){
  full.high.mean.age.new[[i]] <- meanAgeFunction(lam = 0.3,
                                                mu = 0.25,
                                                rho = 1,
                                                v = ages.full_high$Estimated.age[i])
}

##median
full.high.median.age.new <- vector(mode = "numeric",
                                  length = nrow(ages.full_high))

for(i in 1:nrow(ages.full_high)){
  full.high.median.age.new[[i]] <- medianAgeFunction(lam = 0.3,
                                                    mu = 0.25,
                                                    rho = 1,
                                                    v = ages.full_high$Estimated.age[i])
}

##dataframe
ages.full_high_prob <- cbind(ages.full_high,
                            full.high.mean.age.new, full.high.median.age.new) %>% 
  rename(mean.age = full.high.mean.age.new,
         median.age = full.high.median.age.new)

##binding dataframes
ages.full_total_prob <- rbind(ages.full_no,
                              ages.full_int_prob,
                              ages.full_high_prob)

##saving
write_csv(ages.full_total_prob,
      file = "results/data/processed/incomplete_sampling/ages.full_total_prob.csv")

#ages.full_total_prob <- 
  #read_csv(file = "results/data/processed/incomplete_sampling/ages.full_total_prob.csv")

####MAPE

ages.mape.full <- ages.full_total_prob %>% group_by(extinction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Estimated.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)
            
ages.mape.full$fraction <- "full"


##binding mape dataframes
ages.mape.full.sampling <- ages.mape.full %>% 
              pivot_longer(cols = starts_with("mape."),
                           values_to = "mape", names_to = "estimate")

ages.mape.full.sampling$estimate <- factor(ages.mape.full.sampling$estimate,
                                           levels = c("mape.incomplete",
                                                      "mape.mean",
                                                      "mape.median"),
                                           ordered = TRUE)

ages.mape.full.sampling$extinction <- factor(ages.mape.full.sampling$extinction,
                                             levels = c("no_ext",
                                                        "intermediate",
                                                        "high"))
ages.mape.full.sampling$fraction <- "full"




###binding all age mape
ages.mape.total <- rbind(ages.mape.full, ages.mape.no_ext,
                         ages.mape.int, ages.mape.high)

ages.mape.total$extinction <- factor(ages.mape.total$extinction,
                                     levels = c("no_ext",
                                                "intermediate",
                                                "high"),
                                     ordered = TRUE)

ages.mape.total$fraction <- factor(ages.mape.total$fraction,
                                     levels = c("full",
                                                "0.25",
                                                "0.5"),
                                     ordered = TRUE)

write_csv(ages.mape.total, 
          file = "results/data/processed/incomplete_sampling/ages.mape.total.csv")

ages.mape.total <- read_csv("results/data/processed/incomplete_sampling/ages.mape.total.csv")

######binding all ages.mape.sampling df
ages.mape.sampling.tot <- rbind(ages.mape.sampling_no, ages.mape.sampling_int,
                                ages.mape.sampling_high, ages.mape.full.sampling)

##arranging factor
ages.mape.sampling.tot$extinction <- factor(ages.mape.sampling.tot$extinction,
                                            levels = c("no_ext",
                                                       "intermediate",
                                                       "high"))
ages.mape.sampling.tot$fraction <- factor(ages.mape.sampling.tot$fraction,
                                          levels = c("full",
                                                     "0.25",
                                                     "0.5"))
##saving
write_csv(ages.mape.sampling.tot, 
  file = "results/data/processed/incomplete_sampling/ages.mape.sampling.csv")

##read
ages.mape.sampling <- read_csv("results/data/processed/incomplete_sampling/ages.mape.sampling.csv")


# plots -------------------------------------------------------------------
## themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (11),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (9),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (9)),
                      legend.background = element_rect(fill="grey",
                                                       size=.5, linetype="dotted"),
                      legend.position = "bottom")

##########MAPE of phylogenetic age

png("text/figures/MAPE.sampling.extinction.png", 
    width = 17, height = 12, units = "cm", 
    pointsize = 8, res = 300)


ggplot(ages.mape.total, aes(x = fraction, y = mape.incomplete,
                            fill = fraction))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~extinction, labeller = as_labeller(c("no_ext" = "No extinction",
                                                   "intermediate" = "Intermediate extinction",
                                                   "high" = "High extinction")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    name="Missing species",
                    breaks=c("full", "0.25", "0.5"),
                    labels=c("0%", "25%",
                             "50%"))+
  ylim(0,1050)+
  ylab("MAPE")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank(),
        legend.background = element_rect(fill=NA,
                                         size=.5, linetype="dotted"),
        legend.position = c(0.1, 0.8))

dev.off()


##Fully sampled

full.plot <- ages.mape.sampling.tot %>% filter(fraction == "full") %>% 
         ggplot(aes(y = mape, x = estimate,
                                          fill = estimate))+
          geom_boxplot(outlier.shape = NA)+
        #geom_jitter(size = 0.1, alpha = 0.3)+
        facet_grid(fraction~extinction,
                   labeller = as_labeller(c("no_ext" = "No extinction",
                                            "intermediate" = "Intermediate extinction",
                                            "high" = "High extinction",
                                            "full" = "Fully sampled")))+
        scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                          name="Estimation",
                          breaks=c("mape.incomplete", "mape.mean",
                                   "mape.median"),
                          labels=c("Phylogenetic", "Mean",
                                   "Median"))+
        ylim(0, 130)+
        ylab(NULL)+
        xlab(NULL)+
        theme_bw()+
        mynamestheme+
        theme(axis.text.x = element_blank(),
              legend.position = "none")


#### 25% missing species

f25.plot <- ages.mape.sampling.tot %>% filter(fraction == "0.25") %>% 
  ggplot(aes(y = mape, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.25" = "25% missing species")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    name="Estimation",
                    breaks=c("mape.incomplete",
                             "mape.mean", "mape.median"),
                    labels=c("Phylogenetic", "Mean",
                             "Median"))+
  ylim(0, 350)+
  ylab(NULL)+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank(),
        legend.position = "none")

#### 50% missing species

f50.plot <- ages.mape.sampling.tot %>% filter(fraction == "0.5") %>% 
  ggplot(aes(y = mape, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.5" = "50% missing species")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    name="Estimation",
                    breaks=c("mape.incomplete", 
                             "mape.mean", "mape.median"),
                    labels=c("Phylogenetic", "Mean",
                             "Median"))+
  ylim(0, 1050)+
  ylab(NULL)+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())



##plot grid
png("text/figures/MAPE.sampling.new.png", 
    width = 17, height = 20, units = "cm", 
    pointsize = 8, res = 300)


##left
left <- text_grob("MAPE",
                    family = "serif", size = 13,  face = "bold",
                  rot = 90)

##grid for plotting both figures
grid.arrange(full.plot, f25.plot, f50.plot,
             nrow = 3, left = left)

dev.off()


###DELTA FIGURE

ages.delta <- ages.mape.total %>% group_by(fraction, extinction) %>%  
  summarise(delta.mean = mape.mean - mape.incomplete,
            delta.median =  mape.median - mape.incomplete) %>% 
      pivot_longer(cols = starts_with("delta."),
               values_to = "delta", names_to = "estimate")

#arranging factor
ages.delta$estimate <- factor(ages.delta$estimate,
                              levels = c("delta.mean",
                                         "delta.median"),
                              ordered = TRUE)

##plots

##Fully sampled
full.delta <- ages.delta %>% filter(fraction == "full") %>% 
      ggplot(aes(y = delta, x = estimate,
                 fill = estimate))+
        geom_boxplot(outlier.shape = NA)+
        facet_grid(fraction~extinction,
                   labeller = as_labeller(c("no_ext" = "No extinction",
                                            "intermediate" = "Intermediate extinction",
                                            "high" = "High extinction",
                                            "full" = "Fully sampled"
                                           )),
                   scale = "free_y")+
        scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                          name="Estimation",
                          breaks=c(
                                   "delta.mean", "delta.median"),
                          labels=c("Mean",
                                   "Median"))+
        geom_hline(yintercept = 0, linetype = "dashed",
                   size = 1, colour = "red")+
        #ylab(expression(bold(Delta* "MAPE")))+
        ylim(-100, 10)+
        xlab(NULL)+
        ylab(NULL)+
        theme_bw()+
        mynamestheme+
        theme(axis.text.x = element_blank(),
              legend.position = "none")


##25% missing species
f25.delta <- ages.delta %>% filter(fraction == "0.25",
                                   estimate %in% c("delta.mean",
                                                   "delta.median")) %>% 
  ggplot(aes(y = delta, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.25" = "25% missing species"
             )))+
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                    name="Estimation",
                    breaks=c("delta.mean", "delta.median"),
                    labels=c(
                             "Mean",
                             "Median"))+
  geom_hline(yintercept = 0, linetype = "dashed",
             size = 1, colour = "red")+
  #ylab(expression(bold(Delta* "MAPE")))+
  ylim(-350, 20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank(),
        legend.position = "none")

##50% missing species
f50.delta <- ages.delta %>% filter(fraction == "0.5",
                                   estimate %in% c("delta.mean",
                                                   "delta.median")) %>% 
  ggplot(aes(y = delta, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.5" = "50% missing species"
             )))+
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                    name="Estimation",
                    breaks=c(
                             "delta.mean", "delta.median"),
                    labels=c("Mean",
                             "Median"))+
  geom_hline(yintercept = 0, linetype = "dashed",
             size = 1, colour = "red")+
  #ylab(expression(bold(Delta* "MAPE")))+
  ylim(-1000, 20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())

##plot grid

png("text/figures/delta.MAPE.sampling.new.png", 
    width = 17, height = 20, units = "cm", 
    pointsize = 8, res = 300)


##left
left <- text_grob(expression(bold(Delta* "MAPE")),
                  family = "serif", size = 13,  face = "bold",
                  rot = 90)

##grid for plotting both figures
grid.arrange(full.delta, f25.delta, f50.delta,
             nrow = 3, left = left)

dev.off()
