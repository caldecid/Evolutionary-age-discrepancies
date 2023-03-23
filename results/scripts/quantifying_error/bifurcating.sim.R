###############################################################################
###########    Simulation of bifurcating speciation to compare ################
##########            phyletic age and true age                ################
###############################################################################


##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


##setting seed for replication
set.seed(18)

# Trees with lambda = 0.1 -------------------------------------------------

###speciation rate
lam <- 0.1

## turnover rate
turnover <- seq(from = 0, to = 0.99, length.out = 100)

###extinction rate
mu.vector <- lam * turnover

##diversification rates vector
r.vector <- (lam - mu.vector)


##Number of extant species in the simulated tree
Ntip <- 100

##simulating 100 trees with lam = 0.1
Tree.lam.01 <- vector("list", 100)

##for loop due to different extinction rates
for(i in seq_along(mu.vector)){
  Tree.lam.01[i] <- sim.bd.taxa(n = Ntip, lambda = lam, mu = mu.vector[i],
                                complete = TRUE, numbsim = 1)
  print(i)
}

###save Tree with the extinction patterns
save(Tree.lam.01, file = file.path(getwd(), pro, q_error, "Tree.lam.01.RData"))


###simulating species taxonomy; only bifurcating speciation

Taxa.lam.01 <- lapply(Tree.lam.01, sim.taxonomy, beta = 1, lambda.a = 0)

##Calculating estimated and true ages
TaxaExtant.01 <- Map(getAgesExtantSpecies, Taxa.lam.01, Tree.lam.01)

clado.df.01 <- do.call("rbind", TaxaExtant.01)

##diversification
clado.df.01$div <- rep(r.vector, each = Ntip)

##turnover rates
clado.df.01$turnover <- rep(turnover, each = Ntip)

##speciation rates
clado.df.01$lambda <- rep(lam, nrow(clado.df.01))

##number of species
clado.df.01$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.lam.01, tree.max)
clado.df.01$root.age <- rep(root.ages, each = Ntip)
##budding speciation = 0 ; bifurcating = 1
clado.df.01$mode <- rep(1, nrow(clado.df.01))
##anagenetic = 1; not anagenetic = 0
clado.df.01$anag <- rep(0, nrow(clado.df.01))
###extinction rates
clado.df.01$mu <- rep(mu.vector, each = Ntip)


##renaming columns
clado.df.bif.01 <- clado.df.01 %>% rename(True.age = Age,
                                Estimated.age = tip_length_reconstr,
                                N_sisters = n_sisters_reconstr) 


# trees with lam = 0.5 ----------------------------------------------------


###speciation rate
lam <- 0.5

## turnover rate
turnover <- seq(from = 0, to = 0.99, length.out = 100)

###extinction rate
mu.vector <- lam * turnover

##diversification rates vector
r.vector <- (lam - mu.vector)


##Number of extant species in the simulated tree
Ntip <- 100

##simulating 100 trees with lam = 0.1
Tree.lam.05 <- vector("list", 100)

##for loop due to different extinction rates
for(i in seq_along(mu.vector)){
  Tree.lam.05[i] <- sim.bd.taxa(n = Ntip, lambda = lam, mu = mu.vector[i],
                                complete = TRUE, numbsim = 1)
  print(i)
}

###save Tree with the extinction patterns
save(Tree.lam.05, file = file.path(getwd(), pro, q_error, "Tree.lam.05.RData"))

###simulating species taxonomy; only bifurcating speciation

Taxa.lam.0.5 <- lapply(Tree.lam.05, sim.taxonomy, beta = 1, lambda.a = 0)

##Calculating estimated and true ages
TaxaExtant.0.5 <- Map(getAgesExtantSpecies, Taxa.lam.0.5, Tree.lam.05)

clado.df.0.5 <- do.call("rbind", TaxaExtant.0.5)

##diversification
clado.df.0.5$div <- rep(r.vector, each = Ntip)

##turnover rates
clado.df.0.5$turnover <- rep(turnover, each = Ntip)

##speciation rates
clado.df.0.5$lambda <- rep(lam, nrow(clado.df.0.5))

##number of species
clado.df.0.5$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.lam.05, tree.max)
clado.df.0.5$root.age <- rep(root.ages, each = Ntip)
##budding speciation = 0 ; bifurcating = 1
clado.df.0.5$mode <- rep(1, nrow(clado.df.0.5))
##anagenetic = 1; not anagenetic = 0
clado.df.0.5$anag <- rep(0, nrow(clado.df.0.5))
###extinction rates
clado.df.0.5$mu <- rep(mu.vector, each = Ntip)


##renaming columns
clado.df.bif.0.5 <- clado.df.0.5 %>% rename(True.age = Age,
                                          Estimated.age = tip_length_reconstr,
                                          N_sisters = n_sisters_reconstr) 


# trees with lambda = 1 ---------------------------------------------------


###speciation rate
lam <- 1

## turnover rate
turnover <- seq(from = 0, to = 0.99, length.out = 100)

###extinction rate
mu.vector <- lam * turnover

##diversification rates vector
r.vector <- (lam - mu.vector)


##Number of extant species in the simulated tree
Ntip <- 100

##simulating 100 trees with lam = 0.1
Tree.lam.1 <- vector("list", 100)

##for loop due to different extinction rates
for(i in seq_along(mu.vector)){
  Tree.lam.1[i] <- sim.bd.taxa(n = Ntip, lambda = lam, mu = mu.vector[i],
                               complete = TRUE, numbsim = 1)
  print(i)
}

###save Tree with the extinction patterns
save(Tree.lam.1, file = file.path(getwd(), pro, q_error, "Tree.lam.1.RData"))

###simulating species taxonomy; only bifurcating speciation

Taxa.lam.1 <- lapply(Tree.lam.1, sim.taxonomy, beta = 1, lambda.a = 0)

##Calculating estimated and true ages
TaxaExtant.1 <- Map(getAgesExtantSpecies, Taxa.lam.1, Tree.lam.1)

clado.df.1 <- do.call("rbind", TaxaExtant.1)

##diversification
clado.df.1$div <- rep(r.vector, each = Ntip)

##turnover rates
clado.df.1$turnover <- rep(turnover, each = Ntip)

##speciation rates
clado.df.1$lambda <- rep(lam, nrow(clado.df.1))

##number of species
clado.df.1$number.sp <- rep(Ntip, Ntip)

##root ages
root.ages <- sapply(Tree.lam.1, tree.max)
clado.df.1$root.age <- rep(root.ages, each = Ntip)
##budding speciation = 0 ; bifurcating = 1
clado.df.1$mode <- rep(1, nrow(clado.df.1))
##anagenetic = 1; not anagenetic = 0
clado.df.1$anag <- rep(0, nrow(clado.df.1))
###extinction rates
clado.df.1$mu <- rep(mu.vector, each = Ntip)


##renaming columns
clado.df.bif.1 <- clado.df.1 %>% rename(True.age = Age,
                                            Estimated.age = tip_length_reconstr,
                                            N_sisters = n_sisters_reconstr) 
bif.sim <- rbind(clado.df.bif.01, clado.df.bif.0.5, clado.df.bif.1)

write_csv(bif.sim, file = file.path(getwd(), pro, q_error, "bif.sim.csv"))
