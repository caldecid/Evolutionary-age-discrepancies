###############################################################################
##########################        Tree simulations ############################
###############################################################################

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

##setting working directory to the correspondent results carpet
setwd(file.path(getwd(), pro, q_error))


######Don't set a seed cause we are going to use the same parameters for 100 trees


# Trees with mu = 0.9 -------------------------------------------------

###extinction rate
mu <- 0.9

## speciation rate
lambda <- 1


##Number of extant species in the simulated tree
Ntip <- 100


##for loop due to different extinction rates
Tree.mu.09 <- sim.bd.taxa(n = Ntip, numbsim = 1000, lambda = lambda, mu= mu)


###save Tree with the extinction patterns
save(Tree.mu.09, file = "Tree.mu.09.RData")

# Trees with mu = 0.5 -------------------------------------------------

###extinction rate
mu <- 0.5

## speciation rate
lambda <- 1


##Number of extant species in the simulated tree
Ntip <- 100

##Simulating trees mu = 0.5
Tree.mu.05 <- sim.bd.taxa(n = Ntip, numbsim = 1000, lambda = lambda, mu= mu)


###save Tree with the extinction patterns
save(Tree.mu.05, file = "Tree.mu.05.RData")

# Trees with mu = 0 -------------------------------------------------

###extinction rate
mu <- 0

## speciation rate
lambda <- 1

##Number of extant species in the simulated tree
Ntip <- 100

##Simulating trees mu = 0
Tree.mu.0 <- sim.bd.taxa(n = Ntip, numbsim = 1000, lambda = lambda, mu= mu)

###save Tree with the extinction patterns
save(Tree.mu.0, file = "Tree.mu.0.RData")

load(file = file.path(getwd(), pro, q_error, "Tree.mu.09.RData"))
