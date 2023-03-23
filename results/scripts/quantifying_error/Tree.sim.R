
# Simulating trees for simulating conditions ------------------------------

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



# Trees with lambda = 0.5 -------------------------------------------------

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


# Trees with lambda = 1 ---------------------------------------------------

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


