library(ape)
library(PBD)
library(picante)


source('UtilisR/calculate_tip_ages.R')


# Simulate tree under protracted speciation model
#################################################
b_1 <- 0.2 # speciation-initiation rate of good species
la_1 <- 1 # speciation-completion rate
b_2 <- 0.2 # speciation-initiation rate of incipient species 
mu_1 <- 0.1 # extinction rate of good species
mu_2 <- 0.1 # extinction rate of incipient species
P <- c(b_1, la_1, b_2, mu_1, mu_2)
SimTreePBD <- pbd_sim(pars = P, age = 15)

plot(SimTreePBD$igtree.extinct, mar = c(5, 0.5, 0.1, 0.1))
axisPhylo()


# Get age of extant species
###########################
L <- SimTreePBD$L
L <- L[L[, 4] != -1 & L[, 5] == -1, ]


# Rename tips in reconstructed phylogeny according to species in object L
SimTreeRecon <- SimTreePBD$recontree
Lbl <- unlist(lapply(strsplit(SimTreeRecon$tip.label, '-'), function(x) x[1]))
Lbl <- gsub('S', '', Lbl)
SimTreeRecon$tip.label <- Lbl

TipAgesRecon <- calculate_tip_ages(SimTreeRecon)

TipAgesRecon <- TipAgesRecon[match(as.character(L[, 6]), TipAgesRecon$sp), ]

TaxaExtant <- data.frame(sp = L[, 6],
                         Age = L[, 4],
                         Start_Incipient = L[, 3],
                         Tip_Length_Reconstr = TipAgesRecon[, 2])
