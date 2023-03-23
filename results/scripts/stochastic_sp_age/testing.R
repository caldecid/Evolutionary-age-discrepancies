
getAgesSpecies <- function(Taxonomy, Phylo, Tol = NULL) {
  # Assign names of tips and nodes from the phylogeny to the species
  Nt <- Ntip(Phylo)
  Names <- c(Phylo$tip, as.character((Nt + 1):(Nt + Phylo$Nnode)))
  Taxonomy$label <- Names[Taxonomy$edge]
  # What looks like end=0 is not exactly 0, so we need a tolerance 
  # to identify extant species
  if (is.null(Tol)) {
    Tol <- max(branching.times(Phylo)) / 1e6
  }
  # Keep only extant species
  ExtantSpecies <- unique(Taxonomy$sp[Taxonomy$end < Tol])
  TaxaExtantTmp <- Taxonomy[Taxonomy$sp %in% ExtantSpecies, ]
  # Order so that for each species the edge ending at the present is on top of older branches
  TaxaExtantTmp <- TaxaExtantTmp[order(TaxaExtantTmp$sp, TaxaExtantTmp$end), ]
  # A species may appear multiple times and we look for the max of start,
  # look-up the maximum of start
  TaxaExtant <- aggregate(start ~ sp, data = TaxaExtantTmp, FUN = max)
  colnames(TaxaExtant)[2] <- 'Age'
  TaxaExtant <- data.frame(TaxaExtant, TaxaExtantTmp[TaxaExtantTmp$end < Tol, -1])
  return(TaxaExtant)
}


# Simulating trees for experiment -----------------------------------------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

##simulate Tree with 100 Ntips

tree.0 <- sim.bd.taxa(n = 100, lambda = 0.1, mu = 0.05,
                                        complete = TRUE, numbsim = 1)[[1]]

tree.1 <- phytools::rotateNodes(tree.0, nodes = "all")

##nodelabes
node.0 <- (Ntip(tree.0)+1):(Ntip(tree.0)+Nnode(tree.0))

##vector of 100 trees
trees <- vector("list", length = 100)





###only budding speciation
taxas <- vector("list", length(trees))

for(i in 1:length(trees)){
  
  taxas[[i]] <- tryCatch(sim.taxonomy(trees[[i]], beta = 0, lambda.a = 0),
                         error = function(e) NA)
}  


##calculating true age
taxa.age <- Map(getAgesSpecies, taxas, trees)

taxa.df <- do.call("rbind", taxa.age)

##defining the tree
taxa.df$tree <- rep(paste0("tree",1:100), each = 100)

##pivot wider for analysis
taxa.wide <- taxa.df %>% select(sp, Age, start, label, tree) %>% 
            pivot_wider(names_from = tree, values_from = Age)

##determining the variance for each sp
var.trees <- taxa.wide %>% select(-c(sp, start, label)) %>%
              apply(MARGIN = 1, var)

taxa.wide$var <- var.trees

##filtering the higher variances

taxa.wide.high <- taxa.wide %>% filter(var > 150)

##calling them in a graphic
i <- which(trees[[30]]$tip.label %in% taxa.wide.high$label)
##it is the same number among trees
##plotting
plot(tree.0, show.tip.label = FALSE)
tiplabels(tree.0$tip.label[i], i, adj = 0)


##pivoting again
taxa.wide.high$label <- as.factor(taxa.wide.high$label)

taxa.long.high <- taxa.wide.high %>% select(-c(var, start)) %>% 
                  pivot_longer(starts_with("tree"),
                               names_to = "trees", 
                               values_to = "age")

ggplot(taxa.long.high, aes(x = label, y = age))+
   geom_point()+
  xlab("Species")

m1 <- aov(data = taxa.long.high, age ~ label)
anova(m1)
summary(m1)

TukeyHSD(m1, conf.level=.95)

#####working in manipulations

####histogram of species

##simulating the expected distribution of ages
mu <- 0.05
ages <- seq(0, 80, length.out = 100)
expected <- dexp(ages, rate = mu)
dist_ages <- as.data.frame(cbind(ages, expected))

hist(rexp(20, rate = mu), freq = FALSE) # Your species ages
lines(age, expected)


ggplot(taxa.long.high, aes(x = age))+
  geom_histogram(aes(y=..density..), position = "dodge") +
  #geom_histogram(position = "identity")+
  geom_line(data = dist_ages, aes(x = ages, y = expected), color = "red",
            size = 1)+
  facet_wrap(~label)+
  xlab("Species ages")+
  ylab("Density")

#############Rotating 20 % of the nodes########################

##vector of 100 trees
trees.20 <- vector("list", length = 100)


##nested for loop
for(i in 1:length(trees.20)){
  Done <- FALSE 
  while (!Done) {
    x <- sample(node.0, size = round(0.7*length(node.0)), replace = TRUE)
    trees.20[[i]] <- tree.0
    success_rotate <- 0
    for(j in 1:length(x)){
      rotated_tree <- tryCatch(ape::rotate(trees.20[[i]], node = x[j]),
                               error = function(e) NA, warning = function(w) NA)
      if (!all(is.na(rotated_tree))) {
        trees.20[[i]] <- rotated_tree
        success_rotate <- success_rotate + 1
      }
    }
    if (success_rotate == length(x)) {
      png_null_device(4, 4)
      plot_ok <- tryCatch(plot.phylo(trees.20[[i]], plot = FALSE),
                          error = function(e) NA, warning = function(w) NA)
      dev.off()
      if (all(!is.na(plot_ok))) {
        # All attempts to rotate where successfull and phylogeny can be plotted
        Done <- TRUE 
      }
    }
  }
}





###only budding speciation
taxas.20 <- vector("list", length(trees.20))

for(i in 1:length(trees.20)){
  
  taxas.20[[i]] <- tryCatch(sim.taxonomy(trees.20[[i]], beta = 0, lambda.a = 0),
                         error = function(e) NA)
}  


##calculating true age
taxa.age.20 <- Map(getAgesSpecies, taxas.20, trees.20)

taxa.df.20 <- do.call("rbind", taxa.age.20)

##defining the tree
taxa.df.20$tree <- rep(paste0("tree",1:100), each = 100)

##pivot wider for analysis
taxa.wide.20 <- taxa.df.20 %>% select(sp, Age, start, label, tree) %>% 
  pivot_wider(names_from = tree, values_from = Age)

##determining the variance for each sp
var.trees.20 <- taxa.wide.20 %>% select(-c(sp, start, label)) %>%
  apply(MARGIN = 1, var)

taxa.wide.20$var <- var.trees.20

##filtering the higher variances

taxa.wide.high.20 <- taxa.wide.20 %>% filter(var > 150)

##calling them in a graphic
i <- which(trees.20[[30]]$tip.label %in% taxa.wide.high.20$label)
##it is the same number among trees
##plotting
plot(tree.0, show.tip.label = FALSE)
tiplabels(tree.0$tip.label[i], i, adj = 0)


##pivoting again
taxa.wide.high.20$label <- as.factor(taxa.wide.high.20$label)

taxa.long.high.20 <- taxa.wide.high.20 %>% select(-c(var, start)) %>% 
  pivot_longer(starts_with("tree"),
               names_to = "trees", 
               values_to = "age")

##graphics
ggplot(taxa.long.high.20, aes(x = age))+
  geom_histogram(aes(y=..density..), position = "dodge") +
  #geom_histogram(position = "identity")+
  geom_line(data = dist_ages, aes(x = ages, y = expected), color = "red",
            size = 1)+
  facet_wrap(~label)+
  xlab("Species ages")+
  ylab("Density")



#########another exploration

#two histograms (including the red line of expected ages):
  #One histogram with ages of all extinct species and one for the extant species?

Nt <- Ntip(tree.0)
Names <- c(tree.0$tip, as.character((Nt + 1):(Nt + tree.0$Nnode)))
taxa.0$label <- Names[taxa.0$edge]
# What looks like end=0 is not exactly 0, so we need a tolerance 
# to identify extant species
#if (is.null(Tol)) {
Tol <- max(branching.times(tree.0)) / 1e6
#}
# extant species
ExtantSpecies <- unique(taxa.0$sp[taxa.0$end < Tol])
TaxaExtantTmp <- taxa.0[taxa.0$sp %in% ExtantSpecies, ]
#Extinct species  (merely exclude the Extant)
TaxaExtinctTmp <- taxa.0[-which(taxa.0$sp %in% ExtantSpecies), ]
# Order so that for each species the edge ending at the present or most inmediate past
#is on top of older branches
TaxaExtantTmp <- TaxaExtantTmp[order(TaxaExtantTmp$sp, TaxaExtantTmp$end), ]
TaxaExtinctTmp <- TaxaExtinctTmp[order(TaxaExtinctTmp$sp, TaxaExtinctTmp$end), ]
# A species may appear multiple times and we look for the max of start,
# look-up the maximum of start
TaxaExtant <- aggregate(start ~ sp, data = TaxaExtantTmp, FUN = max)
colnames(TaxaExtant)[2] <- 'Age'

TaxaExtant <- data.frame(TaxaExtant, TaxaExtantTmp[TaxaExtantTmp$end < Tol, -1])
TaxaExtant$status <- rep("Extant", nrow(TaxaExtant))
#for extinct species
TaxaExtinct <- aggregate(start ~ sp, data = TaxaExtinctTmp, FUN = max)
colnames(TaxaExtinct[2]) <- 'Age'
Tax.ext <- TaxaExtinctTmp %>% filter(grepl('t', label))
TaxaExtinct.1 <- data.frame(TaxaExtinct, Tax.ext[,-1])
TaxaExtinct.1$Age <- TaxaExtinct.1$start - TaxaExtinct.1$end
TaxaExtinct.1$status <- rep("Extinct", nrow(TaxaExtinct.1))

##dataframes
df.extinct <- TaxaExtinct.1 %>% select(label, Age, status)
df.extant <- TaxaExtant %>% select(label, Age, status)
df.ages <- rbind(df.extinct, df.extant)

##########Figures#####

df.ages %>% 
  ggplot(aes(x = Age))+
  geom_histogram(aes(y=..density..), position = "dodge",
                 bins = 30)+
  geom_line(data = dist_ages, aes(x = ages, y = expected), color = "red",
            size = 1)+
  facet_wrap(~status)+
  xlab("Species ages")+
  ylab("Density")



# rotating 100% of the nodes ----------------------------------------------

##list of 100 trees
trees.100 <- vector("list", length = 100)


for(i in 1:length(trees.100)){
  trees.100[[i]] <- rotateNodes(tree.0)
  print(paste0("tree.",i))
}

##taxa list
taxas.100 <- vector("list", length(trees.100))

for(i in 1:length(trees.100)){
  
  taxas.100[[i]] <- tryCatch(sim.taxonomy(trees.100[[i]], beta = 0, lambda.a = 0),
                             error = function(e) NA)
}  


##calculating true age
taxa.age.100 <- Map(getAgesExtantSpecies, taxas.100, trees.100)

taxa.df.100.2 <- do.call("rbind", taxa.age.100)

##defining the tree
taxa.df.100$tree <- rep(paste0("tree",1:100), each = 100)

##pivot wider for analysis
 taxa.wide.100 <- taxa.df.100 %>% select(Age, start, label, tree) %>% 
                  pivot_wider(names_from = tree, values_from = Age)




##determining the variance for each sp
var.trees.100 <- taxa.wide.100 %>% select(-c(start, label)) %>%
  apply(MARGIN = 1, var)

taxa.wide.100$var <- var.trees.100

###obtaining labels of the highest variance ages
high.sp <- taxa.wide.100 %>% slice_max(var, n = 4) %>% select(label)
high.sp <- as.vector(high.sp$label)

age.high <- df.ages[df.ages$label %in% high.sp,]

##obtaining labels of the lowest variance ages
low.sp <- taxa.wide.100 %>% slice_min(var, n = 4) %>% select(label)
low.sp <- as.vector(low.sp$label)

age.low <- df.ages[df.ages$label %in% low.sp,]

### determining the 4 highest variances

taxa.wide.100 %>% slice_max(var, n = 4) %>%
  select(-c(var, start)) %>% 
  pivot_longer(starts_with("tree"),
               names_to = "trees", 
               values_to = "age") %>% 
  ggplot(aes(x = age))+
  geom_histogram(aes(y=..density..), position = "dodge") +
  #geom_histogram(position = "identity")+
  geom_line(data = dist_ages, aes(x = ages, y = expected), color = "red",
            size = 1)+
  geom_vline(data = age.high, aes(xintercept = Age), size = 1,
             color = "green")+
  facet_wrap(~label)+
  xlab("Species ages")+
  ylab("Density")+
  ggtitle("Species with highest age variance")

#####determining the 4 lowest variances
taxa.wide.100 %>% slice_min(var, n = 4) %>%
  select(-c(var, start)) %>% 
  pivot_longer(starts_with("tree"),
               names_to = "trees", 
               values_to = "age") %>% 
  ggplot(aes(x = age))+
  geom_histogram(aes(y=..density..), position = "dodge") +
  #geom_histogram(position = "identity")+
  geom_line(data = dist_ages, aes(x = ages, y = expected), color = "red",
            size = 1)+
  geom_vline(data = age.low, aes(xintercept = Age), size = 1,
             color = "green")+
  facet_wrap(~label)+
  xlab("Species ages")+
  ylab("Density")+
  ggtitle("Species with lowest age variance")

