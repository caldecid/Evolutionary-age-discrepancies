########nonrandom sampling and extinction

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

##loading trees from the conservation simulation with high, intermediate and low
##extinction rates
load(file = file.path(getwd(), pro, c_status,"trees.conservation.RData"))

##only using intermediate extinction

##prunning trees
##intermediate extinction
trees.extant.int <- lapply(trees.int, prune.fossil.tips) 
##root ages
root.int <- sapply(trees.int, tree.max)


# Bifurcating speciation --------------------------------------------------

##simulating speciation via bifurcation
tax.bif <- lapply(trees.int, sim.taxonomy, beta = 1)

##get ages extant species
ages.bif <- Map(getAgesExtantSpecies, tax.bif, trees.int, Tol = 1e-6) 


###########nonrandom incomple sampling of extant species#######################

##the nonrandom probability increases with age

##get ages extant species with a nonrandom incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.bif, trees.int, SamplingFrac = 0.25,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 75)

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


##get ages extant species with a nonrandom incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.bif, trees.int, SamplingFrac = 0.50,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 50)

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))

##initial list

#collapsing in a DF
ages.bif <- do.call("rbind", ages.bif) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.bif$root.age <- rep(root.int, each = 100)

##adding rTrue.age
ages.bif <- ages.bif %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.bif$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species
ages.bif$species <- paste0(ages.bif$label,".",
                                     ages.bif$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.bif$status <- discretize(ages.bif$rTrue.age, 
                              method = "frequency", 
                              breaks = 5,
                              labels = c("LC", "NT", "VU", "EN", "CR"))

ages.bif$status <- factor(ages.bif$status,
                            levels = c("LC", "NT", "VU", "EN", "CR"))

##selecting ages.bif.0.25 columns
ages.bif.0.25 <- ages.bif.incomp.25 %>% select(species,  
                                               Incomp.age) %>% 
                rename(Incomp.age.25 = Incomp.age)

##merging ages.bif with ages.bif.0.25
ages.total.0.25 <- left_join(ages.bif, ages.bif.0.25, 
                             by = "species")

##selecting ages.bif.0.5 columns
ages.bif.0.5 <- ages.bif.incomp.50 %>% select(species,
                                              Incomp.age) %>% 
                rename(Incomp.age.50 = Incomp.age)

##merging the three dataframes
ages.total <- left_join(ages.total.0.25, ages.bif.0.5,
                        by = "species")

##Calculating the mean of each age for each tree 
ages.mean <-ages.total %>% 
                  group_by(tree, status) %>% 
                  summarise(mean.true = mean(True.age),
                            mean.phy = mean(Estimated.age),
                            mean.0.25 = mean(Incomp.age.25, na.rm = TRUE),
                            mean.0.50 = mean(Incomp.age.50, na.rm = TRUE)) 

#######ages ranking of correct relationship between age and signal
ages.rank <- ages.mean %>% group_by(tree) %>%
              mutate(true.rank = dense_rank(mean.true),
                     phy.rank = dense_rank(mean.phy),
                     f25.rank = dense_rank(mean.0.25),
                     f50.rank = dense_rank(mean.0.50)) %>% 
              mutate(resp_phy = if_else(true.rank == phy.rank, 1, 0),
                     resp_f25 = if_else(true.rank == f25.rank, 1, 0),
                     resp_f50 = if_else(true.rank == f50.rank, 1, 0))

###comparing phylogenetic age
ages.comparison.phy <- ages.rank %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

##comparing incomplete sampling 25%
ages.comparison.f25 <- ages.rank %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25)) %>% 
  filter(sum_f25 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

##comparing incomplete sampling 50%
ages.comparison.f50 <- ages.rank %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50)) %>% 
  filter(sum_f50 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)


# budding speciation ------------------------------------------------------

##simulating speciation via budding
tax.bud <- lapply(trees.int, sim.taxonomy, beta = 0)

##get ages extant species
ages.bud <- Map(getAgesExtantSpecies, tax.bud, trees.int, Tol = 1e-6) 


###########nonrandom incomple sampling of extant species#######################

##the nonrandom probability increases with age

##get ages extant species with a nonrandom incomple sampling (0.25)
ages.bud.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.bud, trees.int, SamplingFrac = 0.25,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bud.incomp.25 <- do.call("rbind", ages.bud.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bud.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 75)

##species name
ages.bud.incomp.25$species <- paste0(ages.bud.incomp.25$label,".",
                                     ages.bud.incomp.25$tree)

##sampling fraction
ages.bud.incomp.25$fraction <- rep(0.25, nrow(ages.bud.incomp.25))


##get ages extant species with a nonrandom incomple sampling (0.50)
ages.bud.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.bud, trees.int, SamplingFrac = 0.50,
                          AgeDependent = TRUE)

##transforming in a DF
ages.bud.incomp.50 <- do.call("rbind", ages.bud.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##tree numbers
ages.bud.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 50)

##species name
ages.bud.incomp.50$species <- paste0(ages.bud.incomp.50$label,".",
                                     ages.bud.incomp.50$tree)

##sampling fraction
ages.bud.incomp.50$fraction <- rep(0.50, nrow(ages.bud.incomp.50))

##initial list

#collapsing in a DF
ages.bud <- do.call("rbind", ages.bud) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.bud$root.age <- rep(root.int, each = 100)

##adding rTrue.age
ages.bud <- ages.bud %>% mutate(rTrue.age = True.age/root.age)
##tree numbers
ages.bud$tree <- rep(paste0("tree.", 1:1000), each = 100)

##species
ages.bud$species <- paste0(ages.bud$label,".",
                           ages.bud$tree)

##assigning the positive ADE through the IUCN status
####intermediate extinction
ages.bud$status <- discretize(ages.bud$rTrue.age, 
                              method = "frequency", 
                              breaks = 5,
                              labels = c("LC", "NT", "VU", "EN", "CR"))

ages.bud$status <- factor(ages.bud$status,
                          levels = c("LC", "NT", "VU", "EN", "CR"))

##selecting ages.bud.0.25 columns
ages.bud.0.25 <- ages.bud.incomp.25 %>% select(species,  
                                               Incomp.age) %>% 
  rename(Incomp.age.25 = Incomp.age)

##merging ages.bud with ages.bud.0.25
ages.bud.total.0.25 <- left_join(ages.bud, ages.bud.0.25, 
                             by = "species")

##selecting ages.bud.0.5 columns
ages.bud.0.5 <- ages.bud.incomp.50 %>% select(species,
                                              Incomp.age) %>% 
  rename(Incomp.age.50 = Incomp.age)

##merging the three dataframes
ages.bud.total <- left_join(ages.total.0.25, ages.bud.0.5,
                        by = "species")

##Calculating the mean of each age for each tree 
ages.bud.mean <-ages.bud.total %>% 
  group_by(tree, status) %>% 
  summarise(mean.true = mean(True.age),
            mean.phy = mean(Estimated.age),
            mean.0.25 = mean(Incomp.age.25, na.rm = TRUE),
            mean.0.50 = mean(Incomp.age.50, na.rm = TRUE)) 

#######ages ranking of correct relationship between age and signal
ages.bud.rank <- ages.bud.mean %>% group_by(tree) %>%
  mutate(true.rank = dense_rank(mean.true),
         phy.rank = dense_rank(mean.phy),
         f25.rank = dense_rank(mean.0.25),
         f50.rank = dense_rank(mean.0.50)) %>% 
  mutate(resp_phy = if_else(true.rank == phy.rank, 1, 0),
         resp_f25 = if_else(true.rank == f25.rank, 1, 0),
         resp_f50 = if_else(true.rank == f50.rank, 1, 0))

###comparing phylogenetic age
ages.bud.comparison.phy <- ages.bud.rank %>% group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

##comparing incomplete sampling 25%
ages.bud.comparison.f25 <- ages.bud.rank %>% group_by(tree) %>% 
  summarise(sum_f25 = sum(resp_f25)) %>% 
  filter(sum_f25 < 5) %>% 
  count() %>% 
  mutate(error = n/1000*100)

##comparing incomplete sampling 25%
ages.bud.comparison.f50 <- ages.bud.rank %>% group_by(tree) %>% 
  summarise(sum_f50 = sum(resp_f50)) %>% 
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



######bifurcating 


##sampling trees
tree.random <- sample(ages.mean$tree, size = 100)

###fake facets for bifurcating
ages.mean$true.age <- "True Age"
ages.mean$f0 <- "0% missing species"
ages.mean$f25 <- "25% missing species"
ages.mean$f50 <- "50% missing species"


###true age plot
f1 <- ages.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.true +1), group = tree))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.4)+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  ggtitle("Bifurcating speciation")+
  facet_wrap(~true.age)+
  mynamestheme+
  theme(legend.position = "none")


##phylogenetic age
f2 <-ages.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.phy +1), group = tree,
                           color = "#1b9e77"))+
  geom_point(size=3, shape = 21, fill="white", color = "#1b9e77")+
  geom_line(alpha = 0.4, color = "#1b9e77")+
  geom_text(data = ages.comparison.phy, aes(x = 2.0, y = 2.22,
                     label = paste0("Error rate:", " ",error,
                                    "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab("log(Phylogenetic age + 1)")+
  xlab(NULL)+
  facet_wrap(~f0)+
  mynamestheme+
  theme(legend.position = "none")

##Incomplete sampling 25%
f3 <-ages.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.0.25 +1), group = tree,
             color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.4, color = "#d95f02")+
  geom_text(data = ages.comparison.f25, aes(x = 2.0, y = 2.22,
                              label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab("log(Phylogenetic age + 1)")+
  xlab(NULL)+
  facet_wrap(~f25)+
  mynamestheme+
  theme(legend.position = "none")


##Incomplete sampling 50%
f4 <-ages.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.0.50 +1), group = tree,
             color = "#7570b3"))+
  geom_point(size=3, shape = 21, fill="white", color = "#7570b3")+
  geom_line(alpha = 0.4, color = "#7570b3")+
  geom_text(data = ages.comparison.f50, aes(x = 2.0, y = 2.35,
                               label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab("log(Phylogenetic age + 1)")+
  xlab(NULL)+
  facet_wrap(~f50)+
  mynamestheme+
  theme(legend.position = "none")


#####budding

##sampling trees
tree.random <- sample(ages.bud.mean$tree, size = 100)

###fake facets for bifurcating
ages.bud.mean$true.age <- "True Age"
ages.bud.mean$f0 <- "0% missing species"
ages.bud.mean$f25 <- "25% missing species"
ages.bud.mean$f50 <- "50% missing species"


###true age plot
f1.bud <- ages.bud.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.true +1), group = tree))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.4)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ggtitle("Budding speciation")+
  facet_wrap(~true.age)+
  mynamestheme+
  theme(legend.position = "none")


##phylogenetic age
f2.bud <-ages.bud.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.phy +1), group = tree,
             color = "#1b9e77"))+
  geom_point(size=3, shape = 21, fill="white", color = "#1b9e77")+
  geom_line(alpha = 0.4, color = "#1b9e77")+
  geom_text(data = ages.bud.comparison.phy, aes(x = 2.0, y = 2.22,
                                       label = paste0("Error rate:", " ",error,
                                                           "%")),
            inherit.aes = FALSE,
            size = 3.5, 
            family = "serif")+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  facet_wrap(~f0)+
  mynamestheme+
  theme(legend.position = "none")

##Incomplete sampling 25%
f3.bud <-ages.bud.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.0.25 +1), group = tree,
             color = "#d95f02"))+
  geom_point(size=3, shape = 21, fill="white", color = "#d95f02")+
  geom_line(alpha = 0.4, color = "#d95f02")+
  geom_text(data = ages.bud.comparison.f25, aes(x = 2.0, y = 1.9,
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


##Incomplete sampling 50%
f4.bud <-ages.bud.mean %>% filter(tree %in% tree.random) %>% 
  ggplot(aes(x = status, y = log(mean.0.50 +1), group = tree,
             color = "#7570b3"))+
  geom_point(size=3, shape = 21, fill="white", color = "#7570b3")+
  geom_line(alpha = 0.4, color = "#7570b3")+
  geom_text(data = ages.bud.comparison.f50, aes(x = 2.0, y = 2.35,
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
x.bif = grid.arrange(f1, f2, f3, f4, nrow = 4)
x.bud = grid.arrange(f1.bud, f2.bud, f3.bud, f4.bud, nrow = 4)

grid.arrange(x.bif, x.bud, ncol = 2, bottom = bottom)
dev.off()


