
# uniform incomplete sampling ---------------------------------------------


##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

###creating new folder
path <- "results/data/processed"
newdir <- "incomplete_sampling"
dir.create(file.path(path, newdir))



# Bifurcating speciation --------------------------------------------------

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


######## 0% bif manipulation
ages.0.bif <- do.call("rbind", ages.0.bif) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.0.bif$root.age <- rep(root.25, each = 100)

##tree numbers
ages.0.bif$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.0.bif <- ages.0.bif %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age)

##extinction
ages.0.bif$extinction <- rep("intermediate", nrow(ages.0.bif))

##species name
ages.0.bif$species <- paste0(ages.0.bif$label,".",
                                     ages.0.bif$tree)

##sampling fraction
ages.0.bif$fraction <- rep(0, nrow(ages.0.bif))


###speciation mode
ages.0.bif$speciation <- rep("bif", nrow(ages.0.bif))



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


###############################budding speciation#########################

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

###0% budding
ages.0.bud <- do.call("rbind", ages.0.bud) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  select(label, True.age, Estimated.age)

##root age
ages.0.bud$root.age <- rep(root.25, each = 100)

##tree numbers
ages.0.bud$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.0.bud <- ages.0.bud %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age)

##extinction
ages.0.bud$extinction <- rep("intermediate", nrow(ages.0.bud))

##species name
ages.0.bud$species <- paste0(ages.0.bud$label,".",
                             ages.0.bud$tree)

##sampling fraction
ages.0.bud$fraction <- rep(0, nrow(ages.0.bud))


###speciation mode
ages.0.bud$speciation <- rep("bud", nrow(ages.0.bud))

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bud.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bud, trees.25.missing,
                          SamplingFrac = 0.25)


ages.bud.incomp.25 <- do.call("rbind", ages.bud.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bud.incomp.25$root.age <- rep(root.25, each = 100)

##tree numbers
ages.bud.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)



##scale ages
ages.bud.incomp.25 <- ages.bud.incomp.25 %>%
                      mutate(rTrue.age = True.age/root.age,
                             rPhylo.age = Estimated.age/root.age,
                             rIncomp.age = Incomp.age/root.age)

##extinction
ages.bud.incomp.25$extinction <- rep("intermediate", nrow(ages.bud.incomp.25))

##species name
ages.bud.incomp.25$species <- paste0(ages.bud.incomp.25$label,".",
                                     ages.bud.incomp.25$tree)

##sampling fraction
ages.bud.incomp.25$fraction <- rep(0.25, nrow(ages.bud.incomp.25))


###speciation mode
ages.bud.incomp.25$speciation <- rep("bud", nrow(ages.bud.incomp.25))

#####sampling fraction = 0.50##########

##get ages extant species with an uniform incomple sampling (0.50)
ages.bud.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bud, trees.50.missing,
                          SamplingFrac = 0.50)


ages.bud.incomp.50 <- do.call("rbind", ages.bud.incomp.50) %>% 
              rename(True.age = Age,
                 Estimated.age = tip_length_reconstr,
                 Incomp.age = IncompletePhyloAge) %>%
          select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bud.incomp.50$root.age <- rep(root.50, each = 100)

##tree numbers
ages.bud.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bud.incomp.50 <- ages.bud.incomp.50 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bud.incomp.50$extinction <- rep("intermediate", nrow(ages.bud.incomp.50))

##species name
ages.bud.incomp.50$species <- paste0(ages.bud.incomp.50$label,".",
                                     ages.bud.incomp.50$tree)

##sampling fraction
ages.bud.incomp.50$fraction <- rep(0.50, nrow(ages.bud.incomp.50))

###speciation mode
ages.bud.incomp.50$speciation <- rep("bud", nrow(ages.bud.incomp.50))

###mergin data frames
ages.incomplete <- rbind(ages.bif.incomp.25,
                         ages.bif.incomp.50,
                         ages.bud.incomp.25,
                         ages.bud.incomp.50)

write.csv(ages.incomplete,
    file = "results/data/processed/incomplete_sampling/ages.incomplete.uniform.csv")

ages.incomplete <- read_csv("results/data/processed/incomplete_sampling/ages.incomplete.uniform.csv")

#########figures################################

# figures -----------------------------------------------------------------

## themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (11),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (11),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (10)),
                      legend.background = element_rect(fill = "gray90",
                                               size = 0.5, linetype = "dotted"),
                      legend.position = "bottom")

###bifurcating 0% incomplete sampling
ages.average.0.bif <- ages.0.bif %>% group_by(tree) %>% 
                      summarise(mape.0 =
                            mean(abs(True.age - Estimated.age)/True.age)*100)
##the mean and sd of the error
ages.average.0.bif %>% summarise(mean.phy = mean(mape.0),
                                    sd.phy = sd(mape.0))


##bifurcating 25% incomplete sampling
ages.average.25.bif <- ages.incomplete %>% filter(fraction == 0.25,
                                                    speciation == "bif") %>% 
        group_by(tree) %>% 
        summarise(mape.0.25 = mean(abs(True.age - Incomp.age)/True.age)*100)

##the mean and sd of the error
ages.average.25.bif %>% summarise(mean.0.25 = mean(mape.0.25),
                                    sd.0.25 = sd(mape.0.25))

##bifurcating 50% incomplete sampling
ages.average.50.bif <- ages.incomplete %>% filter(fraction == 0.50,
                                                   speciation == "bif") %>% 
          group_by(tree) %>% 
          summarise(mape.0.50 = 
                      mean(abs(True.age - Incomp.age)/True.age)*100)

##calculating the mean and sd of the error
ages.average.50.bif %>% summarise(mean.0.5 = mean(mape.0.50),
                                    sd.0.5 = sd(mape.0.50))

###merging mape dataframe
ages.mape.bif <- left_join(ages.average.0.bif,
                           ages.average.25.bif,
                           by = "tree") %>% 
                  left_join(ages.average.50.bif,
                            by = "tree") %>% 
                  pivot_longer(cols = starts_with("mape."),
                               values_to = "mape", names_to = "estimate")

ages.mape.bif$speciation <- "bif"
                          

##budding 

##0% incomplete sampling
ages.average.0.bud <- ages.0.bud %>% group_by(tree) %>% 
                      summarise(mape.0 
                              = mean(abs(True.age - Estimated.age)/True.age)*100)

##summary
ages.average.0.bud %>% summarise(mean.mape = mean(mape.0),
                                 sd.mape = sd(mape.0))

##25% incomplete sampling
ages.average.25.bud <-ages.incomplete %>% filter(fraction == 0.25,
                                                   speciation == "bud") %>% 
                        group_by(tree) %>% 
          summarise(mape.0.25 = mean(abs(True.age - Incomp.age)/True.age)*100)

##calculating the mean and sd of the error
ages.average.25.bud %>% summarise(mean.0.25 = mean(mape.0.25),
                                    sd.0.25 = sd(mape.0.25))
##budding 50% incomplete sampling
ages.average.50.bud <- ages.incomplete %>% filter(fraction == 0.50,
                                                   speciation == "bud") %>% 
        group_by(tree) %>%  
        summarise(mape.0.50 = mean(abs(True.age - Incomp.age)/True.age)*100)

##calculating the mean and sd of the error
ages.average.50.bud %>% summarise(mean.0.5 = mean(mape.0.50),
                                    sd.0.5 = sd(mape.0.50))

###merging mape dataframe
ages.mape.bud <- left_join(ages.average.0.bud, ages.average.25.bud,
                           by = "tree") %>% 
                 left_join(ages.average.50.bud, by = "tree") %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.bud$speciation <- "bud"

##mergind dataframes
ages.mape <- rbind(ages.mape.bif, ages.mape.bud)

ages.mape$estimate <- as.factor(ages.mape$estimate)

ages.mape$estimate <- factor(ages.mape$estimate, levels = c("mape.0",
                                                            "mape.0.25",
                                                            "mape.0.50"))
###PLOT
png("text/figures/MAPE.incomplete.png", 
    width = 17, height = 15, units = "cm", 
    pointsize = 8, res = 300)


ggplot(ages.mape, aes(y = mape, x = estimate, fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  #geom_jitter(size = 0.1, alpha = 0.3)+
  facet_wrap(~speciation,  labeller = as_labeller(c(bif = "Bifurcating",
                                                    bud = "Budding")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                     name="Missing species",
                     breaks=c("mape.0", "mape.0.25", "mape.0.50"),
                     labels=c("0%", "25%",
                              "50%"))+
  ylim(0, 1000)+
  ylab("MAPE (%)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())

dev.off()
