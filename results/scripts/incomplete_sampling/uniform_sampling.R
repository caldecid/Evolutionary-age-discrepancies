
# uniform incomplete sampling ---------------------------------------------


##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

###creating new folder
path <- "results/data/processed"
newdir <- "incomplete_sampling"
dir.create(file.path(path, newdir))


##loading trees from the conservation simulation with high, intermediate and low
##extinction rates
load(file = file.path(getwd(), pro, c_status,"trees.conservation.RData"))

##only using intermediate extinction

##prunning trees
##intermediate extinction
trees.extant.int <- lapply(trees.int, prune.fossil.tips) 
##root ages
root.int <- sapply(trees.int, tree.max)
##simulating speciation via bifurcation
tax.bif <- lapply(trees.int, sim.taxonomy, beta = 1)

##get ages extant species
ages.bif <- Map(getAgesExtantSpecies, tax.bif, trees.int, Tol = 1e-6) 

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.bif, trees.int, SamplingFrac = 0.25)


ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.25$root.age <- rep(root.int, each = 75)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 75)

##lambda
ages.bif.incomp.25$lambda <- rep(0.3, each = 75)

##mu
ages.bif.incomp.25$mu <- rep(0.15, each = 75)

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
                          ages.bif, trees.int, SamplingFrac = 0.50)


ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.50$root.age <- rep(root.int, each = 50)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 50)

##lambda
ages.bif.incomp.50$lambda <- rep(0.3, each = 50)

##mu
ages.bif.incomp.50$mu <- rep(0.15, each = 50)

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

##simulating speciation via budding
tax.bud <- lapply(trees.int, sim.taxonomy, beta = 0)

##get ages extant species
ages.bud <- Map(getAgesExtantSpecies, tax.bud, trees.int, Tol = 1e-6) 

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bud.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.bud, trees.int, SamplingFrac = 0.25)


ages.bud.incomp.25 <- do.call("rbind", ages.bud.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bud.incomp.25$root.age <- rep(root.int, each = 75)

##tree numbers
ages.bud.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 75)

##lambda
ages.bud.incomp.25$lambda <- rep(0.3, each = 75)

##mu
ages.bud.incomp.25$mu <- rep(0.15, each = 75)

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
                          ages.bud, trees.int, SamplingFrac = 0.50)


ages.bud.incomp.50 <- do.call("rbind", ages.bud.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bud.incomp.50$root.age <- rep(root.int, each = 50)

##tree numbers
ages.bud.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 50)

##lambda
ages.bud.incomp.50$lambda <- rep(0.3, each = 50)

##mu
ages.bud.incomp.50$mu <- rep(0.15, each = 50)

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



##bifurcating 0.25 sampling
ages.average.0.25.bif <- ages.bif.incomp.25 %>% 
 ### group_by(speciation, fraction, tree) %>% 
  summarise(mape.phy = mean(abs(True.age - Estimated.age)/True.age)*100,
            mape.0.25.bif = mean(abs(True.age - Incomp.age)/True.age)*100)

##bifurcating 0.50 sampling
ages.average.0.5.bif <- ages.bif.incomp.50 %>% 
          group_by(tree) %>% 
          summarise(mape.phy = mean(abs(True.age - Estimated.age)/True.age)*100,
                  mape.0.50.bif = mean(abs(True.age - Incomp.age)/True.age)*100)

##budding 0.25 sampling
ages.average.0.25.bud <- ages.bud.incomp.25 %>% 
              group_by(tree) %>% 
              summarise(mape.phy = mean(abs(True.age - Estimated.age)/True.age)*100,
                        mape.0.25.bud = mean(abs(True.age - Incomp.age)/True.age)*100)

##budding 0.50 sampling
ages.average.0.5.bud <- ages.bud.incomp.50 %>% 
        group_by(tree) %>%  
        summarise(mape.phy = mean(abs(True.age - Estimated.age)/True.age)*100,
                  mape.0.50.bud = mean(abs(True.age - Incomp.age)/True.age)*100)


%>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

###PLOT
png("text/figures/MAPE.bif.median.png", 
    width = 15, height = 15, units = "cm", 
    pointsize = 8, res = 300)


ggplot(ages.average, aes(y = mape, x = estimate, fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~extinction,  labeller = as_labeller(c(high = "High extinction",
                                                    intermediate = "Intermediate extinction",
                                                    low = "Low extinction")))+
  scale_fill_discrete(name = "Estimation",
                      labels = c("Expected", "Median", "Phylogenetic"))+
  ylim(0, 150)+
  ylab("MAPE (%)")+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())

dev.off()
