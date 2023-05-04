##############################################################################
#################      Taxa simulations budding speciation ###################
###################          for comparisons             ####################
#############################################################################


##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

# mu = 0.9 ----------------------------------------------------------------

load(file.path(getwd(), pro, q_error, "Tree.mu.09.RData"))

##simulating species taxonomy; only budding speciation. Trees mu = 0.9
taxa.mu.09 <- lapply(Tree.mu.09, sim.taxonomy, beta = 0, lambda.a = 0)

save(taxa.mu.09, file = file.path(getwd(), pro, q_error, "taxa.mu.09.RData"))

##Calculating estimated and true ages
taxa.extant.09 <- Map(getAgesExtantSpecies, taxa.mu.09, Tree.mu.09)

##colapsing list into a dataframe
clado.df.09 <- do.call("rbind", taxa.extant.09)

##root ages
root.ages.09 <- sapply(Tree.mu.09, tree.max)

clado.df.09$root.age <- rep(root.ages.09, each = 100)

##extinction rates
clado.df.09$mu <- rep(0.9, nrow(clado.df.09))

##tree
clado.df.09$tree <- rep(paste0("tree.09.", seq(1, 1000, 1)), each = 100)

##budding speciation = 0 ; bifurcating = 1
clado.df.09$mode <- rep("Budding", nrow(clado.df.09))


##manipulating the data frame
clado.df.bud.09 <- clado.df.09 %>% select(label, Age, tip_length_reconstr, mode,
                                                        mu, root.age, tree) %>% 
                                rename(species = label,
                                       True.age = Age,
                                       Estimated.age = tip_length_reconstr) %>%
                                mutate(rTrue.age = True.age/root.age,
                                       rEstimated.age = Estimated.age/root.age,
                                       accuracy = rEstimated.age - rTrue.age)

##defining the intensity of extinction
clado.df.bud.09$extinction <- rep("High", nrow(clado.df.bud.09))

write_csv(taxa.mu.09,#clado.df.bud.09,
          file = file.path(getwd(), pro, q_error, "taxa.mu.09.csv"))

##read data
clado.df.bud.09 <- read_csv(file = file.path(getwd(), pro, q_error, "taxa.mu.09.csv"))


# mu = 0.5 ----------------------------------------------------------------

load(file.path(getwd(), pro, q_error, "Tree.mu.05.RData"))
##simulating species taxonomy; only budding speciation. Trees mu = 0.5

taxa.mu.05 <- lapply(Tree.mu.05, sim.taxonomy, beta = 0, lambda.a = 0)

save(taxa.mu.05, file = file.path(getwd(), pro, q_error, "taxa.mu.05.RData"))

##Calculating estimated and true ages
taxa.extant.05 <- Map(getAgesExtantSpecies, taxa.mu.05, Tree.mu.05)

##colapsing list into a dataframe
clado.df.05 <- do.call("rbind", taxa.extant.05)

##speciation rates
clado.df.05$lambda <- rep(1, nrow(clado.df.05))

##extinction rates
clado.df.05$mu <- rep(0.5, nrow(clado.df.05))

##root ages
root.ages.05 <- sapply(Tree.mu.05, tree.max)

clado.df.05$root.age <- rep(root.ages.05, each = 100)

##tree
clado.df.bud.05$tree <- rep(paste0("tree.05.", seq(1, 1000, 1)), each = 100)
##budding speciation = 0 ; bifurcating = 1
clado.df.05$mode <- rep("Budding", nrow(clado.df.05))


##manipulating the data frame
clado.df.bud.05 <- clado.df.05 %>% select(label, Age, tip_length_reconstr, mode,
                                          tree, mu, root.age) %>% 
                                rename(species = label,
                                       True.age = Age,
                                       Estimated.age = tip_length_reconstr) %>%
                                mutate(rTrue.age = True.age/root.age,
                                        rEstimated.age = Estimated.age/root.age,
                                        accuracy = rEstimated.age - rTrue.age)

##intensity of extinction
clado.df.bud.05$extinction <- rep("Intermediate", nrow(clado.df.bud.05))

##saving
write_csv(clado.df.bud.05, 
          file = file.path(getwd(), pro, q_error, "taxa.mu.05.csv"))

##reading 
clado.df.bud.05 <- read_csv(file = file.path(getwd(), pro, q_error,
                                            "taxa.mu.05.csv"))
                                            

# mu = 0 ----------------------------------------------------------------

load(file.path(getwd(), pro, q_error, "Tree.mu.0.RData"))
##simulating species taxonomy; only budding speciation. Trees mu = 0
taxa.mu.0 <- lapply(Tree.mu.0, sim.taxonomy, beta = 0, lambda.a = 0)

##saving taxonomy list
save(taxa.mu.0, file = file.path(getwd(), pro, q_error, "taxa.mu.0.RData"))

##Calculating estimated and true ages
taxa.extant.0 <- Map(getAgesExtantSpecies, taxa.mu.0, Tree.mu.0)

##colapsing list into a dataframe
clado.df.0 <- do.call("rbind", taxa.extant.0)

##extinction rates
clado.df.0$mu <- rep(0, nrow(clado.df.0))

##root ages
root.ages.0 <- sapply(Tree.mu.0, tree.max)

clado.df.0$root.age <- rep(root.ages.0, each = 100)

##tree
clado.df.0$tree <- rep(paste0("tree.0.", seq(1, 1000, 1)), each = 100)

##budding speciation = 0 ; bifurcating = 1
clado.df.0$mode <- rep("Budding", nrow(clado.df.0))


##manipulating the data frame
clado.df.bud.0 <- clado.df.0 %>% select(label, Age, tip_length_reconstr, mode,
                                          tree ,mu, root.age) %>% 
                                rename(species = label,
                                       True.age = Age,
                                       Estimated.age = tip_length_reconstr) %>%
                                mutate(rTrue.age = True.age/root.age,
                                       rEstimated.age = Estimated.age/root.age,
                                      accuracy = rEstimated.age - rTrue.age)
##intensity of extinction
clado.df.bud.0$extinction <- rep("no_ext", nrow(clado.df.bud.0))

write_csv(clado.df.bud.0,
          file = file.path(getwd(), pro, q_error, "taxa.mu.0.csv"))


##reading 
clado.df.bud.0 <- read_csv(file = file.path(getwd(), pro, q_error,
                                            "taxa.mu.0.csv"))
# Merging data frames -----------------------------------------------------

df.ext.levels <- rbind(clado.df.bud.09, clado.df.bud.05, clado.df.bud.0)

###measuring accuracy
df.ext.levels <- df.ext.levels %>% mutate(accuracy = rEstimated.age - rTrue.age)


##writing
write_csv(df.ext.levels, file = file.path(getwd(), pro, q_error,
                                  "taxa.ext.lev.csv"))

##reading
df.ext.levels <- read_csv(file = file.path(getwd(), pro, q_error,
                                           "taxa.ext.lev.csv"))

df.ext.levels$mode <- rep("Budding", nrow(df.ext.levels))

##needed transformations
df.ext.levels$extinction <- as.factor(df.ext.levels$extinction)

df.ext.levels$extinction <- factor(df.ext.levels$extinction, 
                               levels = c("no_ext", "Intermediate", "High"))



# Figures -----------------------------------------------------------------

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
                      legend.background = element_rect(fill="gray90",
                                                       size=.5, linetype="dotted"))

##histograms for scaled True age and Estimate age
age.ext <- df.ext.levels %>% pivot_longer(cols = c("rTrue.age", "rEstimated.age"),
                                          names_to = "type", values_to = "age")
##type as factor
age.ext$type <- factor(age.ext$type, levels = c("rTrue.age", "rEstimated.age"),
                       labels = c("rTrue.age", "rEstimated.age"))

##png object
png("text/figures/Figure3.facet.True.phylo.age.bud.png",
    width = 14, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f3 <- age.ext %>% 
      ggplot(aes(age))+
      geom_histogram(position = "dodge",breaks = seq(0, 1, by = 0.01))+
      facet_grid(extinction ~ type, scales = "free_y",
             switch = "y", 
             labeller = as_labeller(c(rTrue.age = "True age",
                                      rEstimated.age = "Phylogenetic age",
                                      High = "High",
                                      Intermediate = "Intermediate",
                                      no_ext = "No extinction")))+
             theme_bw()+
  #xlim(0, 1)+
  xlab("Ages (scaled)")+
  ylab("Species Count")+
  ggtitle("Budding speciation")
  


f3 + mynamestheme

dev.off()

####histograms for accuracy (rTrue.Age - rEstimated.age)

##determining how many species are 100% accurate
 accurate.age <- df.ext.levels %>%
  mutate(abs.accuracy = abs(accuracy)) %>% 
  filter(abs.accuracy < 1e-15) %>% group_by(extinction) %>%
  count()

 
######figure displaying a line and number of the accurate species
 ##png object
png("text/figures/Figure4.ratio.age.bud.png",
    width = 12, height = 15, units = "cm", 
     pointsize = 8, res = 300)
 
f4 <- df.ext.levels %>%
  mutate(abs.accuracy = abs(accuracy)) %>% 
  filter(abs.accuracy > 1e-15) %>% 
  ggplot(aes(x = accuracy))+
  geom_histogram(position = "dodge", breaks = seq(-0.3, 0.6, by = 0.02))+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = n),
               data = accurate.age,
               colour = "red")+
  geom_point(aes(x = 0, y = n), data = accurate.age, colour = "red")+
  geom_text(data = accurate.age, aes(x = 0.25, y = 41100,
                     label = paste("Correct estimation:", round(n/1000),"%")),
            size = 3, family = "serif")+
  xlim(-0.4, 0.5)+
  #ylim(0, 5000)+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(High = "High extinction",
                                      Intermediate = "Intermediate extinction",
                                      no_ext = "No extinction")))+
  theme_bw()+
  xlab("Phylogenetic age - True age (scaled)")+
  ylab("Species count")+
  ggtitle("Budding speciation")
  

f4 + mynamestheme

dev.off()


##distribution of ratio between youngest and oldes sp for rphyletic, selecting 
#species and calculate the ratio with the true ages
sp.max <- df.ext.levels %>% select(species, mode, root.age, rTrue.age, rEstimated.age,
                                   extinction) %>% 
                            group_by(root.age, extinction) %>% 
                        filter(rEstimated.age == max(rEstimated.age))     
##eliminate duplicated rows (there are at least two max)
sp.max.1 <- sp.max[-which(duplicated(sp.max$root.age)),] %>% 
                     rename(maxEst = rEstimated.age,
                            maxTrue = rTrue.age)

##defining the sp with the minimun ages
sp.min <- df.ext.levels %>% select(species, root.age, rTrue.age, rEstimated.age,
                                   extinction) %>% 
                            group_by(root.age, extinction) %>% 
                            filter(rEstimated.age == min(rEstimated.age))    

##eliminate the duplicated rows (there are at least two min)
sp.min.1 <- sp.min[-which(duplicated(sp.min$root.age)),] %>% 
                                 rename(sp = species, 
                                        minEst = rEstimated.age,
                                        minTrue = rTrue.age)

##joining both datasets (due to the scale, it is better a minus than a division)
sp.comp <- left_join(sp.max.1, sp.min.1, by = "root.age") %>% 
           mutate(est.ratio = maxEst-minEst,
                  true.ratio = maxTrue-minTrue) %>% 
           pivot_longer(cols = c("est.ratio", "true.ratio"), names_to = "age",
                        values_to = "ratio")

##saving for merging two figures (bifurcating and budding)
write_csv(sp.comp, file = file.path(getwd(), pro, q_error, 
                                    "sp.comp.bud.csv"))

####defining how many trees are qualitatively wrong 
q.wrong <- sp.comp %>% group_by(extinction.x) %>% 
  filter(ratio < 0) %>% count() %>% 
  mutate(q_wrong = n/10)

##figure ratio of oldest and youngest species for the estimated age
png("text/figures/Figure5.old.vs.young.bud.png",
    width = 15, height = 15, units = "cm", 
    pointsize = 8, res = 300)

old.vs.young.bud <-ggplot(sp.comp, aes(x = ratio, fill = age))+ 
               geom_histogram(alpha = 0.5, position = 'identity',
                       breaks = seq(-0.7, 1, by = 0.02))+
              geom_vline(xintercept = 0, color = "red")+
              geom_text(data = q.wrong, aes(x = -0.35, y = 105,
                          label = paste("Qualitatively wrong:", q_wrong,"%")),
                        inherit.aes = FALSE,
                        size = 3, family = "serif")+
              facet_wrap(~ extinction.x, nrow = 3,
              labeller = as_labeller(c(High = "High",
                                      Intermediate = "Intermediate",
                                      no_ext = "No extinction")))+
              scale_fill_discrete(name = NULL,
                     breaks=c("est.ratio", "true.ratio"),
                      labels=c("Phylogenetic ages",
                               "True ages"))+
              xlab("Oldest sp - Youngest sp")+
              ylab("Tree Number")+
              ggtitle("Budding speciation")+
              theme_bw()

old.vs.young + mynamestheme

dev.off()

#############distribution between two random species############################
sp.random.bud <- df.ext.levels %>% select(species, mode,
                                   tree, root.age, rTrue.age, rEstimated.age,
                                   extinction) %>% 
                                group_by(tree, extinction) %>% sample_n(2) %>%
                                  arrange(rEstimated.age,
                                          .by_group = TRUE) %>% 
                                summarise(est.ratio = diff(rEstimated.age),
                                          true.ratio = diff(rTrue.age)) %>% 
          pivot_longer(cols = c("est.ratio", "true.ratio"), names_to = "age",
                       values_to = "ratio")
               
####defining how many trees are qualitatively wrong 
q.wrong.random.bud <- sp.random.bud %>% group_by(extinction) %>% 
  filter(ratio < 0) %>% count() %>% 
  mutate(q_wrong = n/10)                               


##figure ratio of an old and younger random species for the estimated age
png("Figure5.old.vs.young.bud.random.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

old.vs.young.random.bud <-ggplot(sp.random.bud, aes(x = ratio, fill = age))+ 
  geom_histogram(alpha = 0.5, position = 'identity',
                 breaks = seq(-0.7, 1, by = 0.02))+
  geom_vline(xintercept = 0, color = "red")+
  geom_text(data = q.wrong.random, aes(x = -0.38, y = 300,
                                label = paste("Qualitatively wrong:",
                                              q_wrong,"%")),
                inherit.aes = FALSE,
                size = 3, family = "serif")+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(High = "High",
                                      Intermediate = "Intermediate",
                                      no_ext = "No extinction")))+
  scale_fill_discrete(name = NULL,
                      breaks=c("est.ratio", "true.ratio"),
                      labels=c("Phylogenetic ages",
                               "True ages"))+
  xlab("Old sp - Young sp (random)")+
  ylab("Tree Number")+
  ggtitle("Budding speciation")+
  theme_bw()

old.vs.young.random + mynamestheme

dev.off()

