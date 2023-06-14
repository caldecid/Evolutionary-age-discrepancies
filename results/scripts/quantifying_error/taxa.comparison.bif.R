##############################################################################
#################      Taxa simulations bifurcating     #####################
#################             for comparisons           #####################
#############################################################################


##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

# mu = 0.9 ----------------------------------------------------------------

load(file.path(getwd(), pro, q_error, "Tree.mu.09.RData"))

##simulating species taxonomy; only bifurcating speciation. Trees mu = 0.9
taxa.mu.09.bif <- lapply(Tree.mu.09, sim.taxonomy, beta = 1, lambda.a = 0)

save(taxa.mu.09.bif,
       file = file.path(getwd(), pro, q_error, "taxa.mu.09.bif.RData"))

##Calculating estimated and true ages
taxa.extant.09.bif <- Map(getAgesExtantSpecies, taxa.mu.09.bif, Tree.mu.09)

##colapsing list into a dataframe
clado.df.09.bif <- do.call("rbind", taxa.extant.09.bif)

##root ages
root.ages.09 <- sapply(Tree.mu.09, tree.max)

clado.df.09.bif$root.age <- rep(root.ages.09, each = 100)

##extinction rates
clado.df.09.bif$mu <- rep(0.9, nrow(clado.df.09.bif))

##tree
clado.df.09.bif$tree <- rep(paste0("tree.09.", seq(1, 1000, 1)), each = 100)

##budding speciation = 0 ; bifurcating = 1
clado.df.09.bif$mode <- rep("Bifurcating", nrow(clado.df.09.bif))


##manipulating the data frame
clado.df.bif.09 <- clado.df.09.bif %>% select(label, Age, tip_length_reconstr, mode,
                                          mu, root.age, tree) %>% 
  rename(species = label,
         True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  mutate(rTrue.age = True.age/root.age,
         rEstimated.age = Estimated.age/root.age,
         accuracy = rEstimated.age - rTrue.age)

##defining the intensity of extinction
clado.df.bif.09$extinction <- rep("High", nrow(clado.df.bif.09))

write_csv(clado.df.bif.09,
          file = file.path(getwd(), pro, q_error, "taxa.mu.09.bif.csv"))

# mu = 0.5 ----------------------------------------------------------------

load(file.path(getwd(), pro, q_error, "Tree.mu.05.RData"))

##simulating species taxonomy; only bifurcating speciation. Trees mu = 0.5

taxa.mu.05.bif <- lapply(Tree.mu.05, sim.taxonomy, beta = 1, lambda.a = 0)

save(taxa.mu.05.bif, file = file.path(getwd(),
                                      pro, q_error, "taxa.mu.05.bif.RData"))

##Calculating estimated and true ages
taxa.extant.05.bif <- Map(getAgesExtantSpecies, taxa.mu.05.bif, Tree.mu.05)

##colapsing list into a dataframe
clado.df.05.bif <- do.call("rbind", taxa.extant.05.bif)

##speciation rates
clado.df.05.bif$lambda <- rep(1, nrow(clado.df.05.bif))

##extinction rates
clado.df.05.bif$mu <- rep(0.5, nrow(clado.df.05.bif))

##root ages
root.ages.05 <- sapply(Tree.mu.05, tree.max)

clado.df.05.bif$root.age <- rep(root.ages.05, each = 100)


##budding speciation = 0 ; bifurcating = 1
clado.df.05.bif$mode <- rep("Bifurcating", nrow(clado.df.05.bif))

##tree
clado.df.05.bif$tree <- rep(paste0("tree.05.", seq(1, 1000, 1)), each = 100)

##manipulating the data frame
clado.df.bif.05 <- clado.df.05.bif %>% select(label, Age, tip_length_reconstr,
                                            mode, tree, mu, root.age) %>% 
  rename(species = label,
         True.age = Age,
         Estimated.age = tip_length_reconstr) %>%
  mutate(rTrue.age = True.age/root.age,
         rEstimated.age = Estimated.age/root.age,
         accuracy = rEstimated.age - rTrue.age)

##intensity of extinction
clado.df.bif.05$extinction <- rep("Intermediate", nrow(clado.df.bif.05))

##saving
write_csv(clado.df.bif.05, 
          file = file.path(getwd(), pro, q_error, "taxa.mu.05.bif.csv"))

# mu = 0 ----------------------------------------------------------------
load(file.path(getwd(), pro, q_error, "Tree.mu.0.RData"))

##simulating species taxonomy; only bifurcating speciation. Trees mu = 0
taxa.mu.0.bif <- lapply(Tree.mu.0, sim.taxonomy, beta = 1, lambda.a = 0)

##saving taxonomy list
save(taxa.mu.0.bif, file = file.path(getwd(), pro, q_error,
                                     "taxa.mu.0.bif.RData"))

##Calculating estimated and true ages
taxa.extant.0.bif <- Map(getAgesExtantSpecies, taxa.mu.0.bif, Tree.mu.0)

##colapsing list into a dataframe
clado.df.0.bif <- do.call("rbind", taxa.extant.0.bif)

##extinction rates
clado.df.0.bif$mu <- rep(0, nrow(clado.df.0.bif))

##root ages
root.ages.0 <- sapply(Tree.mu.0, tree.max)

clado.df.0.bif$root.age <- rep(root.ages.0, each = 100)

##tree
clado.df.0.bif$tree <- rep(paste0("tree.0.", seq(1, 1000, 1)), each = 100)

##budding speciation = 0 ; bifurcating = 1
clado.df.0.bif$mode <- rep("Bifurcating", nrow(clado.df.0.bif))


##manipulating the data frame
clado.df.bif.0 <- clado.df.0.bif %>% select(label, Age, tip_length_reconstr,
                                      mode, tree ,mu, root.age) %>% 
                              rename(species = label,
                                     True.age = Age,
                                     Estimated.age = tip_length_reconstr) %>%
                              mutate(rTrue.age = True.age/root.age,
                                     rEstimated.age = Estimated.age/root.age,
                                     accuracy = rEstimated.age - rTrue.age)
##intensity of extinction
clado.df.bif.0$extinction <- rep("no_ext", nrow(clado.df.bif.0))

write_csv(clado.df.bif.0,
          file = file.path(getwd(), pro, q_error, "taxa.mu.0.bif.csv"))


# Merging data frames -----------------------------------------------------

df.ext.levels.bif <- rbind(clado.df.bif.09, clado.df.bif.05, clado.df.bif.0)

###needed transformations
df.ext.levels.bif$extinction <- as.factor(df.ext.levels.bif$extinction)

df.ext.levels.bif$extinction <- factor(df.ext.levels.bif$extinction, 
                                   levels = c("no_ext", "Intermediate", "High"))

##writing
write_csv(df.ext.levels.bif, file = file.path(getwd(), pro, q_error,
                                          "taxa.ext.lev.bif.csv"))
##reading
df.ext.levels.bif <- read_csv(file = file.path(getwd(), pro, q_error,
                                                "taxa.ext.lev.bif.csv")) 
##defining levels for facet
df.ext.levels.bif$extinction <- as.factor(df.ext.levels.bif$extinction)

df.ext.levels.bif$extinction <- factor(df.ext.levels$extinction, 
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
age.ext.bif <- df.ext.levels.bif %>% 
                          pivot_longer(cols = c("rTrue.age", "rEstimated.age"),
                                          names_to = "type", values_to = "age")
##type as factor
age.ext.bif$type <- factor(age.ext.bif$type, levels = c("rTrue.age", "rEstimated.age"),
                       labels = c("rTrue.age", "rEstimated.age"))

##png object
png("text/supplementary/FigureSM4.Dist.true.phy.bif.png",
    width = 14, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f3 <- age.ext.bif %>% 
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
  #scale_y_log10()+
  xlab("Ages (scaled)")+
  ylab("Species Count")+
  ggtitle("Bifurcating speciation")



f3 + mynamestheme

dev.off()

####histograms for accuracy (rTrue.Age - rEstimated.age)

##determining how many species are 100% accurate
accurate.age.bif <- df.ext.levels.bif %>%
  mutate(abs.accuracy = abs(accuracy)) %>% 
  filter(abs.accuracy < 1e-15) %>% group_by(extinction) %>%
  count()


######figure displaying a line and number of the accurate species
##png object
png("text/supplementary/FigureSM3.ratio.age.bif.png",
    width = 12, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f4 <- df.ext.levels.bif %>%
  mutate(abs.accuracy = abs(accuracy)) %>% 
  filter(abs.accuracy > 1e-15) %>% 
  ggplot(aes(x = accuracy))+
  geom_histogram(position = "dodge", breaks = seq(-0.1, 0.4, by = 0.02))+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = n),
               data = accurate.age.bif,
               colour = "red")+
  geom_point(aes(x = 0, y = n), data = accurate.age.bif, colour = "red")+
  geom_text(data = accurate.age.bif, aes(x = 0.2, y = 50000,
                   label = paste("Correct estimation:", round(n/1000),"%")),
                   size = 3, family = "serif")+
  xlim(-0.4, 0.4)+
  #ylim(0, 5000)+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(High = "High extinction",
                                      Intermediate = "Intermediate extinction",
                                      no_ext = "No extinction")))+
  theme_bw()+
  xlab("Phylogenetic age - True age (scaled)")+
  ylab("Species count")+
  ggtitle("Bifurcating speciation")
  

f4 + mynamestheme

dev.off()


##distribution of ratio between youngest and oldest sp for rphyletic, selecting 
#species and calculate the ratio with the true ages
sp.max.bif <- df.ext.levels.bif %>% select(species, root.age, rTrue.age, rEstimated.age,
                                   extinction) %>% 
  group_by(root.age, extinction) %>% 
  filter(rEstimated.age == max(rEstimated.age))     
##eliminate duplicated rows (there are at least two max)
sp.max.1.bif <- sp.max.bif[-which(duplicated(sp.max.bif$root.age)),] %>% 
                        rename(maxEst = rEstimated.age,
                               maxTrue = rTrue.age)

##defining the sp with the minimun ages
sp.min.bif <- df.ext.levels.bif %>% select(species, root.age, rTrue.age, rEstimated.age,
                                   extinction) %>% 
  group_by(root.age, extinction) %>% 
  filter(rEstimated.age == min(rEstimated.age))    

##eliminate the duplicated rows (there are at least two min)
sp.min.1.bif <- sp.min.bif[-which(duplicated(sp.min.bif$root.age)),] %>% 
                  rename(sp = species, 
                         minEst = rEstimated.age,
                         minTrue = rTrue.age)

##joining both datasets (due to the scale, it is better a minus than a division)
sp.comp.bif <- left_join(sp.max.1.bif, sp.min.1.bif, by = "root.age") %>% 
  mutate(est.ratio = maxEst-minEst,
         true.ratio = maxTrue-minTrue) %>% 
  pivot_longer(cols = c("est.ratio", "true.ratio"), names_to = "age",
               values_to = "ratio")

####defining how many sp are qualitatively wrong 
q.wrong.bif <- sp.comp.bif %>% group_by(extinction.x) %>% 
  filter(ratio < 0) %>% count() %>% 
  mutate(q_wrong = n/10)

##figure ratio of oldest and youngest species for the estimated age
png("text/figures/Figure5.oldest.vs.young.png", 
    width = 12, height = 15, units = "cm", 
    pointsize = 8, res = 300)

old.vs.young.bif <- sp.comp.bif %>% filter(age == "true.ratio") %>% 
                    ggplot(aes(x = ratio, fill = age))+ 
                    geom_histogram(alpha = 0.5, position = 'identity',
                                   breaks = seq(-0.7, 1, by = 0.02),
                                   fill = "blue")+
                    geom_vline(xintercept = 0, color = "red")+
                    geom_text(aes(x = -0.36, y = 150,
                               label = paste("Wrong:", "0","%")),
                              inherit.aes = FALSE,
                              size = 3, family = "serif")+
                    facet_wrap(~ extinction.x, nrow = 3,
                              labeller = as_labeller(c(High = "High extinction",
                                      Intermediate = "Intermediate extinction",
                                                    no_ext = "No extinction")))+
                    scale_fill_discrete(name = NULL,
                                        breaks=c("est.ratio", "true.ratio"),
                                        labels=c("Phylogenetic ages",
                                                 "True ages"))+
                    xlab(NULL)+
                    ylab("Frequency")+
                    ggtitle("Bifurcating speciation")+
                    theme_bw()+
                    mynamestheme+
                    theme(legend.position = "none")



dev.off()

################merging bifurcating and budding in one figure################

###calling budding data
sp.comp.bud <- read_csv(file = file.path(getwd(), pro, q_error, 
                                         "sp.comp.bud.csv"))

##defining levels for facets
sp.comp.bud$extinction.x <- as.factor(sp.comp.bud$extinction.x)

sp.comp.bud$extinction.x <- factor(sp.comp.bud$extinction.x, 
                             levels = c("no_ext", "Intermediate", "High"))
  

####defining how many sp are qualitatively wrong 
q.wrong.bud<- sp.comp.bud %>% group_by(extinction.x) %>% 
  filter(ratio < 0) %>% count() %>% 
  mutate(q_wrong = n/10)

##budding figure

old.vs.young.bud <- sp.comp.bud %>% filter(age == "true.ratio") %>% 
                    ggplot(aes(x = ratio, fill = age))+ 
                    geom_histogram(alpha = 0.5, position = 'identity',
                                   breaks = seq(-0.7, 1, by = 0.02),
                                   fill = "blue")+
                    geom_vline(xintercept = 0, color = "red")+
                    geom_text(data = q.wrong.bud, aes(x = -0.4, y = 90,
                                        label = paste("Wrong:",
                                                q_wrong,"%")),
                              inherit.aes = FALSE,
                              size = 3, family = "serif")+
                    facet_wrap(~ extinction.x, nrow = 3,
                            labeller = as_labeller(c(High = "High extinction",
                                      Intermediate = "Intermediate extinction",
                                                    no_ext = "No extinction")))+
                    scale_fill_discrete(name = NULL,
                                        breaks=c("est.ratio", "true.ratio"),
                                        labels=c("Phylogenetic ages",
                                                 "True ages"))+
                    xlab(NULL)+
                    ylab(NULL)+
                    ggtitle("Budding speciation")+
                    theme_bw()+
                    mynamestheme+
                    theme(legend.position = "none")
  

######Both figures

png("text/figures/Figure5.oldest.vs.youngest.png", 
    width = 15, height = 15, units = "cm", 
    pointsize = 8, res = 300)


plot <- ggarrange(old.vs.young.bif, old.vs.young.bud, ncol = 2,
          common.legend = FALSE)

annotate_figure(plot, bottom = text_grob("Delta True Age", 
                       family = "serif", face = "bold", size = 13))

dev.off()



#############distribution between two random species############################
sp.random.bif <- df.ext.levels.bif %>% select(species, mode,
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
q.wrong.random.bif <- sp.random.bif %>% group_by(extinction) %>% 
  summarise(q_wrong = length(ratio[ratio <0])/10)                      


##Bifurcating figure
old.vs.young.random.bif <- sp.random.bif %>% filter(age == "true.ratio") %>% 
          ggplot(aes(x = ratio, fill = age))+ 
          geom_histogram(alpha = 0.5, position = 'identity',
                         breaks = seq(-0.7, 1, by = 0.02), fill = "blue")+
          geom_vline(xintercept = 0, color = "red")+
          geom_text(data = q.wrong.random.bif, aes(x = -0.36, y = 400,
                                       label = paste("Wrong:",
                                                             q_wrong,"%")),
                    inherit.aes = FALSE,
                    size = 3, family = "serif")+
          facet_wrap(~ extinction, nrow = 3,
                     labeller = as_labeller(c(High = "High extinction",
                                      Intermediate = "Intermediate extinction",
                                              no_ext = "No extinction")))+
          scale_fill_discrete(name = NULL,
                              breaks=c("est.ratio", "true.ratio"),
                              labels=c("Phylogenetic ages",
                                       "True ages"))+
          xlab(NULL)+
          ylab("Frequency")+
          ggtitle("Bifurcating speciation")+
          theme_bw()+
          mynamestheme
  

##########calling budding figure
##reading
df.ext.levels <- read_csv(file = file.path(getwd(), pro, q_error,
                                           "taxa.ext.lev.csv"))

df.ext.levels$mode <- rep("Budding", nrow(df.ext.levels))

##needed transformations
df.ext.levels$extinction <- as.factor(df.ext.levels$extinction)

df.ext.levels$extinction <- factor(df.ext.levels$extinction, 
                                   levels = c("no_ext", "Intermediate", "High"))

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


old.vs.young.random.bud <- sp.random.bud %>% filter(age == "true.ratio") %>% 
                 ggplot(aes(x = ratio, fill = age))+ 
              geom_histogram(alpha = 0.5, position = 'identity',
                             breaks = seq(-0.7, 1, by = 0.02),
                             fill = "blue")+
              geom_vline(xintercept = 0, color = "red")+
              geom_text(data = q.wrong.random.bud, aes(x = -0.36, y = 200,
                                                   label = paste("Wrong:",
                                                                 q_wrong,"%")),
                        inherit.aes = FALSE,
                        size = 3, family = "serif")+
              facet_wrap(~ extinction, nrow = 3,
                         labeller = as_labeller(c(High = "High extinction",
                                      Intermediate = "Intermediate extinction",
                                        no_ext = "No extinction")))+
              scale_fill_discrete(name = NULL,
                                  breaks=c("est.ratio", "true.ratio"),
                                  labels=c("Phylogenetic ages",
                                           "True ages"))+
              xlab(NULL)+
              ylab(NULL)+
              ggtitle("Budding speciation")+
              theme_bw()+
              mynamestheme

##PNG object

##figure ratio of an old and younger random species for the estimated age
png("text/figures/Figure6.old.vs.young.random.png",
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)


plot.random <- ggarrange(old.vs.young.random.bif, old.vs.young.random.bud, 
          ncol = 2, common.legend = FALSE)

annotate_figure(plot.random, bottom = text_grob("Delta True Age", 
                                  family = "serif", face = "bold", size = 13))

dev.off()


