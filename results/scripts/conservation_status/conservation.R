##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

# Simulations for testing the conservation status  ------------------------

#####################simulating trees####################################

###high extinction scenario
trees.high <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.3, mu = 0.25)

##intermediate extinction scenario
trees.int <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.3, mu = 0.15)

##low extinction scenario
trees.low <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.3, mu = 0.05)


###################prunning extinct species ###########################

####high extinction
##trees
trees.extant.high <- lapply(trees.high, prune.fossil.tips) 
##root ages
root.high <- sapply(trees.high, tree.max)

##intermediate extinction
##trees
trees.extant.int <- lapply(trees.int, prune.fossil.tips) 
##root ages
root.int <- sapply(trees.int, tree.max)

##low extinction
##trees
trees.extant.low <- lapply(trees.low, prune.fossil.tips) 
##root ages
root.low <- sapply(trees.low, tree.max)

###saving trees
save(trees.high, trees.int, trees.low,
      file = file.path(getwd(), pro, c_status,"trees.conservation.RData"))

#############simulating taxonomies only budding######################

##high extinction
tax.high <- lapply(trees.high, sim.taxonomy, beta = 0)

##intermediate extinction
tax.int <- lapply(trees.int, sim.taxonomy, beta = 0)

##low extinction
tax.low <- lapply(trees.low, sim.taxonomy, beta = 0)

#############get ages and variables####################

######high extinction
ages.high <- Map(getAgesExtantSpecies, tax.high, trees.high) 
                                              
ages.high <- do.call("rbind", ages.high) %>% 
                    rename(True.age = Age,
                      Estimated.age = tip_length_reconstr) %>%
                          select(label, True.age, Estimated.age)


##root age
ages.high$root.age <- rep(root.high, each = 100)

##tree numbers
ages.high$tree <- rep(paste0("tree.", 1:100), each = 100)

##scale ages
ages.high <- ages.high %>% mutate(rTrue.age = True.age/root.age,
                                  rPhylo.age = Estimated.age/root.age)

##nesting by tree
ages.high.nest <- ages.high %>% group_by(tree) %>% nest()

##defining intensity of extinction in the tree
ages.high.nest$extinction <- rep("high", 100)

############intermediate extinction
ages.int <- Map(getAgesExtantSpecies, tax.int, trees.int) 

ages.int <- do.call("rbind", ages.int) %>% 
             rename(True.age = Age,
                    Estimated.age = tip_length_reconstr) %>%
             select(label, True.age, Estimated.age)

##root age
ages.int$root.age <- rep(root.int, each = 100)

##tree numbers
ages.int$tree <- rep(paste0("tree.", 1:100), each = 100)

##scale ages
ages.int <- ages.int %>% mutate(rTrue.age = True.age/root.age,
                                rPhylo.age = Estimated.age/root.age)

##nesting by tree
ages.int.nest <- ages.int %>% group_by(tree) %>% nest()

##defining intensity of extinction in the tree
ages.int.nest$extinction <- rep("intermediate", 100)

############low extinction
ages.low <- Map(getAgesExtantSpecies, tax.low, trees.low) 

ages.low <- do.call("rbind", ages.low) %>% 
                    rename(True.age = Age,
                    Estimated.age = tip_length_reconstr) %>%
                    select(label, True.age, Estimated.age)

##root ages
ages.low$root.age <- rep(root.low, each = 100)

##tree numbers
ages.low$tree <- rep(paste0("tree.", 1:100), each = 100)

##scale ages
ages.low <- ages.low %>% mutate(rTrue.age = True.age/root.age,
                                rPhylo.age = Estimated.age/root.age)

##nesting by tree
ages.low.nest <- ages.low %>% group_by(tree) %>% nest()

##defining intensity of extinction in the tree
ages.low.nest$extinction <- rep("low", 100)


# simulating the extinction signal -----------------------------------------

####parameters for assigning the extinction probability to each tree
##slope
b = rep(c(1,2,3), c(33,33,34))
##size of the effect
b.effect = rep(c("low", "intermediate", "high"), c(33,33,34))
##mean error
mean.error = seq(0.5, 4, length.out = 20)
##standard deviation of the error
sd.error = seq(0.5, 2.5, length.out =20)

#####################extinction signal 

#############high extinction

##list for keeping the extinction status
status.high <- vector("list", 100)

##assuming a normal error and sampling the error parameters
for(i in seq_along(status.high)){
  status.high[[i]] <- b[i] * ages.high.nest[[2]][[i]]$rTrue.age
                                + rnorm(nrow(ages.high.nest[[2]][[i]]),
                        mean = sample(mean.error, 1), sd = sample(sd.error, 1)) 
}

##assigning status
for(i in seq_along(status.high)){
  ages.high.nest[[2]][[i]]$status <- status.high[[i]]
}

#############intermediate extinction

##list for keeping the extinction status
status.int <- vector("list", 100)

##assuming a normal error and sampling the error parameters
for(i in seq_along(status.int)){
  status.int[[i]] <- b[i] * ages.int.nest[[2]][[i]]$rTrue.age
  + rnorm(nrow(ages.int.nest[[2]][[i]]),
          mean = sample(mean.error, 1), sd = sample(sd.error, 1)) 
}

##assigning status
for(i in seq_along(status.int)){
  ages.int.nest[[2]][[i]]$status <- status.int[[i]]
}

#############low extinction

##list for keeping the extinction status
status.low <- vector("list", 100)

##assuming a normal error and sampling the error parameters
for(i in seq_along(status.low)){
  status.low[[i]] <- b[i] * ages.low.nest[[2]][[i]]$rTrue.age
  + rnorm(nrow(ages.low.nest[[2]][[i]]),
          mean = sample(mean.error, 1), sd = sample(sd.error, 1)) 
}

##assigning status
for(i in seq_along(status.int)){
  ages.low.nest[[2]][[i]]$status <- status.low[[i]]
}



##saving for further work
ages.nest <- rbind(ages.high.nest, ages.int.nest, ages.low.nest)

##effect size
ages.nest$effect_size <- rep(b.effect, 3)

save(ages.nest, file = file.path(getwd(), pro, c_status, "ages.nest.RData"))

##discretizing

for(i in 1:nrow(ages.nest)) {
  ages.nest[[2]][[i]]$status.dis <- discretize(ages.nest[[2]][[i]]$status, 
                                       method = "fixed", 
                                       breaks = c(0, 0.2, 0.4, 0.6, 0.8, Inf),
                                      labels = c("LC", "NT", "VU", "EN", "CR"))
}


##unnesting 

ages.unnest <- ages.nest %>% unnest(data)

ages.unnest$extinction <- as.factor(ages.unnest$extinction)

ages.unnest$effect_size <- as.factor(ages.unnest$effect_size)

##ordering factors
ages.unnest$effect_size <- factor(ages.unnest$effect_size,
                                  levels=c("low", "intermediate", "high"))


ages.unnest$extinction <- factor(ages.unnest$extinction,
                                 levels=c("low", "intermediate", "high"))

##saving
write_csv(ages.unnest, file = file.path(getwd(), pro, c_status,
                                        "ages.unnest.csv"))

ages.unnest <- read_csv(file = file.path(getwd(), pro, c_status,
                                         "ages.unnest.csv"))


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

########################Effect size
png("text/figures/Figure7.conservation_effect_size.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.unnest, aes(x = status.dis, y = rTrue.age))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~effect_size, labeller = as_labeller(c(high = "High",
                                                intermediate = "Intermediate",
                                                low = "Low")))+
  ylab("True age")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme
  

p2 <- ggplot(ages.unnest, aes(x = status.dis, y = rPhylo.age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~effect_size, labeller = as_labeller(c(high = "High",
                                                intermediate = "Intermediate",
                                                  low = "Low")))+
  ylab("Phylogenetic age")+
  xlab("Conservation status")+
  theme_bw()+
  ylim(0,1)+
  mynamestheme

##setting top
top <- text_grob("Effect size", family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


########################extinction background
png("text/figures/Figure8.conservation_extinction.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.unnest, aes(x = status.dis, y = rTrue.age))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low",
                                                    intermediate = "Intermediate",
                                                    high = "High")))+
  ylab("True age")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme#+
  #coord_flip()


p2 <- ggplot(ages.unnest, aes(x = status.dis, y = rPhylo.age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low",
                                                intermediate = "Intermediate",
                                                high = "High")))+
  ylab("Phylogenetic age")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Extinction Background",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


###########ridgelines

p.ridge.true <- ggplot(ages.unnest, aes(x = rTrue.age, y = status.dis, 
                                        fill = status.dis))+
                geom_density_ridges() +
                theme_ridges() + 
                theme(legend.position = "none")+
                facet_wrap(~extinction, ncol = 3,
                           labeller = as_labeller(c(low = "Low",
                                     intermediate = "Intermediate",
                                                   high = "High")))+
                theme_bw()+
                ylab("Conservation status")+
                xlab("True Age (scale)")+
                mynamestheme


p.ridge.phylo <- ggplot(ages.unnest, aes(x = rPhylo.age, y = status.dis, 
                                         fill = status.dis))+
                 geom_density_ridges() +
                 theme_ridges() + 
                 theme(legend.position = "none")+
                 facet_wrap(~extinction, ncol = 3,
                  labeller = as_labeller(c(low = "Low",
                                      intermediate = "Intermediate",
                                      high = "High")))+
                 theme_bw()+
                 ylab(NULL)+
                xlab("Phylogenetic Age (scale)")+
                ylab("Conservation status")+
                mynamestheme

##setting top
top <- text_grob("Extinction Background",
                 family = "serif", size = 13,  face = "bold")

left <- text_grob("Conservation status",
                  family = "serif", size = 11,  face = "bold")

grid.arrange(p.ridge.true, p.ridge.phylo, nrow = 2, top = top)                


# stochastic function -----------------------------------------------------


##future package conditions 
plan("multisession", workers = 3)

###high extinction
ages.high.stochastic <- vector("list", length = length(trees.extant.high))


for(i in seq_along(trees.extant.high)) {
    ages.high.stochastic[[i]] <- try(stochastic.sp.age(trees.extant.high[[i]]))
} 


  
##intermediate extinction
ages.int.stochastic <- vector("list", length = length(trees.extant.int))

for(i in seq_along(trees.extant.int)) {
  ages.int.stochastic[[i]] <- try(stochastic.sp.age(trees.extant.int[[i]]))
}


##low extinction
ages.low.stochastic <- vector("list", length = length(trees.extant.low))


for(i in seq_along(trees.extant.high)) {
  ages.low.stochastic[[i]] <- try(stochastic.sp.age(trees.extant.low[[i]]))
}

save(ages.high.stochastic, ages.int.stochastic, ages.low.stochastic, 
     file = file.path(getwd(), pro, c_status, "ages.stochastic.nest.RData"))

#####organizing results############

##removing the df with problems

##high extinction
ages.high.stochastic.1 <- ages.high.stochastic[lapply(ages.high.stochastic,
                                                      class) != "try-error"]



##intermediate extinction
ages.int.stochastic.1 <- ages.int.stochastic[lapply(ages.int.stochastic,
                                                    class) != "try-error"]

##low extinction
ages.low.stochastic.1 <- ages.low.stochastic[lapply(ages.low.stochastic,
                                                    class) != "try-error"]


######removing trees which df are causing problems

##high extinction
trees.high.1 <- trees.high[lapply(ages.high.stochastic,
                                                 class) != "try-error"]

##root age
root.high.1 <- sapply(trees.high.1, tree.max)

##assigning tree number and root age
for(i in seq_along(ages.high.stochastic.1)){
  ages.high.stochastic.1[[i]]$tree <- rep(paste0("tree.", i),
                                          nrow(ages.high.stochastic.1[[i]]))
  ages.high.stochastic.1[[i]]$root.age <- rep(root.high.1[i],
                                              nrow(ages.high.stochastic.1[[i]]))
  ages.high.stochastic.1[[i]]$status <- status.high.1[[i]]
  
  ages.high.stochastic.1[[i]]$effect_size <- b.effect.high[i]
}

##unlist
ages.high.st <- do.call("rbind", ages.high.stochastic.1) 

ages.high.st$species <- paste0(ages.high.st$species,".", ages.high.st$tree)

ages.high.st$extinction <- rep("high", nrow(ages.high.st))

ages.high.st <- ages.high.st %>% unnest(dist_ages) %>% 
  filter(Probability == max(Probability)) %>% 
  rename(max.prob.age = Age)

ages.high.st$max.prob.age <- as.numeric(as.character(ages.high.st$max.prob.age))

ages.high.st <- ages.high.st %>% mutate(rStochastic = max.prob.age/root.age)

######intermediate extinction
trees.int.1 <- trees.int[lapply(ages.int.stochastic,
                                               class) != "try-error"]
##root age
root.int.1 <- sapply(trees.int.1, tree.max)

##assigning tree number and root age
for(i in seq_along(ages.int.stochastic.1)){
  ages.int.stochastic.1[[i]]$tree <- rep(paste0("tree.", i),
                                          nrow(ages.int.stochastic.1[[i]]))
  ages.int.stochastic.1[[i]]$root.age <- rep(root.int.1[i],
                                              nrow(ages.int.stochastic.1[[i]]))
  ages.int.stochastic.1[[i]]$status <- status.int.1[[i]]
  
  ages.int.stochastic.1[[i]]$effect_size <- b.effect.int[i]
}

##unlist
ages.int.st <- do.call("rbind", ages.int.stochastic.1) 

ages.int.st$species <- paste0(ages.int.st$species,".", ages.int.st$tree)

ages.int.st$extinction <- rep("intermediate", nrow(ages.int.st))

ages.int.st <- ages.int.st %>% unnest(dist_ages) %>% 
  filter(Probability == max(Probability)) %>% 
  rename(max.prob.age = Age)

ages.int.st$max.prob.age <- as.numeric(as.character(ages.int.st$max.prob.age))

ages.int.st <- ages.int.st %>% mutate(rStochastic = max.prob.age/root.age)

##low extinction
trees.low.1 <- trees.low[lapply(ages.low.stochastic,
                                               class) != "try-error"]
##root age
root.low.1 <- sapply(trees.low.1, tree.max)

##assigning tree number and root age
for(i in seq_along(ages.low.stochastic.1)){
  ages.low.stochastic.1[[i]]$tree <- rep(paste0("tree.", i),
                                          nrow(ages.low.stochastic.1[[i]]))
  ages.low.stochastic.1[[i]]$root.age <- rep(root.high.1[i],
                                              nrow(ages.low.stochastic.1[[i]]))
  
  ages.low.stochastic.1[[i]]$status <- status.low.1[[i]]
  
  ages.low.stochastic.1[[i]]$effect_size <- b.effect.low[i]
  
}


##unlist
ages.low.st <- do.call("rbind", ages.low.stochastic.1) 

ages.low.st$species <- paste0(ages.low.st$species,".", ages.low.st$tree)

ages.low.st$extinction <- rep("low", nrow(ages.low.st))

ages.low.st <- ages.low.st %>% unnest(dist_ages) %>% 
                filter(Probability == max(Probability)) %>% 
                rename(max.prob.age = Age) 

ages.low.st$max.prob.age <- as.numeric(as.character(ages.low.st$max.prob.age))

ages.low.st <- ages.low.st %>% mutate(rStochastic = max.prob.age/root.age)

##############assigning status

##high extinction
status.high.1 <- status.high[lapply(ages.high.stochastic,
                                    class) != "try-error"]


b.effect.high <- b.effect[lapply(ages.high.stochastic,
                                 class) != "try-error"]
##intermediate extinction
status.int.1 <- status.int[lapply(ages.int.stochastic,
                                  class) != "try-error"]

b.effect.int <- b.effect[lapply(ages.int.stochastic,
                                class) != "try-error"]
##low extinction
status.low.1 <- status.low[lapply(ages.low.stochastic,
                                  class) != "try-error"]

b.effect.low <- b.effect[lapply(ages.low.stochastic,
                                class) != "try-error"]

ages.stochastic <- rbind(ages.high.st, ages.int.st, ages.low.st)

write_csv(ages.stochastic, file = file.path(getwd(), pro, c_status,
                                            "ages.stochastic.csv"))

save(ages.high.stochastic.1, ages.int.stochastic.1, ages.low.stochastic.1, 
     file = file.path(getwd(), pro, c_status,
                      "ages.stochastic.list.RData"))

#####Discretizing the status

ages.stochastic$status.dis <- discretize(ages.stochastic$status, 
                                      method = "fixed", 
                                      breaks = c(0, 0.2, 0.4, 0.6, 0.8, Inf),
                                      labels = c("LC", "NT", "VU", "EN", "CR"))

ages.stochastic$extinction <- as.factor(ages.stochastic$extinction)

ages.stochastic$effect_size <- as.factor(ages.stochastic$effect_size)

##ordering factors
ages.stochastic$effect_size <- factor(ages.stochastic$effect_size,
                                  levels=c("low", "intermediate", "high"))


ages.stochastic$extinction <- factor(ages.stochastic$extinction,
                                 levels=c("low", "int", "high"))


# Figures stochastic ages -------------------------------------------------

png("text/figures/Figure9.st.conservation_effect_size.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

##true age
p1 <- ggplot(ages.unnest, aes(x = status.dis, y = rTrue.age))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~effect_size, labeller = as_labeller(c(high = "High",
                                                  intermediate = "Intermediate",
                                                  low = "Low")))+
  ylab("True age")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme

##stochastic age
p2 <- ggplot(ages.stochastic, aes(x = status.dis, y = rStochastic))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~effect_size, labeller = as_labeller(c(high = "High",
                                                  intermediate = "Intermediate",
                                                  low = "Low")))+
  ylab("Stochastic highest probability age")+
  xlab(NULL)+
  ylim(0,1)+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Effect size",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


###########extinction background

png("text/figures/Figure10.st.conservation_extinction.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

##true age
p1 <- ggplot(ages.unnest, aes(x = status.dis, y = rTrue.age))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High",
                                                    intermediate = "Intermediate",
                                                    low = "Low")))+
  ylab("True age")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme

##stochastic age
p2 <- ggplot(ages.stochastic, aes(x = status.dis, y = rStochastic))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High",
                                                    int = "Intermediate",
                                                    low = "Low")))+
  ylab("Stochastic highest probability age")+
  xlab(NULL)+
  ylim(0,1)+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Extinction background",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


# statistical models to evaluate variation --------------------------------

####merging both datasets (ages.unnest and ages.stochastic)

###eliminate trees which gave troubles

##high extinction
ages.high.nest.st <- ages.high.nest[lapply(ages.high.stochastic,
                                           class) != "try-error",]

##assigning tree name
ages.high.nest.st$tree <- paste0("tree.", 1:nrow(ages.high.nest.st))


##intermediate
ages.int.nest.st <- ages.int.nest[lapply(ages.int.stochastic,
                                         class) != "try-error",]

##tree name
ages.int.nest.st$tree <- paste0("tree.", 1:nrow(ages.int.nest.st))

##low
ages.low.nest.st <- ages.low.nest[lapply(ages.low.stochastic,
                                         class) != "try-error",]
##tree name
ages.low.nest.st$tree <- paste0("tree.", 1:nrow(ages.low.nest.st))

##rbind
ages.nest.st <- rbind(ages.high.nest.st, ages.int.nest.st, ages.low.nest.st)


##unnesting 

ages.unnest.st <- ages.nest.st %>% unnest(data)

ages.unnest.st$extinction <- as.factor(ages.unnest.st$extinction)


##create a new column needed for mergin in ages.unnest
ages.unnest.st$species <- paste0(ages.unnest.st$label, ".", ages.unnest.st$tree,
                                 ".",   ages.unnest.st$extinction)

ages.stochastic$species <- paste0(ages.stochastic$species, ".",
                                  ages.stochastic$extinction)
##selecting variables
ages.un <- select(ages.unnest.st, species, rTrue.age, rPhylo.age, status,
                  tree) %>% 
                rename(status.c = status)



#######merging

ages.total <- left_join(ages.stochastic, ages.un, by = "species") %>% 
            select(species, tree.x, rTrue.age, rPhylo.age, rStochastic, 
                     root.age, Probability, effect_size, extinction, status.c)%>% 
              rename(tree = tree.x,
                     status = status.c) %>% 
               drop_na()

ages.total$random.age <- sample(ages.total$rTrue.age)


ages.total %>% mutate(dif.sto = abs(rStochastic-rTrue.age)) %>% 
              ggplot(aes(x = Probability, y = dif.sto))+
              geom_point(alpha = 0.5)

##low level functions
##true.age
model.true <- function(df){
  lm(status ~ rTrue.age, data = df)
}

##phylo age
model.phylo <- function(df){
  lm(status ~rPhylo.age, data = df)
}

##stochastic age
model.st <- function(df){
  lm(status ~ rStochastic, data = df)
}

##stochastic age weighted
model.st.wg <- function(df){
  lm(status ~ rStochastic, data = df, weights=Probability)
}

##random
model.random <- function(df){
  lm(status ~ random.age, data = df)
}

ages.models <- ages.total %>% group_by(tree, extinction, effect_size) %>% 
                nest() %>% 
                mutate(model.tr = map(data, model.true),
                       model.phy = map(data, model.phylo),
                       model.sto = map(data, model.st),
                       model.random = map(data, model.random),
                       model.sto.wg = map(data, model.st.wg))

ages.m.glance <- ages.models %>% mutate(glance.tr = map(model.tr, broom::glance),
                                  glance.phy = map(model.phy, broom::glance),
                          glance.st = map(model.sto, broom::glance),
                    glance.random = map(model.random, broom::glance),
                    glance.st.wg = map(model.sto.wg, broom :: glance)) %>%
         select(-c(model.tr, model.phy, model.sto, model.random, model.sto.wg))


for(i in 1:nrow(ages.m.glance)){
  ages.m.glance$beta.true[i] <- ages.models[[5]][[i]][["coefficients"]][2]
  ages.m.glance$beta.phy[i] <- ages.models[[6]][[i]][["coefficients"]][2]
  ages.m.glance$beta.sto[i] <- ages.models[[7]][[i]][["coefficients"]][2]
  ages.m.glance$beta.random[i] <- ages.models[[8]][[i]][["coefficients"]][2]
  ages.m.glance$beta.sto.wg[i] <- ages.models[[9]][[i]][["coefficients"]][2]
  ages.m.glance$rsq.true[i] <- ages.m.glance[[5]][[i]]$r.squared
  ages.m.glance$rsq.phy[i] <- ages.m.glance[[6]][[i]]$r.squared
  ages.m.glance$rsq.sto[i] <- ages.m.glance[[7]][[i]]$r.squared
  ages.m.glance$rsq.random[i] <- ages.m.glance[[8]][[i]]$r.squared
  ages.m.glance$rsq.sto.wg[i] <- ages.m.glance[[9]][[i]]$r.squared
  ages.m.glance$p.true[i] <- ages.m.glance[[5]][[i]]$p.value
  ages.m.glance$p.phy[i] <- ages.m.glance[[6]][[i]]$p.value
  ages.m.glance$p.sto[i] <- ages.m.glance[[7]][[i]]$p.value
  ages.m.glance$p.random[i] <- ages.m.glance[[8]][[i]]$p.value
  ages.m.glance$p.sto.wg[i] <- ages.m.glance[[9]][[i]]$p.value
}

ages.parameters <- ages.m.glance %>% select(-starts_with("glance"), -data)


# Beta parameter ----------------------------------------------------------
beta.parameters <- ages.parameters %>% pivot_longer(cols = starts_with("beta"),
                                                    names_to = "model",
                              values_to = "beta", names_prefix = "beta") %>%
                              mutate(model = as.factor(model))

beta.parameters$model <- factor(beta.parameters$model,
                                   levels=c(".true", ".phy", ".sto", ".random",
                                            ))


###PNG object
png("text/figures/Figure11.Betas.png", width = 15, height = 17,
    units = "cm", 
    pointsize = 8, res = 300)


##low effect size
p1 <- beta.parameters %>%
     filter(effect_size == "low") %>% 
     ggplot(aes(x = beta, fill = model)) + 
     xlim(0,5)+
     geom_density(alpha = 0.4)+
     geom_vline(xintercept = 1, linewidth = 0.75, colour = "red")+
     scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#bdbdbd"),
                            
                    breaks = c('.true', '.phy', ".sto", ".random"),
                    labels=c('True','Phylogenetic', 'Stochastic', 'Random'))+
     facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                               intermediate = "Intermediate",
                                               low = "Low")))+
     ylab("Low")+
     xlab(NULL)+
     theme_bw()+
     xlim(0,3)+
     labs(fill = NULL)+
     mynamestheme

##intermediate effect size
p2 <- beta.parameters %>%
      filter(effect_size == "intermediate") %>% 
      ggplot(aes(x = beta, fill = model)) + 
      xlim(0,5)+
      geom_density(alpha = 0.5)+
      geom_vline(xintercept = 2, linewidth = 0.75, colour = "red")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#bdbdbd"),
                    
                    breaks = c('.true', '.phy', ".sto", ".random"),
                    labels=c('True','Phylogenetic', 'Stochastic', 'Random'))+
      facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                             intermediate = "Intermediate",
                                                   low = "Low")))+
      ylab("Intermediate")+
      xlab(NULL)+
      xlim(0,4)+
      theme_bw()+
      theme(legend.position='none')+
      labs(fill = NULL)+
      mynamestheme

###high effect size
p3 <- beta.parameters %>%
  filter(effect_size == "high") %>% 
  ggplot(aes(x = beta, fill = model)) + 
  xlim(0,5)+
  geom_density(alpha = 0.7)+
  geom_vline(xintercept = 3, linewidth = 0.75, colour = "red")+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#bdbdbd"),
                    
                    breaks = c('.true', '.phy', ".sto", ".random"),
                    labels=c('True','Phylogenetic', 'Stochastic', 'Random'))+
              
  facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                                   intermediate = "Intermediate",
                                                   low = "Low")))+
  ylab("High")+
  xlab(NULL)+
  theme_bw()+
  labs(fill = NULL)+
  mynamestheme


##setting top
top <- text_grob("Extinction background",
                 family = "serif", size = 13,  face = "bold")

left <- text_grob("Effect size",
                  family = "serif", size = 14,  face = "bold", rot = 90)

bottom <- text_grob("Beta",
                    family = "serif", size = 13,  face = "bold")



##grid for plotting both figures
grid_arrange_shared_legend(p1, p2, p3, nrow = 3, ncol = 1, 
                           top = top, bottom = bottom, left = left)

dev.off()

# R-squared ---------------------------------------------------------------

rsq.parameters <- ages.parameters %>% 
                   pivot_longer(cols = starts_with("rsq"), names_to = "model",
                                      values_to = "rsq", names_prefix = "rsq")
  
rsq.parameters$model <- as.factor(rsq.parameters$model)

rsq.parameters$model <- factor(rsq.parameters$model,
                                levels=c(".true", ".phy", ".sto", ".random"))

###PNG object
png("text/figures/Figure12.RSQ.png", width = 15, height = 17,
    units = "cm", 
    pointsize = 8, res = 300)


##low effect size
p1 <- rsq.parameters %>%
      filter(effect_size == "low",
          model != ".true") %>% 
     ggplot(aes(x = rsq, fill = model)) + 
  #xlim(0,5)+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values = c("#d95f02", "#7570b3", "#bdbdbd"),
                    breaks = c('.phy', ".sto", ".random"),
                    labels=c('Phylogenetic', 'Stochastic', 'Random'))+
    facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                               intermediate = "Intermediate",
                                               low = "Low")))+
  ylab("Low")+
  xlab(NULL)+
  theme_bw()+
  labs(fill = NULL)+
  mynamestheme

##intermediate effect size
p2 <- rsq.parameters %>%
      filter(effect_size == "intermediate",
         model != ".true") %>% 
      ggplot(aes(x = rsq, fill = model)) + 
  #xlim(0,5)+
     geom_density(alpha = 0.5)+
     scale_fill_manual(values = c("#d95f02", "#7570b3", "#bdbdbd"),
                    breaks = c('.phy', ".sto", ".random"),
                    labels=c('Phylogenetic', 'Stochastic', 'Random'))+
     facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                                  intermediate = "Intermediate",
                                                   low = "Low")))+
    ylab("Intermediate")+
    xlab(NULL)+
    theme_bw()+
    labs(fill = NULL)+
    mynamestheme

###high effect size
p3 <- rsq.parameters %>%
      filter(effect_size == "high",
         model != ".true") %>% 
      ggplot(aes(x = rsq, fill = model)) + 
  #xlim(0,5)+
     geom_density(alpha = 0.5)+
     scale_fill_manual(values = c("#d95f02", "#7570b3", "#bdbdbd"),
                    breaks = c('.phy', ".sto", ".random"),
                    labels=c('Phylogenetic', 'Stochastic', 'Random'))+
    facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                                   intermediate = "Intermediate",
                                                   low = "Low")))+
    ylab("High")+
    xlab(NULL)+
    theme_bw()+
    labs(fill = NULL)+
    mynamestheme


##setting top
top <- text_grob("Extinction background",
                 family = "serif", size = 13,  face = "bold")

left <- text_grob("Effect size",
                  family = "serif", size = 14,  face = "bold", rot = 90)

bottom <- text_grob("R-squared",
                    family = "serif", size = 13,  face = "bold")



##grid for plotting both figures
grid_arrange_shared_legend(p1, p2, p3, nrow = 3, ncol = 1, 
                           top = top, bottom = bottom, left = left)

dev.off()



# calculating significance of true age----------------------------------------

t.ages <- ages.parameters %>%
      mutate(p.random_sig = case_when(p.random <= 0.05 ~ "significant",
                                      p.random > 0.05 ~ "no_sig")) %>% 
      group_by(extinction, effect_size, p.random_sig) %>% 
      count() %>% 
      pivot_wider(names_from = p.random_sig, values_from = n, 
                  values_fill = 0) %>% 
      summarise(percentage = (significant/(no_sig + significant))*100)
  

  