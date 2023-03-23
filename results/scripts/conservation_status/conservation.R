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


###################prunning extinct species###########################

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
png("Figure7.conservation_effect_size.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(y1, aes(x = status.dis, y = rTrue.age))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~effect_size, labeller = as_labeller(c(high = "High",
                                                intermediate = "Intermediate",
                                                low = "Low")))+
  ylab("True age (scale)")+
  xlab(NULL)+
  mynamestheme
  

p2 <- ggplot(y1, aes(x = status.dis, y = rPhylo.age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~effect_size, labeller = as_labeller(c(high = "High",
                                                intermediate = "Intermediate",
                                                  low = "Low")))+
  ylab("Phylogenetic age (scale)")+
  xlab("Conservation status")+
  mynamestheme

##setting top
top <- text_grob("Effect size", family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


########################extinction background
png("Figure8.conservation_extinction.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(y1, aes(x = status.dis, y = rTrue.age))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low",
                                                    intermediate = "Intermediate",
                                                    high = "High")))+
  ylab("True age (scale)")+
  xlab(NULL)+
  mynamestheme


p2 <- ggplot(y1, aes(x = status.dis, y = rPhylo.age))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low",
                                                intermediate = "Intermediate",
                                                high = "High")))+
  ylab("Phylogenetic age (scale)")+
  xlab("Conservation status")+
  mynamestheme

##setting top
top <- text_grob("Extinction Background",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()



