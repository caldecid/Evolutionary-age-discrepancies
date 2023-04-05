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
load(file = file.path(getwd(), pro, c_status,"trees.conservation.RData"))
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

##extinction
ages.high$extinction <- rep("high", nrow(ages.high))

##species name
ages.high$species <- paste0(ages.high$label,".", ages.high$tree)


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

##extinction
ages.int$extinction <- rep("intermediThe ate", nrow(ages.int))

##species name
ages.int$species <- paste0(ages.int$label,".", ages.int$tree)


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


##extinction
ages.low$extinction <- rep("low", nrow(ages.low))

##species name
ages.low$species <- paste0(ages.low$lThe abel,".", ages.low$tree)


# stochastic function -----------------------------------------------------

##future package conditions 
plan("multisession", workers = 3)

##high extinction
ages.stochastic.high <- future_lapply(trees.extant.high,
                                               stochastic.sp.age.torsten)

##intermediate extinction
ages.stochastic.int <- future_lapply(trees.extant.int,
                                               stochastic.sp.age.torsten)

##low extinction
ages.stochastic.low <- future_lapply(trees.extant.low,
                                     stochastic.sp.age.torsten)


#####organizing stochastic results############

##assigning tree number and root age
for(i in seq_along(ages.stochastic.high)){
  ages.stochastic.high[[i]]$tree <- rep(paste0("tree.", i),
                                        nrow(ages.stochastic.high[[i]]))
  ages.stochastic.high[[i]]$root.age <The - rep(root.high[i],
                                            nrow(ages.stochastic.high[[i]]))
}

##unlist
ages.high.st <- do.call("rbind", ages.stochastic.high) 

ages.high.st$species <- paste0(ages.high.st$species,".", ages.high.st$tree)

ages.high.st$extinction <- rep("high", nrow(ages.high.st))

ages.high.st <- ages.high.st %>% unnest(dist_ages) %>% 
  filter(Probability == max(Probability)) %>% 
  rename(max.prob.age = Age)

ages.high.st$max.prob.age <- as.numeric(as.character(ages.high.st$max.prob.age))

ages.high.st <- ages.high.st %>% mutate(rStochastic = max.prob.age/root.age)

######intermediate extinction

##assigning tree number and root age
for(i in seq_along(ages.stochastic.intThe )){
  ages.stochastic.int[[i]]$tree <- rep(paste0("tree.", i),
                                         nrow(ages.stochastic.int[[i]]))
  ages.stochastic.int[[i]]$root.age <- rep(root.int[i],
                                             nrow(ages.stochastic.int[[i]]))
  
}

##unlist
ages.int.st <- do.call("rbind", ages.stochastic.int) 

ages.int.st$species <- paste0(ages.int.st$species,".", ages.int.st$tree)

ages.int.st$extinction <- rep("intermediate", nrow(ages.int.st))

ages.int.st <- ages.int.st %>% unnest(dist_ages) %>% 
  filter(Probability == max(Probability)) %>% 
  rename(max.prob.age = Age)

ages.int.st$max.prob.age <- as.numeric(as.character(ages.int.st$max.prob.age))

ages.int.st <- ages.int.st %>% mutate(rStochastic = max.prob.age/root.age)

#####low extinction

##assigning tree number and root age
for(i in seq_along(ages.stochastic.low)){
  ages.stochastic.low[[i]]$tree <- rep(paste0("tree.", i),
                                         nrow(ages.stochastic.low[[i]]))
  ages.stochastic.low[[i]]$root.age <- rep(root.low[i],
                                             nrow(ages.stochastic.low[[i]]))
  
 
  
}


##unlist
ages.low.st <- do.call("rbind", ages.stochastic.low) 

ages.low.st$species <- paste0(ages.low.st$species,".", ages.low.st$tree)

ages.low.st$extinction <- rep("low", nrow(ages.low.st))

ages.low.st <- ages.low.st %>% unnest(dist_ages) %>% 
  filter(Probability == max(Probability)) %>% 
  rename(max.prob.age = Age) 

ages.low.st$max.prob.age <- as.numeric(as.character(ages.low.st$max.prob.age))

ages.low.st <- ages.low.st %>% mutate(rStochastic = max.prob.age/root.age)


# simulating the extinction signal -----------------------------------------

####high extinction
ages.high$status <- discretize(ages.high$rTrue.age, 
                               method = "frequency", 
                               breaks = 5,
                               labels = c("LC", "NT", "VU", "EN", "CR")) 

ages.high.total <- left_join(ages.high, ages.high.st, by = "species")

####intermediate extinction
ages.int$status <- discretize(ages.int$rTrue.age, 
                               method = "frequency", 
                               breaks = 5,
                               labels = c("LC", "NT", "VU", "EN", "CR")) 


ages.int.total <- left_join(ages.int, ages.int.st, by = "species")

#########low extinction
ages.low$status <- discretize(ages.low$rTrue.age, 
                              method = "frequency", 
                              breaks = 5,
                              labels = c("LC", "NT", "VU", "EN", "CR")) 


ages.low.total <- left_join(ages.low, ages.low.st, by = "species")

###rbinding the dfs and cleaning the dataset
ages.total <- rbind(ages.high.total, ages.int.total, ages.low.total) %>% 
              select(-c(extinction.y, tree.y, root.age.y)) %>% 
              rename(extinction = extinction.x,
                     tree = tree.x,
                     root.age = root.age.x)

write_csv(ages.total, file = file.path(getwd(), pro, c_status, "ages.total.csv"))

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
                                                  size=.5, linetype="dotted"),
                      legend.position = "bottom")

########################Effect size
png("text/figures/Figure7.conservation_true.vs.phylo.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High",
                                                intermediate = "Intermediate",
                                                low = "Low")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme
  

p2 <- ggplot(ages.total, aes(x = status, y = log(Estimated.age+1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High",
                                                intermediate = "Intermediate",
                                                  low = "Low")))+
  ylab("log(Phylogenetic age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Extinction background", family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


########################extinction background
png("text/figures/Figure8.conservation_true.vs.sto.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low",
                                                    intermediate = "Intermediate",
                                                    high = "High")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme


p2 <- ggplot(ages.total, aes(x = status, y = log(max.prob.age + 1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low",
                                                intermediate = "Intermediate",
                                                high = "High")))+
  ylab("log(Highest probable age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Extinction Background",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()



ggplot(ages.total, aes(x = abs(log(True.age+1) - log(max.prob.age+1)), y = Probability))+
  geom_point()

# statistical models to evaluate variation --------------------------------

##low level functions
##true.age
model.true <- function(df){
  lm(log(True.age+1) ~ status, data = df)
}

##phylo age
model.phylo <- function(df){
  lm(log(Estimated.age+1) ~ status, data = df)
}

##stochastic age
model.st <- function(df){
  lm(log(max.prob.age+1) ~ status, data = df)
}

###nesting the models

ages.models <- ages.total %>% group_by(tree, extinction) %>% 
                nest() %>% 
                mutate(model.tr = map(data, model.true),
                       model.phy = map(data, model.phylo),
                       model.sto = map(data, model.st))
                    
##glancing the models
ages.m.glance <- ages.models %>% mutate(glance.tr = map(model.tr, broom::glance),
                                  glance.phy = map(model.phy, broom::glance),
                          glance.st = map(model.sto, broom::glance)) %>% 
                 select(-c(model.tr, model.phy, model.sto))


for(i in 1:nrow(ages.m.glance)){
  ages.m.glance$beta.true[i] <- ages.models[[4]][[i]][["coefficients"]][2]
  ages.m.glance$beta.phy[i] <- ages.models[[5]][[i]][["coefficients"]][2]
  ages.m.glance$beta.sto[i] <- ages.models[[6]][[i]][["coefficients"]][2]
  ages.m.glance$rsq.true[i] <- ages.m.glance[[4]][[i]]$r.squared
  ages.m.glance$rsq.phy[i] <- ages.m.glance[[5]][[i]]$r.squared
  ages.m.glance$rsq.sto[i] <- ages.m.glance[[6]][[i]]$r.squared
  ages.m.glance$p.true[i] <- ages.m.glance[[4]][[i]]$p.value
  ages.m.glance$p.phy[i] <- ages.m.glance[[5]][[i]]$p.value
  ages.m.glance$p.sto[i] <- ages.m.glance[[6]][[i]]$p.value
}

###models parameters
ages.parameters <- ages.m.glance %>% select(-starts_with("glance"), -data)


# Beta parameter ----------------------------------------------------------
beta.parameters <- ages.parameters %>% pivot_longer(cols = starts_with("beta"),
                                                    names_to = "model",
                              values_to = "beta", names_prefix = "beta") %>%
                              mutate(model = as.factor(model))

beta.parameters$model <- factor(beta.parameters$model,
                                   levels=c(".true", ".phy", ".sto"))


###PNG object
png("text/figures/Figure11.Betas.png", width = 15, height = 17,
    units = "cm", 
    pointsize = 8, res = 300)


##low effect size
beta.parameters %>%
      ggplot(aes(x = beta, fill = model)) + 
     xlim(0,5)+
     geom_density(alpha = 0.4)+
     scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                            
                    breaks = c('.true', '.phy', ".sto"),
                    labels=c('True','Phylogenetic', 'Stochastic'))+
     facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                               intermediate = "Intermediate",
                                               low = "Low")))+
     ylab(NULL)+
     xlab("Beta")+
     theme_bw()+
     xlim(0,1)+
     labs(fill = "Model")+
    ggtitle("Extinction background")+
     mynamestheme

dev.off()

# R-squared ---------------------------------------------------------------

rsq.parameters <- ages.parameters %>% 
                   pivot_longer(cols = starts_with("rsq"), names_to = "model",
                                      values_to = "rsq", names_prefix = "rsq")
  
rsq.parameters$model <- as.factor(rsq.parameters$model)

rsq.parameters$model <- factor(rsq.parameters$model,
                                levels=c(".true", ".phy", ".sto"))

###PNG object
png("text/figures/Figure12.RSQ.png", width = 15, height = 17,
    units = "cm", 
    pointsize = 8, res = 300)


rsq.parameters %>%
      ggplot(aes(x = rsq, fill = model)) + 
  #xlim(0,5)+
    geom_density(alpha = 0.5)+
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    breaks = c('.true','.phy', ".sto"),
                    labels=c('True', 'Phylogenetic', 'Stochastic'))+
    facet_grid(~extinction, labeller = as_labeller(c(high = "High",
                                               intermediate = "Intermediate",
                                               low = "Low")))+
  ylab(NULL)+
  xlab("R-squared")+
  theme_bw()+
  labs(fill = "Model")+
  ggtitle("Extinction background")+
  mynamestheme


dev.off()




  