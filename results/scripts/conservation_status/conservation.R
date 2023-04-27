##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

# Simulations for testing the conservation status  ------------------------

#####################simulating trees####################################
set.seed(35)

###high extinction scenario
trees.high <- sim.bd.taxa(n = 100, numbsim = 1000, lambda = 0.3, mu = 0.25)

##intermediate extinction scenario
trees.int <- sim.bd.taxa(n = 100, numbsim = 1000, lambda = 0.3, mu = 0.15)

##low extinction scenario
trees.low <- sim.bd.taxa(n = 100, numbsim = 1000, lambda = 0.3, mu = 0.05)


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
ages.high <- Map(getAgesExtantSpecies, tax.high, trees.high, Tol = 1e-6) 
                                              
ages.high <- do.call("rbind", ages.high) %>% 
                    rename(True.age = Age,
                      Estimated.age = tip_length_reconstr) %>%
                          select(label, True.age, Estimated.age)


##root age
ages.high$root.age <- rep(root.high, each = 100)

##tree numbers
ages.high$tree <- rep(paste0("tree.", 1:1000), each = 100)

##scale ages
ages.high <- ages.high %>% mutate(rTrue.age = True.age/root.age,
                                  rPhylo.age = Estimated.age/root.age)

##extinction
ages.high$extinction <- rep("high", nrow(ages.high))

##species name
ages.high$species <- paste0(ages.high$label,".", ages.high$tree)


############intermediate extinction
ages.int <- Map(getAgesExtantSpecies, tax.int, trees.int, Tol = 1e-6) 

ages.int <- do.call("rbind", ages.int) %>% 
             rename(True.age = Age,
                    Estimated.age = tip_length_reconstr) %>%
             select(label, True.age, Estimated.age)

##root age
ages.int$root.age <- rep(root.int, each = 100)

##tree numbers
ages.int$tree <- rep(paste0("tree.", 1:1000), each = 100)

##scale ages
ages.int <- ages.int %>% mutate(rTrue.age = True.age/root.age,
                                rPhylo.age = Estimated.age/root.age)

##extinction
ages.int$extinction <- rep("intermediate", nrow(ages.int))

##species name
ages.int$species <- paste0(ages.int$label,".", ages.int$tree)


############low extinction
ages.low <- Map(getAgesExtantSpecies, tax.low, trees.low, , Tol = 1e-6) 

ages.low <- do.call("rbind", ages.low) %>% 
                    rename(True.age = Age,
                    Estimated.age = tip_length_reconstr) %>%
                    select(label, True.age, Estimated.age)

##root ages
ages.low$root.age <- rep(root.low, each = 100)

##tree numbers
ages.low$tree <- rep(paste0("tree.", 1:1000), each = 100)

##scale ages
ages.low <- ages.low %>% mutate(rTrue.age = True.age/root.age,
                                rPhylo.age = Estimated.age/root.age)


##extinction
ages.low$extinction <- rep("low", nrow(ages.low))

##species name
ages.low$species <- paste0(ages.low$label,".", ages.low$tree)


# geometric function -----------------------------------------------------


##high extinction
ages.geometric.high <- lapply(trees.extant.high, getGeometricAges)

##intermediate extinction
ages.geometric.int <- lapply(trees.extant.int, getGeometricAges)

##low extinction
ages.geometric.low <- lapply(trees.extant.low, getGeometricAges)


#####organizing geometric results############

##assigning tree number and root age
for(i in seq_along(ages.geometric.high)){
  ages.geometric.high[[i]]$tree <- rep(paste0("tree.", i),
                                        nrow(ages.geometric.high[[i]]))
  ages.geometric.high[[i]]$species <- rownames(ages.geometric.high[[i]])
  ages.geometric.high[[i]]$root.age <- rep(root.high[i],
                                            nrow(ages.geometric.high[[i]]))
}

##unlist
ages.high.st <- do.call("rbind", ages.geometric.high) 

ages.high.st$species <- paste0(ages.high.st$species,".", ages.high.st$tree)

ages.high.st$extinction <- rep("high", nrow(ages.high.st))

##calculation scaled ages
ages.high.st <- ages.high.st %>% mutate(rmode = modal/root.age,
                                        rmean = mean/root.age)

######intermediate extinction

##assigning tree number and root age
for(i in seq_along(ages.geometric.int)){
  ages.geometric.int[[i]]$tree <- rep(paste0("tree.", i),
                                         nrow(ages.geometric.int[[i]]))
  ages.geometric.int[[i]]$species <- rownames(ages.geometric.int[[i]])
  ages.geometric.int[[i]]$root.age <- rep(root.int[i],
                                             nrow(ages.geometric.int[[i]]))
  
}

##unlist
ages.int.st <- do.call("rbind", ages.geometric.int) 

ages.int.st$species <- paste0(ages.int.st$species,".", ages.int.st$tree)

ages.int.st$extinction <- rep("intermediate", nrow(ages.int.st))


##scaled ages
ages.int.st <- ages.int.st %>% mutate(rmode = modal/root.age,
                                      rmean = mean/root.age)

#####low extinction

##assigning tree number and root age
for(i in seq_along(ages.geometric.low)){
  ages.geometric.low[[i]]$tree <- rep(paste0("tree.", i),
                                         nrow(ages.geometric.low[[i]]))
  ages.geometric.low[[i]]$species <- rownames(ages.geometric.low[[i]])
  ages.geometric.low[[i]]$root.age <- rep(root.low[i],
                                             nrow(ages.geometric.low[[i]]))
  
 
  
}


##unlist
ages.low.st <- do.call("rbind", ages.geometric.low) 

ages.low.st$species <- paste0(ages.low.st$species,".", ages.low.st$tree)

ages.low.st$extinction <- rep("low", nrow(ages.low.st))

ages.low.st <- ages.low.st %>% mutate(rmode = modal/root.age,
                                      rmean = mean/root.age)


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

##saving
write_csv(ages.total, file = file.path(getwd(), pro,
                                       c_status, "ages.total.csv"))

###calling the csv file
ages.total <- read_csv(file = file.path(getwd(),
                                        pro, c_status, "ages.total.csv"))

##as factor
ages.total$status <- as.factor(ages.total$status)

ages.total$status <- factor(ages.total$status,
                                 levels = c("LC", "NT", "VU", "EN", "CR"))

ages.total$extinction <- as.factor(ages.total$extinction)

ages.total$extinction <- factor(ages.total$extinction,
                                levels = c("low", "intermediate", "high"))

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

#####True age vs Phylo age and Conservation status
png("text/figures/Figure7.conservation_true.vs.phylo.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High extinction",
                                    intermediate = "Intermediate extinction",
                                    low = "Low extinction")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme
  

p2 <- ggplot(ages.total, aes(x = status, y = log(Estimated.age+1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High extinction",
                                    intermediate = "Intermediate extinction",
                                      low = "Low extinction")))+
  ylab("log(Phylogenetic age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Positive effect", family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


######True age vs Highest prob age and Conservation status
png("text/figures/Figure8.conservation_true.vs.mode.png", width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                  intermediate = "Intermediate extinction",
                                  high = "High extinction")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme


p2 <- ggplot(ages.total, aes(x = status, y = log(modal + 1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  ylab("log(Most probable age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Positive effect",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()



######True age vs Mean prob age and Conservation status
png("text/figures/Figure8.1.conservation.true.vs.mean.png", 
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                       intermediate = "Intermediate extinction",
                                       high = "High extinction")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme


p2 <- ggplot(ages.total, aes(x = status, y = log(mean+ 1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                     intermediate = "Intermediate extinction",
                                     high = "High extinction")))+
  ylab("log(Mean probable age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("Positive effect",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()

##########################No effects################################

##randomize rows and then assign the status

##extracting the status column 
status <- ages.total$status

##removing the status column from the dataframe and shuffling the rows
ages.total.random <- ages.total %>% select(-status) %>% sample_frac()
 
##joining the status column                          
ages.total.random$status <- status

########Figures

#####True age vs Phylo age and Conservation status
png("text/figures/Figure9.conservation_true.vs.phylo.no.effects.png", 
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total.random, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High extinction",
                                     intermediate = "Intermediate extinction",
                                     low = "Low extinction")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme


p2 <- ggplot(ages.total.random, aes(x = status, y = log(Estimated.age+1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(high = "High extinction",
                                     intermediate = "Intermediate extinction",
                                     low = "Low extinction")))+
  ylab("log(Phylogenetic age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("No effect", family = "serif", 
                 size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()


######True age vs Highest prob age and Conservation status
png("text/figures/Figure9.1.conservation_true.vs.mode.no.effect.png",
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total.random, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                       intermediate = "Intermediate extinction",
                                       high = "High extinction")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme


p2 <- ggplot(ages.total.random, aes(x = status, y = log(modal + 1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                       intermediate = "Intermediate extinction",
                                       high = "High extinction")))+
  ylab("log(Most probable age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("No effect",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()



######True age vs Mean prob age and Conservation status
png("text/figures/Figure9.2.conservation.true.vs.mean.no.effect.png", 
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

p1 <- ggplot(ages.total.random, aes(x = status, y = log(True.age+1)))+
  geom_boxplot(outlier.shape  = NA)+
  geom_jitter(color="black", size=0.4, alpha = 0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                       intermediate = "Intermediate extinction",
                                       high = "High extinction")))+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme


p2 <- ggplot(ages.total.random, aes(x = status, y = log(mean + 1)))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(color="black", size=0.4, alpha=0.1)+
  facet_wrap(~extinction, labeller = as_labeller(c(low = "Low extinction",
                                       intermediate = "Intermediate extinction",
                                       high = "High extinction")))+
  ylab("log(Mean probable age + 1)")+
  xlab("Conservation status")+
  theme_bw()+
  mynamestheme

##setting top
top <- text_grob("No effect",
                 family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(p1, p2, nrow = 2, top = top)

dev.off()

##########line plots##########################################################

##Calculating the mean of each age for each tree and extinction background
ages.mean <-ages.total%>% 
                 group_by(extinction, tree, status) %>%
                 summarise(mean.true = mean(True.age),
                           mean.phy = mean(Estimated.age),
                           mean.mode = mean(modal),
                           mean.mean = mean(mean)) 

##renaming tree for grouping
ages.mean$tree <- paste0(ages.mean$tree, ".", 
                              ages.mean$extinction)

#######ages ranking of correct relationship between age and signal
ages.rank <- ages.mean %>% group_by(tree, extinction) %>%
                     mutate(true.rank = dense_rank(mean.true),
                         phy.rank = dense_rank(mean.phy),
                         mode.rank = dense_rank(mean.mode),
                         mean.rank = dense_rank(mean.mean)) %>% 
                      mutate(resp_phy = if_else(true.rank == phy.rank, 1, 0),
                         resp_mode = if_else(true.rank == mode.rank, 1, 0),
                         resp_mean = if_else(true.rank == mean.rank, 1, 0)) 
  
  
  ###comparing phylogenetic age
  ages.comparison.phy <- ages.rank %>% group_by(tree, extinction) %>% 
  summarise(sum_phy = sum(resp_phy)) %>% 
  filter(sum_phy >= 3) %>% 
  group_by(extinction) %>% 
  count()


###comparing mode age from geometric function
ages.comparison.mode <- ages.rank %>% group_by(tree, extinction) %>% 
  summarise(sum_mode = sum( resp_mode)) %>% 
  filter(sum_mode >= 3) %>% 
  group_by(extinction) %>% 
  count()


####comparing mean age from geometric function
ages.comparison.mean <- ages.rank %>% group_by(tree, extinction) %>% 
  summarise(sum_mean = sum(resp_mean)) %>% 
  filter(sum_mean >= 3) %>% 
  group_by(extinction) %>% 
  count()

##mean true age vs mean phylo age

png("text/figures/Figure10.line.true.vs.phy.png", 
    width = 17, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)



f1 <-ggplot(ages.mean, aes(x = status, y = log(mean.true +1), group = tree,
                           color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.1)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                    intermediate = "Intermediate extinction",
                                             high = "High extinction")))+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  


##mean phylo age
f2 <-ggplot(ages.mean, aes(x = status, y = log(mean.phy +1), group = tree,
                           colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape = 21, fill="white")+
  geom_line(alpha = 0.1)+
  geom_text(data = ages.comparison.phy, aes(x = 2.05, y = 2.22,
                                label = paste0("Correct estimation:", " ",
                             round(n/10),"%")),
            inherit.aes = FALSE,
            size = 3, 
            family = "serif")+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                        intermediate = "Intermediate extinction",
                                        high = "High extinction")))+
  theme_bw()+
  ylab("log(Phylogenetic age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  

##setting top
top <- text_grob("Positive effect",
                 family = "serif", size = 13,  face = "bold")

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(f1, f2, nrow = 2, top = top, bottom= bottom)

dev.off()


##############highest probabiltiy age
png("text/figures/Figure10.1.line.true.vs.mode.png", 
    width = 17, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)


##mean true age
f1 <-ggplot(ages.mean, aes(x = status, y = log(mean.true +1), group = tree,
                           color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.1)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  


##mean mode age
f2 <-ggplot(ages.mean, aes(x = status, y = log(mean.mode +1), group = tree,
                           colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.1)+
  geom_text(data = ages.comparison.mode, aes(x = 2.05, y = 2.22,
                                  label = paste0("Correct estimation:", " ",
                                                         round(n/10),"%")),
            inherit.aes = FALSE,
            size = 3, 
            family = "serif")+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(Most probable age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")

##setting top
top <- text_grob("Positive effect",
                 family = "serif", size = 13,  face = "bold")

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(f1, f2, nrow = 2, top = top, bottom= bottom)

dev.off()

#########mean probability age
png("text/figures/Figure10.2.line.true.vs.mean.png", 
    width = 17, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)


##mean true age
f1 <-ggplot(ages.mean, aes(x = status, y = log(mean.true +1), group = tree,
                           color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.1)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  


##mean mean age
f2 <-ggplot(ages.mean, aes(x = status, y = log(mean.mean +1), group = tree,
                           colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.1)+
  geom_text(data = ages.comparison.mean, aes(x = 2.05, y = 3,
                                    label = paste0("Correct estimation:", " ",
                                                           round(n/10),"%")),
            inherit.aes = FALSE,
            size = 3, 
            family = "serif")+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(Mean probable age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  

##setting top
top <- text_grob("Positive effect",
                 family = "serif", size = 13,  face = "bold")

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(f1, f2, nrow = 2, top = top, bottom= bottom)

dev.off()




#####################line plots no effect#################################

##Calculating the mean of each age for each tree and extinction background
ages.mean.random <-ages.total.random %>% 
                                  group_by(extinction, tree, status) %>%
                                  summarise(mean.true = mean(True.age),
                                  mean.phy = mean(Estimated.age),
                                  mean.mode = mean(modal),
                                  mean.mean = mean(mean)) 

##renaming tree for grouping
ages.mean.random$tree <- paste0(ages.mean.random$tree, ".", 
                         ages.mean.random$extinction)



###mean true age
png("text/figures/Figure11.line.true.vs.phy.no.effect.png", 
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)


##mean true age
f1 <-ggplot(ages.mean.random, aes(x = status, y = log(mean.true +1), group = tree,
                           color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.5)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  theme(legend.position = "none")+
  mynamestheme


##mean phylo age
f2 <-ggplot(ages.mean.random, aes(x = status, y = log(mean.phy +1), group = tree,
                           colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.5)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(Phylogenetic age + 1)")+
  xlab(NULL)+
  theme(legend.position = "none")+
  mynamestheme

##setting top
top <- text_grob("No effect",
                 family = "serif", size = 13,  face = "bold")

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(f1, f2, nrow = 2, top = top, bottom= bottom)

dev.off()


##############highest probabiltiy age
png("text/figures/Figure11.1.line.true.vs.mode.no.effect.png", 
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)


##mean true age
f1 <-ggplot(ages.mean.random, aes(x = status, y = log(mean.true +1), group = tree,
                           color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.5)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  


##mean mode age
f2 <-ggplot(ages.mean.random, aes(x = status, y = log(mean.mode +1), group = tree,
                           colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.5)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(Most probable age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  

##setting top
top <- text_grob("No effect",
                 family = "serif", size = 13,  face = "bold")

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(f1, f2, nrow = 2, top = top, bottom= bottom)

dev.off()

#########mean probability age
png("text/figures/Figure11.2.line.true.vs.mean.no.effect.png", 
    width = 15, height = 15,
    units = "cm", 
    pointsize = 8, res = 300)


##mean true age
f1 <-ggplot(ages.mean.random, aes(x = status, y = log(mean.true +1), group = tree,
                           color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.5)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(True age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  


##mean phylo age
f2 <-ggplot(ages.mean.random, aes(x = status, y = log(mean.mean +1), group = tree,
                           colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.5)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("log(Mean probable age + 1)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")
  

##setting top
top <- text_grob("No effect",
                 family = "serif", size = 13,  face = "bold")

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(f1, f2, nrow = 2, top = top, bottom= bottom)

dev.off()




                                                         





  