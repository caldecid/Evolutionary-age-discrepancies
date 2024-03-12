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
#############simulating taxonomies only bifurcating######################

##high extinction
tax.high <- lapply(trees.high, sim.taxonomy, beta = 1)

##intermediate extinction
tax.int <- lapply(trees.int, sim.taxonomy, beta = 1)

##low extinction
tax.low <- lapply(trees.low, sim.taxonomy, beta = 1)

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
ages.low <- Map(getAgesExtantSpecies, tax.low, trees.low,  Tol = 1e-6) 

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


# New function -----------------------------------------------------

#mean

high.mean.age <- vector(mode = "numeric",
                        length = nrow(ages.high))

for(i in 1:nrow(ages.high)){
  high.mean.age[[i]] <- meanAgeFunction(lam = 0.3,
                                        mu = 0.25,
                                        rho = 1,
                                        v = ages.high$Estimated.age[i])
}



##median
high.median.age <- vector(mode = "numeric",
                          length = nrow(ages.high))

for(i in 1:nrow(ages.high)){
  high.median.age[[i]] <- medianAgeFunction(lam = 0.3,
                                            mu = 0.25,
                                            rho = 1,
                                            v = ages.high$Estimated.age[i])
}


##joining dataframes
ages.high.bif <- cbind(ages.high, high.mean.age) %>% 
                    rename(mean_age = high.mean.age)


#########intermediate extinction scenario
#mean

int.mean.age  <- vector(mode = "numeric",
                        length = nrow(ages.int))

for(i in 1:nrow(ages.int)){
  int.mean.age [[i]] <- meanAgeFunction(lam = 0.3,
                                        mu = 0.15,
                                        rho = 1,
                                        v = ages.int$Estimated.age[i])
}

##median
int.median.age <- vector(mode = "numeric",
                          length = nrow(ages.int))


for(i in 1:nrow(ages.int)){
  int.median.age[[i]] <- medianAgeFunction(lam = 0.3,
                                            mu = 0.15,
                                            rho = 1,
                                            v = ages.int$Estimated.age[i])
}



##joining dataframes
ages.int.bif <- cbind(ages.int, int.mean.age) %>% 
                rename(mean_age = int.mean.age)


#######low extinction scenario
#mean

low.mean.age <- vector(mode = "numeric",
                        length = nrow(ages.low))

for(i in 1:nrow(ages.low)){
  low.mean.age[[i]] <- meanAgeFunction(lam = 0.3,
                                        mu = 0.05,
                                        rho = 1,
                                        v = ages.low$Estimated.age[i])
}

##median
low.median.age <- vector(mode = "numeric",
                          length = nrow(ages.low))


for(i in 1:nrow(ages.low)){
  low.median.age[[i]] <- medianAgeFunction(lam = 0.3,
                                            mu = 0.05,
                                            rho = 1,
                                            v = ages.low$Estimated.age[i])
}


##joining dataframes
ages.low.bif <- cbind(ages.low, low.mean.age) %>% 
  rename(mean_age = low.mean.age)

##################merging the dataframes with different extinction#####
ages.total.bif <- rbind(ages.high.bif,ages.int.bif, ages.low.bif)

ages.total.bif$extinction <- factor(ages.total.bif$extinction,
                                    levels = c("low", "intermediate", "high"),
                                    ordered = TRUE)
##saving
write_csv(ages.total.bif, file = file.path(pro, c_status, "ages.total.bif.csv"))

#ages.total.bif <- read_csv(file = file.path(pro, c_status, "ages.total.bif.csv"))

# simulating the extinction signal -----------------------------------------

####high extinction
ages.high.bif <- ages.high.bif %>% group_by(tree) %>%
                mutate(status = discretize(True.age, 
                                     method = "frequency", 
                                     breaks = 5,
                                     labels = c("LC", "NT", "VU", "EN", "CR")))


ages.high.bif$status <- factor(ages.high.bif$status,
                            levels = c("LC", "NT", "VU", "EN", "CR"),
                            ordered= TRUE)

####intermediate extinction
ages.int.bif <- ages.int.bif %>% group_by(tree) %>% 
                          mutate(status = discretize(True.age, 
                                     method = "frequency", 
                                     breaks = 5,
                                     labels = c("LC", "NT", "VU", "EN", "CR")))


ages.int.bif$status <- factor(ages.int.bif$status,
                                 levels = c("LC", "NT", "VU", "EN", "CR"),
                                 ordered= TRUE)


#########low extinction
ages.low.bif <- ages.low.bif %>% group_by(tree)  %>% 
                                     mutate(status = discretize(True.age, 
                                      method = "frequency", 
                                      breaks = 5,
                                      labels = c("LC", "NT", "VU", "EN", "CR")))


ages.low.bif$status <- factor(ages.low.bif$status,
                                 levels = c("LC", "NT", "VU", "EN", "CR"),
                                 ordered= TRUE)


###rbinding the dfs and cleaning the dataset
ages.total.status <- rbind(ages.high.bif, ages.int.bif, ages.low.bif)
      

##saving
write_csv(ages.total.status, file = file.path(getwd(), pro,
                                       c_status, "ages.total.status.csv"))

#ages.total.status = read_csv(file = file.path(getwd(), pro,
 #                                             c_status, "ages.total.status.csv"))



##as factor
ages.total.status$status <- factor(ages.total.status$status,
                       levels = c("LC", "NT", "VU", "EN", "CR"),
                       ordered = TRUE)

ages.total.status$extinction <- factor(ages.total.status$extinction,
                        levels = c("low", "intermediate", "high"),
                        ordered = TRUE)

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


##########line plots##########################################################

##Calculating the mean of each age for each tree and extinction background

##High extinction
ages.mean <-ages.total.status %>% 
                 group_by(extinction,tree, status) %>%
                 summarise(mean.true = mean(True.age),
                           mean.phy = mean(Estimated.age),
                           mean.mean = mean(mean_age)) 



##factors
ages.mean$status <- factor(ages.mean$status,
                           levels = c("LC", "NT", "VU", "EN", "CR"),
                           ordered = TRUE)
##factor
ages.mean$extinction <- factor(ages.mean$extinction, levels = c("low",
                                                                "intermediate",
                                                                "high"),
                               ordered = TRUE)
##changing tree names for grouping
ages.mean$tree <- paste0(ages.mean$tree, ".", ages.mean$extinction)

#######ages ranking of correct relationship between age and signal
ages.rank <- ages.mean %>% group_by(tree, extinction) %>%
                     mutate(true.rank = dense_rank(mean.true),
                         phy.rank = dense_rank(mean.phy),
                         mean.rank = dense_rank(mean.mean)) %>% 
                      mutate(resp_phy = if_else(true.rank == phy.rank, 1, 0),
                         resp_mean = if_else(true.rank == mean.rank, 1, 0)) 
  
  
  ###comparing phylogenetic age
  ages.comparison.phy <- ages.rank %>% 
    group_by(extinction,tree) %>% 
      summarise(sum_phy = sum(resp_phy)) %>% 
      filter(sum_phy < 5) %>% 
      count()%>% 
      mutate(error = n/1000*100)


###comparing mean estimation from the probabilistic function
ages.comparison.mean <- ages.rank %>% group_by(extinction, tree) %>% 
  summarise(sum_mean = sum(resp_mean))%>% 
  filter(sum_mean < 5) %>% 
  count()%>% 
  mutate(error = n/1000*100)



###True age plot

##low
tree.low <- ages.mean %>%
  filter(extinction == "low") %>% pull(tree)

tree.low <- sample(tree.low, size = 100)

##intermediate
tree.int <- ages.mean %>%
  filter(extinction == "intermediate") %>% pull(tree)

tree.int <- sample(tree.int, size = 100)

##high
tree.high <- ages.mean %>%
  filter(extinction == "high") %>% pull(tree)

tree.high <- sample(tree.high, size = 100)  

##binding
trees.true <- c(tree.low, tree.int, tree.high)

###plot
ftrue <- ages.mean %>% filter(tree %in% trees.true) %>% 
  ggplot(aes(x = status, y = log(mean.true +1), group = tree,
             color = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.3)+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("True age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")

##Phylogenetic 

##low extinction

##correct estimation 
phy.low.correct <- ages.rank %>% filter(extinction == "low") %>%
  group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy, na.rm = TRUE)) %>% 
  filter(sum_phy == 5) %>%
  pull(tree)

phy.low.correct <- sample(phy.low.correct, 
                          size = round(length(phy.low.correct)*0.1))


##correct estimation 
phy.low.incorrect <- ages.rank %>% filter(extinction == "low") %>%
  group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy, na.rm = TRUE)) %>% 
  filter(sum_phy < 5) %>%
  pull(tree)

phy.low.incorrect <- sample(phy.low.incorrect,
                            size = length(phy.low.incorrect)*0.1)


##int extinction

##correct estimation 
phy.int.correct <- ages.rank %>% filter(extinction == "intermediate") %>%
  group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy, na.rm = TRUE)) %>% 
  filter(sum_phy == 5) %>%
  pull(tree)

phy.int.correct <- sample(phy.int.correct, 
                          size = round(length(phy.int.correct)*0.1))


##incorrect estimation 
phy.int.incorrect <- ages.rank %>% filter(extinction == "intermediate") %>%
  group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy, na.rm = TRUE)) %>% 
  filter(sum_phy < 5) %>%
  pull(tree)

phy.int.incorrect <- sample(phy.int.incorrect,
                            size = round(length(phy.int.incorrect)*0.1))

##high extinction

##correct estimation 
phy.high.correct <- ages.rank %>% filter(extinction == "high") %>%
  group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy, na.rm = TRUE)) %>% 
  filter(sum_phy == 5) %>%
  pull(tree)

phy.high.correct <- sample(phy.high.correct, 
                           size = round(length(phy.high.correct)*0.1))


##incorrect estimation 
phy.high.incorrect <- ages.rank %>% filter(extinction == "high") %>%
  group_by(tree) %>% 
  summarise(sum_phy = sum(resp_phy, na.rm = TRUE)) %>% 
  filter(sum_phy < 5) %>%
  pull(tree)

phy.high.incorrect <- sample(phy.high.incorrect,
                             size = round(length(phy.high.incorrect)*0.1))

##binding correct
phy.correct <- c(phy.low.correct, phy.int.correct, phy.high.correct)

##filtering correct dataset
ages.phy.correct <- ages.mean %>% filter(tree %in% phy.correct)

##binding incorrect
phy.incorrect <- c(phy.low.incorrect, phy.int.incorrect, phy.high.incorrect)

##filtering incorrect dataset
ages.phy.incorrect <- ages.mean %>% filter(tree %in% phy.incorrect)

##Plot
fphylo <-ggplot(ages.phy.correct, aes(x = status, y = log(mean.phy +1), group = tree,
                                      colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape = 21, fill="white")+
  geom_line(alpha = 0.3)+
  geom_point(data = ages.phy.incorrect,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.phy.incorrect,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.phy, aes(x = 2.05, y = 2.22,
                                            label = paste0("Error rate:", " ",
                                                           error,"%")),
            inherit.aes = FALSE,
            size = 3, 
            family = "serif")+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("Phylogenetic age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")

####Mean point estimate probabilistic function

##low extinction
##correct estimation 
mean.low.correct <- ages.rank %>% filter(extinction == "low") %>%
  group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean == 5) %>%
  pull(tree)

mean.low.correct <- sample(mean.low.correct, 
                           size = round(length(mean.low.correct)*0.1))


##incorrect estimation 
mean.low.incorrect <- ages.rank %>% filter(extinction == "low") %>%
  group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>%
  pull(tree)

mean.low.incorrect <- sample(mean.low.incorrect,
                             size = round(length(mean.low.incorrect)*0.1))


##int extinction

##correct estimation 
mean.int.correct <- ages.rank %>% filter(extinction == "intermediate") %>%
  group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean == 5) %>%
  pull(tree)

mean.int.correct <- sample(mean.int.correct, 
                           size = round(length(mean.int.correct)*0.1))


##incorrect estimation 
mean.int.incorrect <- ages.rank %>% filter(extinction == "intermediate") %>%
  group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>%
  pull(tree)

mean.int.incorrect <- sample(mean.int.incorrect,
                             size = round(length(mean.int.incorrect)*0.1))

##high extinction

##correct estimation 
mean.high.correct <- ages.rank %>% filter(extinction == "high") %>%
  group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean == 5) %>%
  pull(tree)

mean.high.correct <- sample(mean.high.correct, 
                            size = round(length(mean.high.correct)*0.1))


##incorrect estimation 
mean.high.incorrect <- ages.rank %>% filter(extinction == "high") %>%
  group_by(tree) %>% 
  summarise(sum_mean = sum(resp_mean, na.rm = TRUE)) %>% 
  filter(sum_mean < 5) %>%
  pull(tree)

mean.high.incorrect <- sample(mean.high.incorrect,
                              size = round(length(mean.high.incorrect)*0.1))

##binding correct
mean.correct <- c(mean.low.correct, mean.int.correct, mean.high.correct)

##filtering correct dataset
ages.mean.correct <- ages.mean %>% filter(tree %in% mean.correct)

##binding incorrect
mean.incorrect <- c(mean.low.incorrect, mean.int.incorrect, mean.high.incorrect)

##filtering incorrect dataset
ages.mean.incorrect <- ages.mean %>% filter(tree %in% mean.incorrect)


##plot
fprob <-ggplot(ages.mean.correct,
               aes(x = status, y = log(mean.mean +1), group = tree,
                   colour = extinction))+
  scale_color_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  geom_point(size=3, shape=21, fill="white")+
  geom_line(alpha = 0.3)+
  geom_point(data = ages.mean.incorrect,
             size=3, shape = 21, fill="white",
             alpha = 0.4, color = "#969696")+
  geom_line(data = ages.mean.incorrect,
            alpha = 0.4, color = "#969696")+
  geom_text(data = ages.comparison.mean, aes(x = 2.05, y = 1.7,
                                             label = paste0("Error rate:", " ",
                                                            error,"%")),
            inherit.aes = FALSE,
            size = 3, 
            family = "serif")+
  facet_wrap(~extinction,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  ylab("Corrected age (log)")+
  xlab(NULL)+
  mynamestheme+
  theme(legend.position = "none")



###png object
png("text/figures/Figure.c_status.png", 
    width = 17, height = 17,
    units = "cm", 
    pointsize = 8, res = 300)

##bottom
bottom <- text_grob("Conservation status",
                    family = "serif", size = 13,  face = "bold")

##grid for plotting both figures
grid.arrange(ftrue, fphylo, fprob, nrow = 3, bottom= bottom)

dev.off()

