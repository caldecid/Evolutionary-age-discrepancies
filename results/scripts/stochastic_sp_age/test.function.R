##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


##future package conditions 
plan("multisession", workers = 4)

# testing the function for high turnover -------------------------------------

#########simulating 100 trees with high turnover

tree.0.high <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.2, mu =0.19)

####taxonomy simulation only budding for obtaining true ages

taxa.0.high <- lapply(tree.0.high, sim.taxonomy, beta = 0)

###obtaining true ages
taxa.true.high <- Map(getAgesExtantSpecies, taxa.0.high, tree.0.high) 
                  
##prunning the simulated tree
tree.extant.high <- lapply(tree.0.high, prune.fossil.tips)


##Stochastic function for the trees dataset
ages.test.high <- vector("list", length = 100)

for(i in seq_along(ages.test.high)){
  ages.test.high[[i]] <- try(stochastic.sp.age(tree.extant.high[[i]]))
}

##assigning tree numbers
for(i in seq_along(ages.test.high)){
  ages.test.high[[i]]$tree <- paste0("tree.", i)
}

##filtering elements with error which became a list due to tree assignation
ages.test.high.1 <- ages.test.high[lapply(ages.test.high, class) != "list"]

##colapsing list into a dataframe
ages.high <- do.call("rbind", ages.test.high.1)

##assigning the tree to the species name
ages.high$species <- paste0(ages.high$species, ".", ages.high$tree)


##taxonomies
taxa.true.high.1 <- taxa.true.high[lapply(ages.test.high, class) != "list"]


##colapsing list into a dataframe and renaming variables
taxa.high <- do.call("rbind", taxa.true.high.1) %>% rename(True.age = Age,
                                            Estimated.age = tip_length_reconstr,
                                            N_sisters = n_sisters_reconstr) %>% 
                                select(label, True.age, Estimated.age)

taxa.high$label <- paste0(taxa.high$label, ".", ages.high$tree)


##adding the true age of species and unnesting for obtaining the high prob age
ages.high.1 <- left_join(ages.high, taxa.high, by = c("species" = "label")) %>% 
                unnest(dist_ages) %>% 
                filter(Probability == max(Probability)) %>%
               select(species, tree, True.age, Estimated.age, Probability, Age) %>% 
               rename(max.prob.age = Age) 



##organizing the data
ages.high.1$max.prob.age <-
                         as.numeric(as.character(ages.high.1$max.prob.age))

##extinction scenario
ages.high.1$extinction <- rep("high", nrow(ages.high.1))

##assigning root age
##root age
root.age.high <- sapply(tree.extant.high, tree.max)
##tibble
root.high <- tibble(root.age.high, tree = paste0("tree.", 1:100))


##joining dataframes and scaling ages
ages.high.1 <- left_join(ages.high.1, root.high, by = "tree") %>% 
                  mutate(rTrue.age = True.age/root.age.high,
                      rmax.age = max.prob.age/root.age.high) %>% 
                  rename(root.age = root.age.high)

# Intermediate turnover ---------------------------------------------------

########simulating a tree with int turnover

tree.0.int <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.2, mu =0.1)

####taxonomy simulation only budding 

taxa.0.int <- lapply(tree.0.int, sim.taxonomy, beta = 0)

###obtainig true ages
taxa.true.int <- Map(getAgesExtantSpecies, taxa.0.int, tree.0.int) 

##colapsing list into a dataframe and renaming variables
taxa.int <- do.call("rbind", taxa.true.int) %>% rename(True.age = Age,
                                           Estimated.age = tip_length_reconstr,
                                          N_sisters = n_sisters_reconstr) %>% 
                                    select(label, True.age, Estimated.age)

##assigning trees to the label
taxa.int$label <- paste0(taxa.int$label, ".tree.",
                          rep(1:length(taxa.0.int), each =100))

 
##prunning the simulated tree
tree.extant.int <- lapply(tree.0.int, prune.fossil.tips)

##root age
root.age.int <- sapply(tree.extant.int, tree.max)

##stochastic function using tree extant as input
ages.test.int <- lapply(tree.extant.int, stochastic.sp.age)

##assigning tree numbers
for(i in seq_along(ages.test.int)){
  ages.test.int[[i]]$tree <- paste0("tree.", i)
}

##collapsing list into df
ages.int <- do.call("rbind", ages.test.int)

##assigning the tree to the species name
ages.int$species <- paste0(ages.int$species,".", ages.int$tree)

##adding the true age of species
ages.int.1 <- left_join(ages.int, taxa.int, by = c("species"="label")) %>%
   unnest(dist_ages) %>% 
   filter(Probability == max(Probability)) %>%
   select(species, tree, True.age, Estimated.age, Probability, Age) %>% 
   rename(max.prob.age = Age) 

##assigning variables
ages.int.1$max.prob.age <-
  as.numeric(as.character(ages.int.1$max.prob.age))

ages.int.1$extinction <- rep("intermediate", nrow(ages.int.1))

##root age
root.int <- tibble(root.age = root.age.int, tree = paste0("tree.", 1:100))
##merging data frames and scaling ages
ages.int.1 <- ages.int.1 %>% left_join(root.int, by = "tree") %>% 
                           mutate(rTrue.age = True.age/root.age,
                                  rmax.age = max.prob.age/root.age)

# low extinction -----------------------------------------------------------

tree.0.low <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.2, mu =0.01)

####taxonomy simulation only budding 

taxa.0.low <- lapply(tree.0.low, sim.taxonomy, beta = 0)

###obtainig true ages

taxa.true.low <- Map(getAgesExtantSpecies, taxa.0.low, tree.0.low) 

##filtering list elements with error
taxa.true.low.1 <- taxa.true.low[lapply(ages.test.low, class) != "try-error"]


##colapsing list to dataframe
taxa.low <- do.call("rbind", taxa.true.low.1) %>% rename(True.age = Age,
                                            Estimated.age = tip_length_reconstr,
                                           N_sisters = n_sisters_reconstr) %>% 
                                        select(label, True.age, Estimated.age)

##assigning trees to the label
taxa.low$label <- paste0(taxa.low$label, ".tree.",
                         rep(1:length(taxa.true.low.1), each =100))


  
##prunning the simulated tree
tree.extant.low <- lapply(tree.0.low, prune.fossil.tips)

root.age.low <- sapply(tree.extant.low, tree.max)


##stochastic function using tree extant as input

ages.test.low <- vector("list", length = 100)
for(i in seq_along(ages.test.low)){
  ages.test.low[[i]] <- try(stochastic.sp.age(tree.extant.low[[i]]))
}

##assigning tree numbers
for(i in seq_along(ages.test.low)){
  ages.test.low[[i]]$tree <- paste0("tree.", i)
}

##filtering elements with error
ages.test.low.1 <- ages.test.low[lapply(ages.test.low, class) != "list"]

##colapsing list into a dataframe
ages.low <- do.call("rbind", ages.test.low.1) 

##assigning tree number to species
ages.low$species <- paste0(ages.low$species,".", ages.low$tree)


##adding the true age of species
ages.low.1 <- left_join(ages.low, taxa.low, by = c("species" = "label")) %>% 
               unnest(dist_ages) %>% 
              filter(Probability == max(Probability)) %>%
            select(species, tree, True.age, Estimated.age, Probability, Age) %>% 
            rename(max.prob.age = Age) 

##organizing ages 
ages.low.1$max.prob.age <-
  as.numeric(as.character(ages.low.1$max.prob.age))

##extinction intensity
ages.low.1$extinction <- rep("low", nrow(ages.low.1))

##root ages
root.low <- tibble(root.age = root.age.low, tree = paste0("tree.", 1:100))

##merging dataframes and calculating scale ages
ages.low.1 <- ages.low.1 %>% left_join(root.low, by = "tree") %>% 
              mutate(rTrue.age = True.age/root.age,
                     rmax.age = max.prob.age/root.age)

# saving results ------------------------------------------------------------

##saving the nested dataframes
save(ages.high, ages.int, ages.low,
     file = file.path(getwd(), pro, st_sp, "ages.test.RData"))

##building dataframes with max probability ages
ages.test.max <- rbind(ages.high.1, ages.int.1, ages.low.1)

##calculating accuracy for stochastic and phylogenetic scale calculation
ages.test.max <- ages.test.max %>%
                 mutate(accuracy.sto = rmax.age - rTrue.age)
                       

##saving
write_csv(ages.test.max, file = file.path(getwd(), pro, st_sp,
                                          "ages.test.max.csv"))


# figures -----------------------------------------------------------------

##True age vs max prob age
ages.test.max %>% ggplot(aes(rTrue.age, rmax.age, extinction))+
  geom_point(alpha = 0.7, aes(color = extinction))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  xlab("True age (scaled)")+
  ylab("Highest age probability (scaled)")

##accuracy of the max prob age regarding true age

##determining how many species are 100% accurate regarding the stochastic fun
accurate.st <- ages.test.max %>%
  mutate(abs.accuracy.st = abs(accuracy.sto)) %>% 
  filter(abs.accuracy.st < 1e-8) %>% group_by(extinction) %>%
  count()



png("Figure6.sto.function.png", width = 10, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f6 <- ages.test.max %>%
  mutate(abs.accuracy = abs(accuracy.sto)) %>% 
  filter(abs.accuracy > 1e-8) %>% 
  ggplot(aes(x = accuracy.sto))+
  geom_histogram(position = "dodge", breaks = seq(0, 0.5, 0.02))+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = n),
               data = accurate.st,
               colour = "red")+
  geom_point(aes(x = 0, y = n), data = accurate.st, colour = "red")+
  geom_text(data = accurate.st, aes(x = 0.2, y = 7500,
            label = paste("Correct estimation:", round(n/100),"%")),
            size = 3.4, family = "serif")+
  xlim(-0.10, 0.5)+
  #ylim(0, 5000)+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(high = "High extinction",
                                      intermediate = "Intermediate extinction",
                                      low = "Low extinction")))+
  theme_bw()+
  xlab("Highest probability age - True age (scale)")+
  ylab("Species count")+
  ggtitle("Stochastic ages function accuracy")
  
f6 + mynamestheme

dev.off()

# confidence intervals ----------------------------------------------------


x1 <- ages.test.low.1[[35]][[2]][[1]]

x1$freq <- x1$Probability*1000



x1.vec <- as.numeric(as.character(rep(x1$Age, x1$freq)))

hist(x1.vec)

##for all the list

##low extinction
freq.ages.low <- vector("list", length(ages.test.low.1))

for(i in 1:length(ages.test.low.1)){
  freq.ages.low[[i]] <- as.numeric(as.character(rep(ages.test.low.1[[i]][[2]][[1]]$Age,
                                                    ages.test.low.1[[i]][[2]][[1]]$Probability*1000)))
}

freq.ages.low.un <- unlist(freq.ages.low)

####confidence intervals

CI95 = hdi(freq.ages.low.un, 0.95)

hist(freq.ages.low.un)

##high extinction
freq.ages.high <- vector("list", length(ages.test.high.1))

for(i in 1:length(ages.test.high.1)){
  freq.ages.high[[i]] <- as.numeric(as.character(rep(ages.test.high.1[[i]][[2]][[1]]$Age,
                            ages.test.high.1[[i]][[2]][[1]]$Probability*1000)))
}

freq.ages.high.un <- unlist(freq.ages.high)

hist(freq.ages.high.un)

######Torsten function

ages.torst <- stochastic.sp.age.torsten(tree.0)

ages.mine <- stochastic.sp.age(tree.0)
