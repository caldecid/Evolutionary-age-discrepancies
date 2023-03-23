##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


##future package conditions 
plan("multisession", workers = 4)

# testing Torsten function for high turnover -------------------------------------


##Stochastic function (Torsten) for the trees dataset
##high
ages.torst.high <- future_lapply(tree.extant.high, stochastic.sp.age.torsten)

##assigning tree numbers
for(i in seq_along(ages.torst.high)){
  ages.torst.high[[i]]$tree <- paste0("tree.", i)
}

##colapsing list into a dataframe
ages.tor.high.df <- do.call("rbind", ages.torst.high)


##assigning the tree to the species name
ages.tor.high.df$species <- paste0(ages.tor.high.df$species, ".",
                                  ages.tor.high.df$tree)






##colapsing list into a dataframe and renaming variables
taxa.high <- do.call("rbind", taxa.true.high) %>% rename(True.age = Age,
                                          Estimated.age = tip_length_reconstr,
                                           N_sisters = n_sisters_reconstr) %>% 
                                      select(label, True.age, Estimated.age)

taxa.high$label <- paste0(taxa.high$label, ".tree.",
                          rep(1:length(ages.torst.high), each =100))


##root age
root.age.high <- sapply(tree.extant.high, tree.max)


##adding the true age of species and extracting max prob age
ages.tor.high.df <- left_join(ages.tor.high.df, taxa.high,
                              by = c("species" = "label")) %>%
                unnest(dist_ages) %>% 
                filter(Probability == max(Probability)) %>%
                select(species, tree, True.age, Estimated.age, Probability, Age) %>% 
                rename(max.prob.age = Age) 




##organizing the data
ages.tor.high.df$max.prob.age <-
  as.numeric(as.character(ages.tor.high.df$max.prob.age))

##extinction scenario
ages.tor.high.df$extinction <- rep("high", nrow(ages.tor.high.df))

##assigning root age
ages.tor.high.df$root.age <- rep(root.age.high, each = 100)

##scaling ages to root age
ages.tor.high.df <- ages.tor.high.df %>% mutate(rTrue.age = True.age/root.age,
                                              rmax.age = max.prob.age/root.age)

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

##intermediate
ages.torst.int <- future_lapply(tree.extant.int, stochastic.sp.age.torsten)

##assigning tree numbers
for(i in seq_along(ages.torst.int)){
  ages.torst.int[[i]]$tree <- paste0("tree.", i)
}

##collapsing list into df
ages.tor.int.df <- do.call("rbind", ages.torst.int)

##assigning the tree to the species name
ages.tor.int.df$species <- paste0(ages.tor.int.df$species, ".",
                                  ages.tor.int.df$tree)

##adding the true age of species
ages.tor.int.df <- left_join(ages.tor.int.df, taxa.int,
                                                 by = c("species"="label"))%>% 
                                           unnest(dist_ages) %>% 
                                  filter(Probability == max(Probability)) %>%
                    select(species, tree, True.age,
                           Estimated.age, Probability, Age) %>% 
                    rename(max.prob.age = Age)

##transforming and adding data 
ages.tor.int.df$max.prob.age <-
  as.numeric(as.character(ages.tor.int.df$max.prob.age))

ages.tor.int.df$extinction <- rep("intermediate", nrow(ages.tor.int.df))

## root age and scaling
root.age.inter <- tibble(root.age.int, tree = paste0("tree.", 1:100))

ages.tor.int.df <- left_join(ages.tor.int.df, root.age.inter, by = "tree") %>% 
                   rename(root.age = root.age.int)

###error
ages.tor.int.df <- ages.tor.int.df %>% mutate(rTrue.age = True.age/root.age,
                                              rmax.age = max.prob.age/root.age)

# low extinction -----------------------------------------------------------

tree.0.low <- sim.bd.taxa(n = 100, numbsim = 100, lambda = 0.2, mu =0.01)

####taxonomy simulation only budding 

taxa.0.low <- lapply(tree.0.low, sim.taxonomy, beta = 0)

###obtainig true ages

taxa.true.low <- Map(getAgesExtantSpecies, taxa.0.low, tree.0.low) 




##colapsing list to dataframe
taxa.low <- do.call("rbind", taxa.true.low) %>% rename(True.age = Age,
                                    Estimated.age = tip_length_reconstr,
                                   N_sisters = n_sisters_reconstr) %>% 
  select(label, True.age, Estimated.age)

##assigning trees to the label
taxa.low$label <- paste0(taxa.low$label, ".tree.",
                         rep(1:length(taxa.true.low), each =100))



##prunning the simulated tree
tree.extant.low <- lapply(tree.0.low, prune.fossil.tips)

##root.age
root.age.low <- sapply(tree.extant.low, tree.max)


##stochastic function using tree extant as input

ages.torst.low <- future_lapply(tree.extant.low, stochastic.sp.age.torsten)

##assigning tree numbers
for(i in seq_along(ages.torst.low)){
  ages.torst.low[[i]]$tree <- paste0("tree.", i)
}

##collapsing list into df
ages.tor.low.df <- do.call("rbind", ages.torst.low)

##assigning the tree to the species name
ages.tor.low.df$species <- paste0(ages.tor.low.df$species, ".",
                                  ages.tor.low.df$tree)

##adding the true age of species
ages.tor.low.df <- left_join(ages.tor.low.df, taxa.low,
                             by = c("species"="label"))%>% 
                  unnest(dist_ages) %>% 
                  filter(Probability == max(Probability)) %>%
                  select(species, tree, True.age,
                         Estimated.age, Probability, Age) %>% 
                          rename(max.prob.age = Age)

##arranging and assigning values
ages.tor.low.df$max.prob.age <-
  as.numeric(as.character(ages.tor.low.df$max.prob.age))

##extinction scenario
ages.tor.low.df$extinction <- rep("low", nrow(ages.tor.low.df))

##root age
ages.tor.low.df$root.age <- rep(root.age.low, each = 100)

##calculating scaled ages
ages.tor.low.df <- ages.tor.low.df %>% 
                                     mutate(rTrue.age = True.age/root.age,
                                            rmax.age = max.prob.age/root.age)



# saving results ------------------------------------------------------------

##saving the nested dataframes
save(ages.torst.high, ages.torst.int, ages.torst.low,
     file = file.path(getwd(), pro, st_sp, "ages.torst.RData"))

##building dataframes with max probability ages
ages.tor.max <- rbind(ages.tor.high.df, ages.tor.int.df, ages.tor.low.df)

##calculating accuracy for stochastic
ages.tor.max <- ages.tor.max %>%
  mutate(accuracy.sto = rmax.age - rTrue.age)

##saving
write_csv(ages.tor.max, file = file.path(getwd(), pro, st_sp,
                                          "ages.torst.max.csv"))


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
                      legend.background = element_rect(fill="gray90",
                                                    size=.5, linetype="dotted"))


##scale True age vs max prob age
ages.tor.max %>%
  ggplot(aes(rTrue.age, rmax.age))+
  geom_point(alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  facet_wrap(~ extinction, nrow = 3)+  
  xlab("scale True age")+
  ylab("scale highest probability age")
 

##accuracy of the max prob age regarding true age

##determining how many species are 100% accurate regarding the stochastic fun
accurate.st <- ages.tor.max %>%
  mutate(abs.accuracy.st = abs(accuracy.sto)) %>% 
  filter(abs.accuracy.st < 1e-8) %>% group_by(extinction) %>%
  count()



png("Torsten.function.png", width = 10, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f6 <- ages.tor.max %>%
  mutate(abs.accuracy = abs(accuracy.sto)) %>% 
  filter(abs.accuracy > 1e-8) %>% 
  ggplot(aes(x = accuracy.sto))+
  geom_histogram(position = "dodge", breaks = seq(-0.8, 0.4, 0.02))+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = n),
               data = accurate.st,
               colour = "red")+
  geom_point(aes(x = 0, y = n), data = accurate.st, colour = "red")+
  geom_text(data = accurate.st, aes(x = 0.35, y = 4000,
            label = paste("Correct estimation:", round(n/100),"%")),
            size = 3.3, family = "serif")+
  xlim(-0.75,0.6)+
  #ylim(0, 5000)+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(high = "High extinction",
                                      intermediate = "Intermediate extinction",
                                      low = "Low extinction")))+
  theme_bw()+
  xlab("Highest probability age - True age")+
  ylab("Species count")+
  ggtitle("Torsten function accuracy")

f6 + mynamestheme

dev.off()

