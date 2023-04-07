##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))



# Testing the stochastic function -----------------------------------------

###using the same data generated in the simulation of the conservation status

load(file = file.path(getwd(), pro, c_status, "ages.high.stochastic.RData"))

load(file = file.path(getwd(), pro, c_status, "ages.int.stochastic.RData"))

load(file = file.path(getwd(), pro, c_status, "ages.low.stochastic.RData"))

ages.total <- read_csv(file = file.path(getwd(), pro, c_status,
                                        "ages.total.csv"))

load(file = file.path(getwd(), pro, c_status, "trees.conservation.RData"))


##high extinction

##assigning tree number 
for(i in seq_along(ages.high.stochastic)){
  ages.high.stochastic[[i]]$tree <- rep(paste0("tree.", i),
                                       nrow(ages.high.stochastic[[i]]))
}

##unlist
ages.high.st <- do.call("rbind", ages.high.stochastic) 

ages.high.st$species <- paste0(ages.high.st$species,".", ages.high.st$tree, ".",
                               "high")

ages.high.st <- select(ages.high.st, species, dist_ages)


##intermediate extinction

##assigning tree number 
for(i in seq_along(ages.int.stochastic)){
  ages.int.stochastic[[i]]$tree <- rep(paste0("tree.", i),
                                        nrow(ages.int.stochastic[[i]]))
}

##unlist
ages.int.st <- do.call("rbind", ages.int.stochastic) 

ages.int.st$species <- paste0(ages.int.st$species,".", ages.int.st$tree, ".",
                              "intermediate")

ages.int.st <- select(ages.int.st, species, dist_ages)

##low extinction

##assigning tree number 
for(i in seq_along(ages.low.stochastic)){
  ages.low.stochastic[[i]]$tree <- rep(paste0("tree.", i),
                                       nrow(ages.low.stochastic[[i]]))
}

##unlist
ages.low.st <- do.call("rbind", ages.low.stochastic) 

ages.low.st$species <- paste0(ages.low.st$species,".", ages.low.st$tree, ".",
                              "low")

ages.low.st <- select(ages.low.st, species, dist_ages)

##rbinding
ages.df <- rbind(ages.high.st, ages.int.st, ages.low.st)

##changing species name for joining dfs
ages.total$species <- paste0(ages.total$species, '.', ages.total$extinction)

##ages.total.join
ages.total.join <- left_join(ages.total, ages.df, by = "species")

ages.total.join$extinction <- as.factor(ages.total.join$extinction)

ages.total.join$extinction <- factor(ages.total.join$extinction,
                                     levels=c("low", "intermediate", "high"))

##saving
save(ages.total.join, file = file.path(getwd(), pro, st_sp,
                                       "ages.total.join.RData"))

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

##accuracy of the highest prob age regarding true age

##determining how many species are 100% accurate regarding the stochastic fun
accurate.st <- ages.total.join %>%
  mutate(abs.accuracy.st = abs(rStochastic - rTrue.age)) %>% 
  filter(abs.accuracy.st < 1e-8) %>% group_by(extinction) %>%
  count()



png("text/figures/Figure6.accuracy.mode.stochastic.png",
    width = 10, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f6 <- ages.total.join %>%
  mutate(abs.accuracy = abs(rStochastic - rTrue.age)) %>% 
  filter(abs.accuracy > 1e-8) %>% 
  ggplot(aes(x = rStochastic - rTrue.age))+
  geom_histogram(position = "dodge", breaks = seq(-0.8, 0.4, 0.02))+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = n),
               data = accurate.st,
               colour = "red")+
  geom_point(aes(x = 0, y = n), data = accurate.st, colour = "red")+
  geom_text(data = accurate.st, aes(x = 0.35, y = 4000,
                                    label = paste("Correct estimation:",
                                                  round(n/100),"%")),
            size = 3.3, family = "serif")+
  xlim(-0.75,0.6)+
  #ylim(0, 5000)+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  xlab("Highest probability age - True age")+
  ylab("Species count")+
  ggtitle("Stochastic function accuracy")

f6 + mynamestheme

dev.off()

##accuracy of the mean prob age regarding true age

##creating column mean.age
ages.total.join$mean.age <- rep(NA, nrow(ages.total.join))

####calculating the mean of the prob ages
list.ages <- vector("list", length = nrow(ages.total.join))

for(i in 1:nrow(ages.total.join)){
  list.ages[[i]] <- ages.total.join[i,]$dist_ages
  list.ages[[i]] <- as.tibble(list.ages[[i]][[1]])
  list.ages[[i]]$Age <- as.numeric(as.character(list.ages[[i]]$Age))
  list.ages[[i]]$Probability <- round(list.ages[[i]]$Probability*100)
  ages.total.join[i,]$mean.age <-  mean(rep(list.ages[[i]]$Age,
                                              list.ages[[i]]$Probability))
  
}

ages.total.join <- ages.total.join %>% mutate(rMean = mean.age/root.age)



##determining how many species are 100% accurate regarding the stochastic fun


accurate.mean <- ages.total.join %>%
  mutate(abs.accuracy.mean = abs(rMean - rTrue.age)) %>% 
  filter(abs.accuracy.mean < 1e-8) %>% 
  group_by(extinction) %>%
  count() ##there are no correct estimations



png("text/figures/Figure6.1.accuracy.mean.stochastic.png",
    width = 10, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f6.1 <- ages.total.join %>%
  mutate(abs.accuracy.mean = abs(rMean - rTrue.age)) %>% 
  ggplot(aes(x = rMean - rTrue.age))+
  geom_histogram(position = "dodge", breaks = seq(-1, 1, 0.02))+
  geom_text(x = 0.4, y = 1500,
              label = "No correct estimation",
            size = 3.3, family = "serif")+
  xlim(-0.75,0.75)+
  facet_wrap(~ extinction, nrow = 3,
             labeller = as_labeller(c(low = "Low extinction",
                                      intermediate = "Intermediate extinction",
                                      high = "High extinction")))+
  theme_bw()+
  xlab("Mean probability age - True age")+
  ylab("Species count")+
  ggtitle("Stochastic function accuracy")

f6.1 + mynamestheme

dev.off()

######################estimating coverage using the Confidence intervals

coverage.st <- ages.total.join %>%
                  filter(abs(True.age - max.prob.age) < 1e-08 | 
                           True.age >= CI_low & True.age <= CI_up) %>% 
                 count(extinction) %>% 
                 mutate(cov = n/10000 * 100) 

png("text/figures/Figure6.2.coverage.CI.stochastic.png",
    width = 10, height = 12, units = "cm", 
    pointsize = 8, res = 300)


coverage.st %>% ggplot(aes(extinction, cov, fill = extinction))+
      geom_col(width=0.75)+
      geom_text(aes(label = paste(round(cov), "%")),
                size = 3.5, family = "serif",
                position = position_stack(vjust = 1.05))+
      ylim(0, 100) + 
      scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
      scale_x_discrete(labels = c("low" = "Low","intermediate" = "Intermediate",
                                  "high" = "High"))+
      ylab("Coverage (%)") + 
      xlab("Extinction Background")+
      theme_bw()+
      theme(legend.position = "none")+
      mynamestheme
        
     
dev.off()






