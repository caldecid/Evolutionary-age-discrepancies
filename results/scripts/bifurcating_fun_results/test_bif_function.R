########testing the bifurcating function

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

###creating new folder
path <- "results/data/processed"
newdir <- "bif_function"
dir.create(file.path(path, newdir))


##loading trees from the conservation simulation with high, intermediate and low
##extinction rates
load(file = file.path(getwd(), pro, c_status,"trees.conservation.RData"))


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

#############simulating taxonomies only bifurcating######################

##high extinction
tax.high <- lapply(trees.high, sim.taxonomy, beta = 1)

##intermediate extinction
tax.int <- lapply(trees.int, sim.taxonomy, beta = 1)

##low extinction
tax.low <- lapply(trees.low, sim.taxonomy, beta = 1)

save(tax.high, tax.int, tax.low,
     file = file.path(getwd(), pro, newdir, "tax_bif.RData"))

load(file = file.path(getwd(), pro, newdir, "tax_bif.RData"))


# Get ages and variables --------------------------------------------------


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

##lambda
ages.high$lambda <- rep(0.3, each = 100)

##mu
ages.high$mu <- rep(0.25, each = 100)

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

##lambda
ages.int$lambda <- rep(0.3, each = 100)

##mu
ages.int$mu <- rep(0.15, each = 100)

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

##lambda
ages.low$lambda <- rep(0.3, each = 100)

##mu
ages.low$mu <- rep(0.05, each = 100)

##scale ages
ages.low <- ages.low %>% mutate(rTrue.age = True.age/root.age,
                                rPhylo.age = Estimated.age/root.age)


##extinction
ages.low$extinction <- rep("low", nrow(ages.low))

##species name
ages.low$species <- paste0(ages.low$label,".", ages.low$tree)


# bifurcating function ----------------------------------------------------

#empty list
high.list <- list()

for(i in 1:nrow(ages.high)){
  high.list[[i]] <- get_sp_age_prob(lam = ages.high$lambda[i],
                                   mu = ages.high$mu[i],
                                   node_age = ages.high$Estimated.age[i])
}


##expected age
high.expected.age <- as.vector(do.call("rbind",
                     lapply(high.list, function(x) sum(x$time * x$prob))))

##median age
high.median.age <- as.vector(do.call("rbind", lapply(high.list, function(x)
                                 weighted.median(x = x$time,
                                   w = x$prob))))

##confidence interval

#low
lwr.high <- sapply(high.list, function(x) weightedQuantile(x$time, x$prob,
                                                             0.025))
#upper
upr.high <- sapply(high.list, function(x) weightedQuantile(x$time, x$prob,
                                                           0.975))

##joining dataframes
ages.high.bif <- cbind(ages.high, high.expected.age, high.median.age,
                             lwr.high, upr.high) %>% 
                    rename(ex_age = high.expected.age,
                           median_age = high.median.age,
                           lwr.ci = lwr.high,
                           upr.ci = upr.high)
  
#########intermediate extinction scenario
int.list <- list()

for(i in 1:nrow(ages.int)){
  int.list[[i]] <- get_sp_age_prob(lam = ages.int$lambda[i],
                                   mu = ages.int$mu[i],
                                   node_age = ages.int$Estimated.age[i])
}
##expected age
int.expected.age <- sapply(int.list, function(x) sum(x$time * x$prob))

##median age
int.median.age <- sapply(int.list, function(x)
                              weighted.median(x = x$time,
                                            w = x$prob))

##confidence interval

#low
lwr.int <- sapply(int.list, function(x) weightedQuantile(x$time, x$prob,
                                                           0.025))
#upper
upr.int <- sapply(int.list, function(x) weightedQuantile(x$time, x$prob,
                                                           0.975))
##joining dataframes
ages.int.bif <- cbind(ages.int, int.expected.age, int.median.age,
                       lwr.int, upr.int) %>% 
                    rename(ex_age = int.expected.age,
                           median_age = int.median.age,
                           lwr.ci = lwr.int,
                           upr.ci = upr.int)


#######low extinction scenario
low.list <- list()

for(i in 1:nrow(ages.low)){
  low.list[[i]] <- get_sp_age_prob(lam = ages.low$lambda[i],
                                   mu = ages.low$mu[i],
                                   node_age = ages.low$Estimated.age[i])
}

##expected age
low.expected.age <- sapply(low.list, function(x) sum(x$time * x$prob))

##median age
low.median.age <- sapply(low.list, function(x)
                           weighted.median(x = x$time,
                                w = x$prob))

##confidence interval

#low
lwr.low <- sapply(low.list, function(x) weightedQuantile(x$time, x$prob,
                                                         0.025))
#upper
upr.low <- sapply(low.list, function(x) weightedQuantile(x$time, x$prob,
                                                         0.975))
##joining dataframes
ages.low.bif <- cbind(ages.low, low.expected.age, low.median.age,
                      lwr.low, upr.low) %>% 
                rename(ex_age = low.expected.age,
                       median_age = low.median.age,
                       lwr.ci = lwr.low,
                       upr.ci = upr.low)

##################merging the dataframes with different extinction#####
ages.total.bif <- rbind(ages.high.bif, ages.int.bif, ages.low.bif) %>% 
                  mutate(rexp_age = ex_age/root.age,
                         rmedian_age = median_age/root.age,
                         rlow_ci = lwr.ci/root.age,
                         rupr_ci = upr.ci/root.age)

ages.total.bif$extinction <- as.factor(ages.total.bif$extinction)

ages.total.bif$extinction <- factor(ages.total.bif$extinction,
                                    levels = c("low", "intermediate", "high"),
                                    ordered = TRUE)

##saving
write_csv(ages.total.bif, file = file.path(pro, newdir, "ages.total.bif.csv"))

ages.total.bif <- read_csv(file = file.path(pro, newdir, "ages.total.bif.csv"))

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
                      legend.background = element_rect(fill = "gray90",
                                             size = 0.5, linetype = "dotted"),
                      legend.position = "bottom")



######################estimating function coverage##################

#############coverage for a specific tree
##pivot longer for having in the same column phylo and mean age

ages.filt <- ages.total.bif %>% filter(tree == "tree.13") %>% 
  pivot_longer(cols = c("rPhylo.age", "rmedian_age"),
               names_to = "Estimated",
               values_to = "ages")

ages.filt$Estimated <- as.factor(ages.filt$Estimated)

ages.filt$Estimated <- factor(ages.filt$Estimated, 
                              levels = c("rPhylo.age", "rmedian_age"))

#####ploting 

cover.specific <- ggplot(ages.filt,
                         aes(x = rTrue.age, y = ages, color = Estimated))+
  geom_segment(aes(x = rTrue.age, y = rlow_ci, xend = rTrue.age,
                   yend = rupr_ci),
               colour = "gray", alpha = 0.5)+
  geom_point(alpha = 0.7)+
  scale_color_discrete(labels = c(rmedian_age = "Median",
                                  rPhylo.age = "Phylogenetic"))+
  labs(color = NULL)+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  
  facet_wrap(~extinction, nrow = 3,
             labeller = as_labeller(c(high = "High extinction",
                                      intermediate = "Intermediate extinction",
                                      low = "Low extinction")))+
  theme_bw()+
  
  xlab("True age")+
  ylab("Estimated age")+
  mynamestheme



dev.off()

#######total coverage##############
coverage.bif <- ages.total.bif%>%
  filter(round(True.age, digits = 4) >= round(lwr.ci, digits = 4) &
           round(True.age, digits = 4) <= round(upr.ci, digits = 4))%>% 
  count(extinction) %>% 
  mutate(cov = n/1000) 


cover.total <- ggplot(coverage.bif, aes(extinction, cov, fill = extinction))+
  geom_col(width=0.75)+
  geom_text(aes(label = paste(round(cov), "%")),
            size = 3.5, family = "serif",
            position = position_stack(vjust = 0.98))+
  ylim(0, 100) + 
  scale_fill_manual(values = c("#7fc97f","#beaed4","#fdc086"))+
  scale_x_discrete(labels = c("low" = "Low","intermediate" = "Intermediate",
                              "high" = "High"))+
  ylab("Coverage (%)") + 
  xlab("Extinction Background")+
  theme_bw()+
  mynamestheme+
  theme(legend.position = "none")

######binding both figures
png("text/figures/FigureSM3.coverage.png", 
    width = 15, height = 15, units = "cm", 
    pointsize = 8, res = 300)


ggarrange(cover.total, cover.specific, ncol = 2,
          common.legend = FALSE)

dev.off()



####################MAPE###############


##For the MAPE for each category
ages.average <- ages.total.bif %>% 
               slice(sample(1:n())) %>% 
            group_by(extinction, tree) %>% 
  summarise(mape.phy = mean(abs(True.age - Estimated.age)/True.age)*100,
            mape.mean = mean(abs(True.age - ex_age)/True.age)*100,
            mape.median = mean(abs(True.age - median_age)/True.age)*100) %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.average$estimate <- factor(ages.average$estimate,
                            levels = c("mape.phy",
                                       "mape.ex",
                                       "mape.median"))

##for the delta MAPE
ages.delta <- ages.total.bif %>% 
      slice(sample(1:n())) %>% 
          group_by(extinction, tree) %>% 
          summarise(mape.phy = mean(abs(True.age - Estimated.age)/True.age)*100,
                    mape.mean = mean(abs(True.age - ex_age)/True.age)*100,
                    mape.median = mean(abs(True.age - median_age)/True.age)*100,
                    delta.mean = mape.mean - mape.phy,
                    delta.median =  mape.median - mape.phy) %>% 
          select(-c(mape.phy, mape.mean, mape.median)) %>% 
          pivot_longer(cols = starts_with("delta."),
                       values_to = "delta", names_to = "estimate")

##Summaries for results
delta.summary <- ages.delta %>% 
                   pivot_wider(names_from = estimate, values_from = delta) %>%  
                   group_by(extinction) %>% 
                   summarise(mean.delta.mean = mean(delta.mean),
                           sd.delta.mean = sd(delta.mean),
                           mean.delta.median = mean(delta.median),
                           sd.delta.median = sd(delta.median))
###MAPE
png("text/figures/MAPE.bif.prob.png", 
    width = 15, height = 15, units = "cm", 
    pointsize = 8, res = 300)


  ggplot(ages.average, aes(y = mape, x = estimate, fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~extinction,  labeller = as_labeller(c(high = "High extinction",
                                  intermediate = "Intermediate extinction",
                                  low = "Low extinction")))+
  scale_fill_discrete(name = "Estimation",
                        labels = c("Phylogenetic","Mean", "Median"))+
  ylim(0, 150)+
  ylab("MAPE (%)")+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())
  
  dev.off()

  
###DELTA MAPE
 png("text/figures/Delta.MAPE.png", 
      width = 15, height = 15, units = "cm", 
      pointsize = 8, res = 300)
 
 ggplot(ages.delta, aes(y = delta, x = estimate, fill = estimate))+
    geom_boxplot(outlier.shape = NA)+
    facet_wrap(~extinction,  labeller = as_labeller(c(high = "High extinction",
                                      intermediate = "Intermediate extinction",
                                                      low = "Low extinction")))+
    scale_fill_discrete(name = "Estimation",
                        labels = c("Mean", "Median"))+
    geom_hline(yintercept = 0, linetype = "dashed",
               size = 1, colour = "red")+
    ylim(-80, 10)+
    ylab(expression(bold(Delta* "MAPE (%)")))+
    xlab(NULL)+
    theme_bw()+
    mynamestheme+
    theme(axis.text.x = element_blank())
  
  dev.off()
  