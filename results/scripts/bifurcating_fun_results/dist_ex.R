# Simulating experiment across extinction and sampling scenarios ----------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

###loops for generating the simualted data

##different extinction rates
mu.s <- c(0.05, 0.15, 0.25)
##different sampling rates
rho.s <- c(1, 0.75, 0.50)
##fix parameters
lambda <- 0.3
branching_time <- 3
reps <- 100000

##list for keeping datasets
##data general for
data.general <- vector("list", length = 3)

for(i in 1:3){
  
  data.general[[i]] <- vector("list", length = 3)
}



##dent general for density curve
dent.list <- vector("list", length = 3)

for(i in 1:3){
  dent.list[[i]] <- vector("list", length = 3)
}

for(i in seq_along(mu.s)){
  
  for(j in seq_along(rho.s)) {
    ##calculating sp_ages
    sp_ages <- simSpAge(branching_time, lambda, mu.s[i], rho.s[j], reps)
    df <- data.frame(sp_ages = sp_ages,
                     mu = rep(mu.s[i], reps),
                     rho = rep(rho.s[j], reps))
    ##data general list
    data.general[[i]][[j]] = df
    
    prob_node_age <- ageFunction(lambda, mu.s[i], rho.s[j], branching_time, branching_time) 
    # the probability that t = branching_time
    prob_not_node_age <- integrate(ageDensityFunction,
                                   lambda = lambda, mu = mu.s[i],
                                   rho = rho.s[j], 
                                   lower = 0, upper = branching_time)$value
    t <- seq(0, branching_time, length.out = 1001)
    d <- sapply(t, ageDensityFunction, lambda = lambda, mu = mu.s[i],
                rho = rho.s[j])
    
    
    ##now extracting the probability of the branching time and obtaining a scalar
    scalar = 1 - mean(sp_ages == branching_time)
    
    ##the density is multiplied by scalar (the probability of not being the branching times)
    ##and by 0.2 due to defined histogram breaks
    
    den_t <- data.frame(t = t, 
                        den = d* 0.2 *scalar/  prob_not_node_age,
                        mu = rep(mu.s[i], length(t)),
                        rho = rep(rho.s[j], length(t)))
    
    ###den_general
    dent.list[[i]][[j]] <- den_t
    
  }  
  
}



for(i in seq_along(data.general)){
  data.general[[i]] <- do.call("rbind", data.general[[i]])
  dent.list[[i]] <- do.call("rbind", dent.list[[i]])
}

##general simulations datasets
data.general <- do.call("rbind", data.general)
dent.list <- do.call("rbind", dent.list)


mean(sp_ages == branching_time) # the frequency of simulations with age equal to v


##########Figures############

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

# Fully sampled -----------------------------------------------------------

######low extinction; fully sampled
data.low.full <- data.general %>% filter(mu ==0.05, 
                                         rho == 1)
data.low.full$extinction <- "Low extinction"

##probability of the branching time
prob.mean.low.full <- mean(data.low.full$sp_ages == branching_time)


###analytical probability
prob.an.low.full = ageFunction(lambda =0.3, mu = 0.05, rho = 1, t = 3,v = 3) # the analytical probability

den.low.full <- dent.list %>% filter(mu ==0.05, 
                                     rho == 1)

##ggplot

low.full <- data.low.full %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.low.full)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.low.full, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.low.full,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.low.full, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8, alpha = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.low.full,
                   yend = prob.an.low.full), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.015, 0.8), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.93)+
  mynamestheme



######Intermediate extinction; fully sampled
data.Intermediate.full <- data.general %>% filter(mu ==0.15, 
                                                  rho == 1)

data.Intermediate.full$extinction <- "Intermediate extinction"

##probability of the branching time
prob.mean.int.full <- mean(data.Intermediate.full$sp_ages == branching_time)

###analytical probability
prob.an.int.full = ageFunction(lambda = 0.3, mu = 0.15, rho = 1,
                                           t = 3,v = 3)

den.Intermediate.full <- dent.list %>% filter(mu ==0.15, 
                                              rho == 1)



##ggplot
Intermediate.full <- data.Intermediate.full %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.int.full)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.int.full, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.int.full,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.Intermediate.full, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.int.full,
                   yend = prob.an.int.full), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.030, 0.65), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.8)+
  mynamestheme



######high extinction; fully sampled
data.high.full <- data.general %>% filter(mu ==0.25, 
                                          rho == 1)
data.high.full$extinction <- "High extinction"
data.high.full$samp <- "Fully sampled"

##probability of the branching time
prob.mean.hig.full <- mean(data.high.full$sp_ages == branching_time)

###analytical probability
prob.an.high.full = ageFunction(lambda =0.3, mu = 0.25, rho = 1, t = 3,v = 3)

den.high.full <- dent.list %>% filter(mu ==0.25, 
                                      rho == 1)



##ggplot
high.full <- data.high.full %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.hig.full)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.hig.full, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.hig.full,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.high.full, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.high.full,
                   yend = prob.an.high.full), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.045, 0.55), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.66)+
  mynamestheme



# 25% missing species -----------------------------------------------------

######low extinction; 75% sampled
data.low.f25 <- data.general %>% filter(mu ==0.05, 
                                        rho == 0.75)
data.low.f25$extinction <- "Low extinction"

##probability of the branching time
prob.mean.low.f25 <- mean(data.low.f25$sp_ages == branching_time)

###analytical probability
prob.an.low.f25 = ageFunction(lambda =0.3, mu = 0.05, rho = 0.75, t = 3,v = 3)

den.low.f25 <- dent.list %>% filter(mu ==0.05, 
                                    rho == 0.75)



##ggplot
low.f25 <- data.low.f25 %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.low.f25)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.low.f25, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.low.f25,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.low.f25, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.low.f25,
                   yend = prob.an.low.f25), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.035, 0.60), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.70)+
  mynamestheme



######intermediate extinction; 75% sampled
data.intermediate.f25 <- data.general %>% filter(mu ==0.15, 
                                                 rho == 0.75)
data.intermediate.f25$extinction <- "Intermediate extinction"

##probability of the branching time
prob.mean.int.f25 <- mean(data.intermediate.f25$sp_ages == branching_time)

###analytical probability
prob.an.int.f25 = ageFunction(lambda =0.3, mu = 0.15, rho = 0.75, t = 3,v = 3)

den.intermediate.f25 <- dent.list %>% filter(mu ==0.15, 
                                             rho == 0.75)



##ggplot
intermediate.f25 <- data.intermediate.f25 %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.int.f25 )),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.int.f25 , 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.int.f25 ,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.intermediate.f25, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.int.f25,
                   yend = prob.an.int.f25), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.035, 0.53), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.60)+
  mynamestheme



######high extinction; 75% sampled
data.high.f25 <- data.general %>% filter(mu ==0.25, 
                                         rho == 0.75)
data.high.f25$extinction <- "High extinction"
data.high.f25$samp <- "25% missing species"

##probability of the branching time
prob.mean.high.f25 <- mean(data.high.f25$sp_ages == branching_time)

###analytical probability
prob.an.high.f25 = ageFunction(lambda =0.3, mu = 0.25, rho = 0.75, t = 3,v = 3)

den.high.f25 <- dent.list %>% filter(mu ==0.25, 
                                     rho == 0.75)



##ggplot
high.f25 <- data.high.f25 %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.high.f25)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.high.f25, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.high.f25,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.high.f25, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.high.f25,
                   yend = prob.an.high.f25), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.05, 0.45), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.51)+
  mynamestheme



# sampling 50% ------------------------------------------------------------

######low extinction; 75% sampled
data.low.f50 <- data.general %>% filter(mu ==0.05, 
                                        rho == 0.50)
data.low.f50$extinction <- "Low extinction"

##probability of the branching time
prob.mean.low.f50 <- mean(data.low.f50$sp_ages == branching_time)

###analytical probability
prob.an.low.f50 = ageFunction(lambda =0.3, mu = 0.05, rho = 0.50, t = 3,v = 3)

den.low.f50 <- dent.list %>% filter(mu ==0.05, 
                                    rho == 0.50)



##ggplot
low.f50 <- data.low.f50 %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.low.f50)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.low.f50, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.low.f50,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.low.f50, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.low.f50,
                   yend = prob.an.low.f50), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.07, 0.435), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.48)+
  mynamestheme



######intermediate extinction; 50% sampled
data.intermediate.f50 <- data.general %>% filter(mu ==0.15, 
                                                 rho == 0.50)
data.intermediate.f50$extinction <- "Intermediate extinction"

##probability of the branching time
prob.mean.int.f50 <- mean(data.intermediate.f50$sp_ages == branching_time)

###analytical probability
prob.an.int.f50 = ageFunction(lambda =0.3, mu = 0.15, rho = 0.50, t = 3,v = 3)

den.intermediate.f50 <- dent.list %>% filter(mu ==0.15, 
                                             rho == 0.50)



##ggplot
intermediate.f50 <- data.intermediate.f50 %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.int.f50)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.int.f50, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.int.f50,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.intermediate.f50, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.int.f50,
                   yend = prob.an.int.f50), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.07, 0.365), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.41)+
  mynamestheme


######high extinction; 50% sampled
data.high.f50 <- data.general %>% filter(mu ==0.25, 
                                         rho == 0.50)
data.high.f50$extinction <- "High extinction"
data.high.f50$samp <- "50% missing species"

##probability of the branching time
prob.mean.high.f50 <- mean(data.high.f50$sp_ages == branching_time)

###analytical probability
prob.an.high.f50 = ageFunction(lambda =0.3, mu = 0.25, rho = 0.50, t = 3,v = 3)

###density
den.high.f50 <- dent.list %>% filter(mu ==0.25, 
                                     rho == 0.50)



##ggplot
high.f50 <- data.high.f50 %>% filter(sp_ages < branching_time) %>% 
  ggplot(aes(sp_ages))+
  geom_histogram(aes(y=..count../sum(..count..)*(1 - prob.mean.high.f50)),
                 colour="black", fill = "white", 
                 breaks=seq(0, 3, by = 0.2)) +
  geom_segment(x = 3, xend = 3, y = 0, yend = prob.mean.high.f50, 
               colour = "darkgrey", linewidth = 0.8)+
  geom_point(x = 3, y =prob.mean.high.f50,
             colour = "darkgrey",
             size = 3)+
  geom_line(data = den.high.f50, aes(x = t, y = den),
            colour = "red",
            linewidth = 0.8)+
  geom_segment(aes(x = 2.9, xend = 3.1, y = prob.an.high.f50,
                   yend = prob.an.high.f50), colour = "red",
               linewidth = 0.8)+
  scale_y_break(c(0.07, 0.32), scales= "free")+
  facet_grid(~extinction)+
  theme_bw()+
  ylab(NULL)+
  xlab(NULL)+
  ylim(0,0.37)+
  mynamestheme



#############plotting 

##plot grid
pdf("text/supplementary/dist.ex.rho.2.pdf", width = 8.5)


##left
left <- text_grob("Relative Probability",
                  family = "serif", size = 15,  face = "bold",
                  rot = 90)

bottom <- text_grob("Ages",
                    family = "serif", size = 15,  face = "bold",
                    rot = 0)

##grid for plotting both figures
full <- grid.arrange(print(low.full), print(Intermediate.full),
                     print(high.full),
                     ncol = 3)

f25 <- grid.arrange(print(low.f25), print(intermediate.f25), print(high.f25),
                    ncol = 3)
f50 <- grid.arrange(print(low.f50), print(intermediate.f50),
                    print(high.f50),
                    ncol = 3)
grid.arrange(full, f25, f50, nrow = 3,
             left = left, bottom = bottom)

dev.off()
