##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))



# Testing the geometric function -----------------------------------------

###using the same data generated in the simulation of the conservation status

ages.total <- read_csv(file = file.path(getwd(), pro, c_status,
                                        "ages.total.csv"))

load(file = file.path(getwd(), pro, c_status, "trees.conservation.RData"))



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

##determining how many species are 100% accurate regarding the geometric fun
accurate.st <- ages.total %>%
  mutate(abs.accuracy.st = abs(rmode - rTrue.age)) %>% 
  filter(abs.accuracy.st < 1e-8) %>% group_by(extinction) %>%
  count()



png("text/figures/Figure6.accuracy.mode.png",
    width = 10, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f6 <- ages.total %>%
  mutate(abs.accuracy = abs(rmode - rTrue.age)) %>% 
  filter(abs.accuracy > 1e-8) %>% 
  ggplot(aes(x = rmode - rTrue.age))+
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
  xlab("Most probability age - True age")+
  ylab("Species count")+
  ggtitle("Geometric function accuracy")

f6 + mynamestheme

dev.off()

##accuracy of the mean prob age regarding true age

##determining how many species are 100% accurate regarding the geometric fun


accurate.mean <- ages.total %>%
  mutate(abs.accuracy.mean = abs(rmean - rTrue.age)) %>% 
  filter(abs.accuracy.mean < 1e-8) %>% 
  group_by(extinction) %>%
  count() ##there are no correct estimations



png("text/figures/Figure6.1.accuracy.mean.png",
    width = 10, height = 15, units = "cm", 
    pointsize = 8, res = 300)

f6.1 <- ages.total %>%
  mutate(abs.accuracy.mean = abs(rmean - rTrue.age)) %>% 
  ggplot(aes(x = rmean - rTrue.age))+
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
  ggtitle("Geometric function accuracy")

f6.1 + mynamestheme

dev.off()

######################estimating coverage using the Confidence intervals

coverage.st <- ages.total%>%
                  filter(abs(True.age - modal) < 1e-08 | 
                           True.age >= lwr_ci & True.age <= upr_ci) %>% 
                 count(extinction) %>% 
                 mutate(cov = n/10000 * 100) 

png("text/figures/Figure6.2.coverage.CI.geometric.png",
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






