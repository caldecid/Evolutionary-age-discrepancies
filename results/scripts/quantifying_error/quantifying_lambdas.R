# Quantifying the error. Figure 3 and 4 -----------------------------------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

###calling lambda speciation data set
q.error.lam <- read_csv(file = file.path(getwd(), 
                                         pro, q_error, "q.error.lambda.csv"))

q.error.scale <- q.error.lam %>% mutate(rTrue.age = True.age/root.age,
                                        rEstimated.age = Estimated.age/root.age)
##mean absolute error, ratio error and MAPE, by mode of speciation
q.mode.error <- q.error.lam %>% 
  group_by(speciation) %>%
  summarise(mean_abs = mean(abs(True.age - Estimated.age)),
            sd_m = sd(abs(True.age - Estimated.age)),
            ratio = mean(log10(Estimated.age/True.age)),
            ratio_sd = sd(log10(Estimated.age/True.age)),
            mape = mean(abs(True.age - Estimated.age)/True.age*100),
            sd_mape = sd(abs(True.age - Estimated.age)/True.age*100))


##overall mean absolute error and ratio error
q.overal.error <- q.error.lam %>% 
  summarise(mean_abs = mean(abs(True.age - Estimated.age)),
            sd_m = sd(abs(True.age - Estimated.age)),
            ratio = mean(log10(Estimated.age/True.age)),
            ratio_sd = sd(log10(Estimated.age/True.age)),
            mape = mean(abs(True.age - Estimated.age)/True.age*100),
            sd_mape = sd(abs(True.age - Estimated.age)/True.age*100)) %>% 
  mutate(speciation = "overall", .before = mean_abs)

##merging the two dataframes  
error.quantification <- rbind(q.overal.error, q.mode.error)          

##saving the error quantification
write_csv(error.quantification, file = file.path(getwd(), pro, q_error,
                                                 "error_quantification.csv"))


##Figure faceting True age vs Phylogenetic age by budding and bifurcating 

##defining fonts, faces and sizes for all the graphics

#themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (11)),
                      axis.title = element_text(family = "serif", size = (12),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (9),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (9)),
                      legend.position = c(0.1, 0.8),
                      legend.background = element_rect(fill="#ffffffaa",
                                                size=.3, linetype="dotted"))

# scaled ages -------------------------------------------------------------



# fig 3 scale true vs phylo without intermediate turnover -----------------
png("text/figures/Figure3.True.vs.Phylo.turnover.png", width = 17,
    height = 12, units = "cm", 
    pointsize = 8, res = 300)



scale.turnover <- q.error.scale %>% filter(speciation %in% 
                                             c("Bifurcating", "Budding")) %>% 
  slice(sample(1:n())) %>%
  mutate(turnover_bin = cut_interval(turnover,length = 0.25)) %>% 
  filter(turnover <= 0.25 | turnover > 0.75) %>% 
  ggplot(aes(rTrue.age, rEstimated.age, turnover_bin))+
  geom_point(alpha = 0.3, aes(color = turnover_bin))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  xlim(0,1)+
  ylim(0,1)+
  xlab("True age")+
  ylab("Phylogenetic age")+
  scale_color_manual(values = c("#d95f02", "#7570b3"),
                     name="Extinction fraction",
                     breaks=c("[0,0.25]", "(0.75,1]"),
                     labels=c("Low [0-0.25]",
                              "High (0.75-0.99]"))+
  facet_wrap(~speciation)+
  theme_bw()

scale.turnover + mynamestheme 

dev.off()

####anagenetic True vs phylogenetic

##SM Figure True vs Phylo anagenetic without intermediate turnover

png("text/supplementary/FigureSM1.scale.anagenetic.turnover.png", width = 17,
    height = 12, units = "cm", 
    pointsize = 8, res = 300)

scale.turnover.ana <- q.error.scale %>% filter(speciation %in% 
                         c("Anagenetic-budding", "Anagenetic-bifurcating")) %>% 
  slice(sample(1:n())) %>%
  mutate(turnover_bin = cut_interval(turnover,length = 0.25)) %>% 
  filter(turnover <= 0.25 | turnover > 0.75) %>% 
  ggplot(aes(rTrue.age, rEstimated.age, turnover_bin))+
  geom_point(alpha = 0.3, aes(color = turnover_bin))+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  xlim(0,1)+
  ylim(0,1)+
  xlab("True age")+
  ylab("Phylogenetic age")+
  scale_color_manual(values = c("#d95f02", "#7570b3"),
                     name="Extinction fraction",
                     breaks=c("[0,0.25]", "(0.75,1]"),
                     labels=c("Low [0-0.25]",
                              "High (0.75-0.99]"))+
  facet_wrap(~speciation)+
  theme_bw()

scale.turnover.ana + mynamestheme

dev.off()


# MAPE Figure 4 --------------------------------------------------------------

##MAPE grouped by mu, div, turnover and speciation modes
q.ext.bud.bif <- q.error.lam %>% filter(speciation %in% 
                              c("Bifurcating", "Budding")) %>% 
  slice(sample(1:n())) %>% 
  group_by(mu,speciation, div,
           turnover, lambda) %>% 
  summarise(mape = mean(abs(True.age - Estimated.age)/True.age)*100) 

q.ext.bud.bif$lambda <- as.factor(q.ext.bud.bif$lambda)


###model
log.model = lm(turnover ~ log(mape), q.ext.bud.bif)

###Figure MAPE vs turnover
png("text/figures/Figure4.MAPE.vs.turnover.png", width = 15, height = 10, units = "cm", 
    pointsize = 8, res = 300)

mape.turn <- ggplot(q.ext.bud.bif, aes(x = turnover, y = mape, color = lambda))+
  geom_point(alpha = 0.4, aes(color = lambda))+
  geom_smooth(aes(color = lambda), linewidth = 1.1, 
              linetype = 2, method = "lm", se = FALSE,
              formula = y ~ exp(x))+
  scale_color_manual("Speciation rates",
                     values = c("#66c2a5",
                                "#fc8d62",
                                "#8da0cb"))+
  facet_wrap(~speciation)+
 
  xlab("Extinction fraction")+
  ylab("MAPE")+
  theme_bw()

mape.turn + mynamestheme

dev.off()

##Supplementary-Anagenetic
q.ext.ana <- q.error.lam %>% filter(speciation %in% 
                         c("Anagenetic-budding", "Anagenetic-bifurcating")) %>%
                              group_by(mu,speciation, div,
                                              turnover, lambda) %>% 
  summarise(mape = mean(abs(True.age - Estimated.age)/True.age)*100)  
  

q.ext.ana$lambda <- as.factor(q.ext.ana$lambda)
q.ext.ana$speciation <- as.factor(q.ext.ana$speciation)
levels(q.ext.ana$speciation) <- c("ana_bif", "ana_bud")

###labeler for ploting
ana <- c("ana_bif" = "Anagenetic-Bifurcating",
         "ana_bud" = "Anagenetic-Budding")

##MAPE vs Turnover anagenetic
png("text/supplementary/SM2.mape.vs.turnover.png", width = 17, height = 10, units = "cm", 
    pointsize = 8, res = 300)

mape.turn.ana <- q.ext.ana %>% filter(mape < 500) %>% 
  ggplot(aes(x = turnover, y = mape, lambda))+
  geom_point(alpha = 0.4, aes(color = lambda))+
  geom_smooth(aes(color = lambda), linewidth = 1.1, 
              linetype = 2, method = "lm", se = FALSE,
              formula = y ~ exp(x))+
  scale_color_manual("Speciation rates",
                     values = c("#66c2a5",
                                "#fc8d62",
                                "#8da0cb"))+
  facet_wrap(~speciation, labeller = as_labeller(ana))+
  xlab("Extinction fraction")+
  ylab("MAPE (%)")+
  theme_bw()

mape.turn.ana + mynamestheme

dev.off()




