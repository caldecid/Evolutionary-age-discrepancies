
# Quantifying the error ---------------------------------------------------

##libraries
library(tidyverse)
library(FossilSim)

###calling training data

train.predictors <- read_csv("ages.training.predictors.txt")

train.response <- read_csv("ages.training.response.txt")

train_data <- cbind(train.response, train.predictors)

##############mean absolute error######################################

train.error <- train_data %>%
                 mutate(speciation = case_when(mode ==0 & anag == 0 ~ "Budding",
                              mode == 1 & anag == 1 ~ "Anagenetic-cladogenetic",
                                    mode == 0 & anag == 1~ "Anagenetic-budding",
                                        mode == 1 & anag == 0 ~ "Cladogenetic"))

train.error$speciation <- factor(train.error$speciation,
                                 levels = c("Budding","Anagenetic-budding",
                                  "Cladogenetic", "Anagenetic-cladogenetic"))

##mean absolute error and ratio error by mode of speciation
train.mode.error <- train.error %>% 
                                 group_by(speciation) %>%
                      summarise(mean_abs = mean(abs(True.age - Estimated.age)),
                                sd_m = sd(abs(True.age - Estimated.age)),
                                 ratio = mean(log10(Estimated.age/True.age)),
                                 ratio_sd = sd(log10(Estimated.age/True.age)))


 ##overall mean absolute error and ratio error
train.overal.error <- train.error %>% 
                     summarise(mean_abs = mean(abs(True.age - Estimated.age)),
                      sd_m = sd(abs(True.age - Estimated.age)),
                      ratio = mean(log10(Estimated.age/True.age)),
                      ratio_sd = sd(log10(Estimated.age/True.age))) %>% 
                      mutate(speciation = "overall", .before = mean_abs)

##merging the two dataframes  
error.quantification <- rbind(train.overal.error, train.mode.error)          

##saving the error quantification
write_csv(error.quantification, file = "error_quantification.csv")

##graphics faceting True age vs Estimated age by speciation mode

##defining fonts, faces and sizes for all the graphics

#themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)))

##pdf object
png("true_age.vs.estimated_age.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)
##main
true.vs.est <-ggplot(data = train.error, aes(x = True.age, y = Estimated.age))+
              geom_point(alpha = 0.1)+
              geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
              xlim(0,50)+
              ylim(0,50)+
              xlab("True age (myrs)")+
              ylab("Estimated age (myrs)")+
              facet_wrap(~speciation)+
              theme_bw()

true.vs.est + mynamestheme

dev.off()



