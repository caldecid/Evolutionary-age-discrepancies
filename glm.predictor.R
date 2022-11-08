
# Models for correcting species age -------------------------------------------

##packages
library(tidyverse)
library(readr)
library(FossilSim)

##calling training data

train.predictors <- read_delim("ages.training.predictors.txt")

train.response <- read_csv("ages.training.response.txt")

train_data <- cbind(train.response, train.predictors)

##transforming variables and eliminating the transformed variables 
train_data_log <- train_data %>% transmute(True.age_log = log10(True.age),
                                      Estimated.age_log = log10(Estimated.age),
                                      root.age_log = log10(root.age),
                                      div = div,
                                      turnover = turnover,
                                      number.sp_log = log10(number.sp),
                                      mode = mode,
                                      N_sisters_log = log10(N_sisters),
                                      anag = anag) 
             
                                        

##fitting the Linear model

formula <- True.age_log~.


lm_corrector <- lm(formula, data = train_data_log)

summary(lm_corrector)

##calling test data
test.predictors <- read_csv("ages.test.predictors.txt")

test.response <- read_csv("ages.test.response.txt")

test_data <- cbind(test.response, test.predictors)

##transforming variables and eliminating the transformed variables
##also using the div and turnover rates estimated from the phylogenies as 
#the ones that enter in the model

test_data_log <- test_data %>%  transmute(True.age_log = log10(True.age),
                                       Estimated.age_log = log10(Estimated.age),
                                          root.age_log = log10(root.age),
                                          div = div.est,
                                          turnover = turnover.est,
                                          number.sp_log = log10(number.sp),
                                          mode = mode,
                                          N_sisters_log = log10(N_sisters),
                                          anag = anag) 




##predicting age with the lm 
pred.lm <- predict(lm_corrector, newdata = test_data_log,
                   interval = "prediction")

pred.lm <- as.data.frame(pred.lm)




##adding the lm prediction to the test data set
test_data_log$predict.lm <- pred.lm$fit


##adding the lm Confidence intervals
test_data_log$lm_low <- pred.lm$lwr

test_data_log$lm_up <- pred.lm$upr


