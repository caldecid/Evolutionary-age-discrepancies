
# GLM correcting estimated age --------------------------------------------

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
             
                                        

##fitting the glm model

formula <- True.age_log~.

glm_corrector <- glm(formula, data = train_data_log)

summary(glm_corrector)

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


##predicting the age
predict.age <- predict.glm(glm_corrector, newdata = test_data_log,
                           type = "response")

##adding the prediction to the test data set
test_data_log$predict <- predict.age

##calculating the residuals between the prediction and the true age (log)
test_data_log$residuals <- test_data_log$predict - test_data_log$True.age_log


##MSE between the corrected age and the true age
mean((test_data_log$residuals^2))

##raw MSE (between the estimated age and the true age)
mean((test_data_log$True.age_log - test_data_log$Estimated.age_log)^2)

############graphics #############
##grouping the test set by mode of speciation and trees, pivoting longer for ggplot
results.glm <- test_data_log %>% group_by(mode, anag, div) %>%
              summarise(mse.glm = mean(residuals^2),
              mse.raw = mean((True.age_log- Estimated.age_log)^2)) %>% 
              mutate(speciation = case_when(mode ==0 & anag == 0 ~ "Budding",
                              mode == 1 & anag == 1 ~ "Anagenetic-cladogenetic",
                              mode == 0 & anag == 1~ "Anagenetic-budding",
                              mode == 1 & anag == 0 ~ "Cladogenetic"))%>% 
                  pivot_longer(cols = starts_with("mse"), names_to = "model",
                                    names_prefix = "mse.", values_to = "MSE")

results.glm$speciation <- as.factor(results.glm$speciation)
results.glm$model <- as.factor(results.glm$model)

##ordering the factors for perfoming graphics
results.glm$model <- factor(results.glm$model, levels = c("raw", "glm"))
results.glm$speciation <- factor(results.glm$speciation,
                          levels = c("Budding","Anagenetic-budding",
                                     "Cladogenetic", "Anagenetic-cladogenetic"))
##graphics
ggplot(results.glm, aes(x = speciation, y = MSE, color = model))+
  geom_boxplot()
  
