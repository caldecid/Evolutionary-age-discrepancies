
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
predict.glm <- predict.glm(glm_corrector, newdata = test_data_log,
                           type = "response")

##adding the prediction to the test data set
test_data_log$predict.glm <- predict.glm


###adding the BNN results

results.bnn <- read_csv("results_bnn.csv")

test_data_log$predict.bnn <- results.bnn$bnn_mean

##calculating the residuals between the prediction and the true age (log)
test_data_log$residuals.glm <- test_data_log$predict.glm - test_data_log$True.age_log

test_data_log$residuals.bnn <- test_data_log$predict.bnn - test_data_log$True.age_log

##MSE between the corrected age and the true age
mean((test_data_log$residuals.glm^2))
mean((test_data_log$residuals.bnn^2))
##raw MSE (between the estimated age and the true age)
mean((test_data_log$True.age_log - test_data_log$Estimated.age_log)^2)

############graphics #############
##grouping the test set by mode of speciation and trees, pivoting longer for ggplot
results.general <- test_data_log %>% group_by(mode, anag, div) %>%
              summarise(mse.glm = mean(residuals.glm^2),
              mse.bnn = mean(residuals.bnn^2),
              mse.raw = mean((True.age_log- Estimated.age_log)^2)) %>% 
              mutate(speciation = case_when(mode ==0 & anag == 0 ~ "Budding",
                              mode == 1 & anag == 1 ~ "Anagenetic-cladogenetic",
                              mode == 0 & anag == 1~ "Anagenetic-budding",
                              mode == 1 & anag == 0 ~ "Cladogenetic"))%>% 
                  pivot_longer(cols = starts_with("mse"), names_to = "model",
                                    names_prefix = "mse.", values_to = "MSE")

results.general$speciation <- as.factor(results.general$speciation)
results.general$model <- as.factor(results.general$model)

##ordering the factors for perfoming graphics
results.general$model <- factor(results.general$model,
                                levels = c("raw", "glm", "bnn"))

results.general$speciation <- factor(results.general$speciation,
                          levels = c("Budding","Anagenetic-budding",
                                     "Cladogenetic", "Anagenetic-cladogenetic"))
##graphics
ggplot(results.general, aes(x = speciation, y = MSE, color = model))+
  geom_boxplot()+
  ylim(0,1)



