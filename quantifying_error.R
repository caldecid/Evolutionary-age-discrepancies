
# Quantifying the error ---------------------------------------------------

###using the training data

##############mean absolute error######################################

train.error <- train_data %>%
                 mutate(speciation = case_when(mode ==0 & anag == 0 ~ "Budding",
                              mode == 1 & anag == 1 ~ "Anagenetic-cladogenetic",
                                    mode == 0 & anag == 1~ "Anagenetic-budding",
                                        mode == 1 & anag == 0 ~ "Cladogenetic"),
                        True.log = log10(True.age),
                        Est.log = log10(Estimated.age))


##general error
train.mode.error <- train.error %>% 
  group_by(speciation) %>% summarise(mean_abs = mean(abs(True.age - Estimated.age)),
                                     sd_m = sd(abs(True.age - Estimated.age)),
                                     mean_log_abs = mean(abs(True.log - Est.log)),
                                     sd_log = sd(abs(True.log - Est.log)),
                                     mse = mean((True.age - Estimated.age)^2),
                                     sd_mse = sd((True.age - Estimated.age)^2),
                                     mse_log = mean((True.log - Est.log)^2),
                                     sd_mse_log = mean((True.log - Est.log)^2),
                                     ratio = mean(log10(Estimated.age/True.age)),
                                     ratio_sd = sd(log10(Estimated.age/True.age)))
  
  
  
          



##########error by mode of speciation

train.mode.error <- train.error %>% group_by(div, speciation) %>% 
           summarise(mean_absolute_error = mean(abs(True.age - Estimated.age)),
            sd_absolute_error = sd(abs(True.age - Estimated.age))) %>% 
             group_by(speciation) %>% 
            summarise(mean_abs = mean(mean_absolute_error),
                      sd_abs = mean(sd_absolute_error))
