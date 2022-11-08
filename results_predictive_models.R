
# Results predictive models -----------------------------------------------

##calling the lm results
results <- read_csv("results_lm.csv")

###adding the BNN results

results.bnn <- read_csv("results_bnn.csv")


##BNN prediction
results$bnn.mean <- results.bnn$bnn_mean

##BNN confidence intervals
#low
results$bnn.lowCI <- results.bnn$bnn_lowCI

#high 
results$bnn.upCI<- results.bnn$bnn_upCI

##calculating the residuals between the prediction and the true age (log)


results$residuals.lm <- results$lm.mean - results$true.age.log

results$residuals.bnn <- results$bnn.mean - results$true.age.log



##MSE between the corrected age and the true age
mean((results$residuals.bnn^2))
mean((results$residuals.lm^2))

##raw MSE (between the estimated age and the true age)
mean((results$true.age.log - results$Estimated.age.log)^2)

############graphics #############

##for the boxplot

##grouping the test set by mode of speciation and trees, pivoting longer for ggplot
results.general <- results %>% group_by(mode, anag, div) %>%
  summarise(mse.lm = mean(residuals.lm^2),
            mse.bnn = mean(residuals.bnn^2),
            mse.raw = mean((true.age.log- Estimated.age.log)^2)) %>% 
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
                                levels = c("raw", "lm", "bnn"))

results.general$speciation <- factor(results.general$speciation,
                                     levels = c("Budding","Anagenetic-budding",
                                                "Cladogenetic", "Anagenetic-cladogenetic"))
##graphics
#themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)))

##pdf object
png("boxplot.mse.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

box.mse <- ggplot(results.general, aes(x = speciation, y = MSE, color = model))+
           geom_boxplot()+
           xlab("Speciation mode")+
           ylab("Mean Squared Error")+
           #geom_jitter(size = 1, alpha=0.5)+
            ylim(0,1)+
            theme_bw()

box.mse + mynamestheme

dev.off()

###for the species age coverage graphic
results.species <- results %>% select(true.age.log, mode, anag,
                                      lm.mean, lm.lowCI, lm.upCI, bnn.mean, 
                                      bnn.lowCI, bnn.upCI)%>% 
  mutate(lm_coverage = if_else(true.age.log >= lm.lowCI &
                                 true.age.log <= lm.upCI, 1, 0),
         bnn_coverage = if_else(true.age.log >= bnn.lowCI &
                                  true.age.log <= bnn.upCI, 1, 0)) %>% 
  mutate(speciation = case_when(mode ==0 & anag == 0 ~ "Budding",
                                mode == 1 & anag == 1 ~ "Anagenetic-cladogenetic",
                                mode == 0 & anag == 1~ "Anagenetic-budding",
                                mode == 1 & anag == 0 ~ "Cladogenetic")) %>% 
  group_by(speciation) %>% 
  summarise(lm_cov_sum = sum(lm_coverage),
            bnn_cov_sum = sum(bnn_coverage),
            n = n()) %>% 
  mutate(per_cov_lm = lm_cov_sum/n * 100,
         per_cov_bnn = bnn_cov_sum/n * 100,
         per_cov_raw = lm_cov_sum - lm_cov_sum) %>% 
  pivot_longer(cols = starts_with("per"),
               names_to = "model",
               names_prefix = "per_cov_",
               values_to = "sp_coverage")

##ordering the factors for perfoming graphics
results.species$model <- factor(results.species$model,
                                levels = c("raw", "lm", "bnn"))

results.species$speciation <- factor(results.species$speciation,
                                     levels = c("Budding","Anagenetic-budding",
                                                "Cladogenetic", "Anagenetic-cladogenetic"))


###graphics
png("coverage.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

coverage <- ggplot(results.species, aes(x = speciation, y = sp_coverage, fill = model))+ 
           geom_bar(stat = "identity", position="dodge") +
            xlab("Speciation mode")+
           ylab("Species age coverage (%)")+
           theme_classic()

coverage + mynamestheme

dev.off()




