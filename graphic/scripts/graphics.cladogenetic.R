# Graphics for cladogenetic simulations ------------------------------------------------

##defining fonts, faces and sizes for all the graphics
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)))



# turnover rates vs ratio error -------------------------------------------

png("turnover.vs.ratio.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

tvr <- ggplot(data = cladogenetic.data, aes(x = turnover, y= ratio_error))+
  
  geom_point(alpha = 0.1)+
  #geom_smooth(method = "loess", se = F)+
  #geom_abline(slope = 1, intercept = 0, size= 1, color = "red")+
  #geom_abline(slope = 1, intercept = 0)+
  ggtitle("Cladogenetic speciation")+
  xlab("Turnover rates (Extinction/Speciation)")+
  #xlim(0, 20)+
  #ylim(0, 30)+
  ylab("log10(Estimated age/ True age)")+
  theme_minimal()




tvr + mynamestheme

dev.off()

# True age vs estimated age -----------------------------------------------

png("true_age.vs.estimated_age.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)


##ggplot basic information

true.vs.phylo <- ggplot(data = cladogenetic.data, aes(x = True.age,
                                                       y= Estimate.age))+
  
  geom_point(alpha = 0.1)+
  #geom_smooth(method = "loess", se = F)+
  #geom_abline(slope = 1, intercept = 0, size= 1, color = "red")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  ggtitle("Cladogenetic speciation")+
  xlab("True age")+
  #xlim(0, 50)+
  #ylim(0, 50)+
  ylab("Estimated age")




true.vs.phylo + mynamestheme

dev.off()


# true age vs error -------------------------------------------------------
png("true_age.vs.error.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

true.vs.error <- ggplot(data = cladogenetic.data, aes(x = True.age,
                                                       y= relative_error))+
  geom_point(alpha = 0.1)+
  geom_abline(slope = 0, intercept = 0)+
  ggtitle("Cladogenetic speciation")+
  ylim(0, 25)+
  xlab("True age")+
  ylab("relative error")+
  theme_minimal()

true.vs.error + mynamestheme

dev.off()

# Estimated age vs error -------------------------------------------------------
png("Estimated_age.vs.error.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

est.vs.error <- ggplot(data = cladogenetic.data, aes(x = Estimate.age,
                                                      y= relative_error))+
  geom_point(alpha = 0.1)+
  #geom_abline(slope = 0, intercept = 0)+
  ggtitle("Cladogenetic speciation")+
  ylim(0, 25)+
  xlab("Estimated age")+
  ylab("relative error")+
  theme_minimal()

est.vs.error + mynamestheme

dev.off()
# true age vs error faceted by extinction rates ---------------------------

png("true_age.vs.error.extinction.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.extinction <- ggplot(data = cladogenetic.data,
                           aes(x = True.age, y= relative_error))+
  geom_point(alpha = 0.1)+
  ylim(-1, 10)+
  facet_wrap(~mu.int)+
  ggtitle("Extinction rates intervals")+
  xlab("True age")+
  ylab("Relative error")+
  theme_minimal()

facet.extinction + mynamestheme

dev.off()

# true age vs ratio faceted by extinction rates ---------------------------

png("true_age.vs.ratio.error.extinction.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.extinction <- ggplot(data = cladogenetic.data,
                           aes(x = True.age, y= ratio_error))+
  geom_point(alpha = 0.1)+
  #ylim(0, 2.5)+
  xlim(0, 7.5)+
  facet_wrap(~mu.int)+
  ggtitle("Cladogenetic speciation; Extinction rates intervals")+
  xlab("True age")+
  ylab("log10(Estimated age/ True age)")+
  theme_minimal()

facet.extinction + mynamestheme

dev.off()
# estimated age vs ratio faceted by extinction rates ---------------------------

png("estimated_age.vs.ratio.extinction.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.extinction <- ggplot(data = cladogenetic.data,
                           aes(x = Estimate.age, y= ratio_error))+
  geom_point(alpha = 0.1)+
  #ylim(0, 2.5)+
  xlim(0, 15)+
  facet_wrap(~mu.int)+
  ggtitle("Cladogenetic speciation; Extinction rates intervals")+
  xlab("Estimated age")+
  ylab("log10(Estimated age/ True age)")+
  theme_minimal()

facet.extinction + mynamestheme

dev.off()

# Estimated age vs eror faceted by turnover rates -----------------------

png("turnover.vs.abs_error.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.div <- ggplot(data = cladogenetic.data,
                    aes(x = turnover, y= absolute_error))+
  geom_point(alpha = 0.1)+
  #ylim(0, 2.5)+
  #facet_wrap(~div.int)+
  ggtitle("Cladogenetic speciation")+
  xlab("Turnover rates")+
  ylab("Absolute error")+
  theme_minimal()

facet.div + mynamestheme

dev.off()


# mean absolute error vs extinction rate

png("mean.abs.error.vs.extinction.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.ext <- ggplot(clado.mean.data, aes(x = mu, y = mean_abs_error))+
  geom_point()+
  geom_smooth(method = "loess", se = F)+
  #ggtitle("Cladogenetic speciation")+
  xlab("extinction rates")+
  ylab("mean absolute error")+
  theme_minimal()

mean.abs.ext + mynamestheme

dev.off()

#absolute error vs extinction rate
png("abs.error.vs.extinction.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

abs.ext <- ggplot(cladogenetic.data, aes(x = mu, y = absolute_error))+
  geom_point(alpha = 0.2)+
  ggtitle("Cladogenetic speciation")+
  xlab("extinction rates")+
  ylab("absolute error")+
  theme_minimal()

abs.ext + mynamestheme

dev.off()



########absolute error along different extinctions rates
png("absolute.error.hist.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

ta<- ggplot(data = cladogenetic.data, aes(absolute_error))+
  geom_histogram(bins = 20)+
  xlim(0,10)+
  ylim(0,100)+
  facet_wrap(~mu.int)+
  ggtitle("Faceted by extinction rates")+
  theme_minimal()

ta + mynamestheme

dev.off()

#########phylo age along different extinction rates

png("Estimate.age.hist.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

ea<- ggplot(data = cladogenetic.data, aes(Estimate.age))+
  geom_histogram(bins = 45)+
  facet_wrap(~mu.int)+
  ggtitle("Faceted by extinction rates")+
  xlab("Estimated age")+
  theme_minimal()

ea + mynamestheme

dev.off()


# histogram relative errors faceted by extinction rates -------------------

png("hist.relative.error.facet.extinction.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

re<- ggplot(data = cladogenetic.data, aes(relative_error))+
  geom_histogram()+
  facet_wrap(~mu.int)+
  xlim(0,10)+
  ylim(0, 125)+
  #geom_vline(xintercept = 0)+
  ggtitle("Faceted by extinction rates")+
  xlab("Relative Error")+
  theme_minimal()

re + mynamestheme

dev.off()



# histogram relative errors faceted by diversification rates --------------

png("hist.relative.error.facet.diversification.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

rd<- ggplot(data = disc, aes(relative_error_esc))+
  geom_histogram()+
  facet_wrap(~div.int, scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("Faceted by diversification rates")+
  xlab("Relative Error")+
  theme_minimal()

rd + mynamestheme

dev.off()


# mean ratio vs extinction rates ----------------------------

png("mean.ratio.vs.ext.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.div <- ggplot(clado.mean.data, aes(x = mu, y = mean_ratio))+
  geom_point()+
  ggtitle("Cladogenetic speciation")+
  xlab("Extinction rates")+
  ylab("mean log10(estimated age/true age)")+
  theme_minimal()

mean.abs.div + mynamestheme

dev.off()


# median relative error vs extinction rates ---------------------------------

png("median.rel.error.vs.ext.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

median.rel.ext <- ggplot(clado.mean.data, aes(x = mu, y = median_relative_sc_error))+
  geom_point()+
  xlab("Extinction rates")+
  ylab("median relative error")+
  ylim(-1,0)+
  theme_minimal()

median.rel.ext + mynamestheme

dev.off()


# median relative error vs diversification rates --------------------------

png("median.rel.error.vs.div.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

median.rel.div <- ggplot(clado.mean.data, aes(x = div, y = median_relative_sc_error))+
  geom_point()+
  xlab("Diversification rates")+
  ylab("median relative error")+
  ylim(-1,0)+
  theme_minimal()

median.rel.div + mynamestheme

dev.off()


# histogram relative error faceted by estimated age ----------------------

png("hist.relative.error.facet.Estimate.age.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

ra<- ggplot(data = cladogenetic.data %>% filter(relative_error_esc < 5)
            , aes(relative_error_esc))+
  geom_histogram()+
  facet_wrap(~phylo.int, scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("Faceted by estimated age")+
  xlab("Relative Error")+
  theme_minimal()

ra + mynamestheme

dev.off()



# histogram relative error faceted by estimated age without extrem --------

##filtering extreme values
disc.2 <- cladogenetic.data %>% filter(Estimate.age < 0.45)
disc.2$phylo.int <- cut_interval(disc.2$Estimate.age, 4)

png("hist.relative.error.facet.age.without.extreme.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

rext<- ggplot(data = disc.2 %>% filter(relative_error_esc < 5)
              , aes(relative_error_esc))+
  geom_histogram()+
  facet_wrap(~phylo.int, scales = "free")+
  geom_vline(xintercept = 0)+
  ggtitle("Faceted by estimated age, without extreme values")+
  xlab("Relative Error")+
  theme_minimal()

rext + mynamestheme

dev.off()


# median relative error vs turnover rate ------------------------------------

png("median.rel.error.vs.turnover.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

median.rel.tur <- ggplot(clado.mean.data, aes(x = turnover,
                                                 y = median_relative_sc_error))+
  geom_point()+
  xlab("Turnover rates")+
  ylab("Median relative error")+
  ylim(-1,0)+
  theme_minimal()

median.rel.tur + mynamestheme

dev.off()


# mean absolute error vs turnover rates -----------------------------------
png("mean.abs.error.vs.turnover.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.tur <- ggplot(clado.mean.data, aes(x = turnover,
                                               y = mean_abs_esc_error))+
  geom_point()+
  xlab("Turnover rates")+
  ylab("Mean absolute error")+
  theme_minimal()

mean.abs.tur + mynamestheme

dev.off()

