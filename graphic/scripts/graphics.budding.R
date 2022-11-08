# Graphics for budding simulations ------------------------------------------------

##defining fonts, faces and sizes for all the graphics
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)))


# turnover rates vs ratio -----------------------------------------------
png("turnover.vs.ratio.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)


##ggplot basic information

tvr <- ggplot(data = budding.df, aes(x = turnover, y= ratio_error))+
  
  geom_point(alpha = 0.1)+
  #geom_smooth(method = "loess", se = F)+
  #geom_abline(slope = 1, intercept = 0, size= 1, color = "red")+
  #geom_abline(slope = 1, intercept = 0)+
  ggtitle("Budding speciation")+
  xlab("Turnover rates (Extinction/Speciation)")+
  #xlim(0, 20)+
  #ylim(0, 30)+
  ylab("log10(Estimated age/ True age)")+
  theme_minimal()




tvr + mynamestheme

dev.off()


# true age vs estimated age -----------------------------------------------


png("true_age.vs.estimated_age.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)


##ggplot basic information

true.vs.phylo <- ggplot(data = budding.df, aes(x = True.age,
                                                      y= Estimate.age))+
  
  geom_point(alpha = 0.1)+
  #geom_smooth(method = "loess", se = F)+
  #geom_abline(slope = 1, intercept = 0, size= 1, color = "red")+
  geom_abline(slope = 1, intercept = 0)+
  ggtitle("Budding speciation")+
  xlab("True age")+
  xlim(0, 20)+
  ylim(0, 30)+
  ylab("Estimated age")




true.vs.phylo + mynamestheme

dev.off()


# true age vs error -------------------------------------------------------
png("true_age.vs.error.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

true.vs.error <- ggplot(data = budding.df, aes(x = True.age,
                                                      y= relative_error))+
  geom_point(alpha = 0.1)+
  geom_abline(slope = 0, intercept = 0)+
  ggtitle("Budding speciation")+
  ylim(-1, 20)+
  xlab("True age")+
  ylab("relative error")+
  theme_minimal()

true.vs.error + mynamestheme

dev.off()

# Estimated age vs error -------------------------------------------------------
png("Estimated_age.vs.error.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

est.vs.error <- ggplot(data = budding.df, aes(x = Estimate.age,
                                                     y= relative_error))+
  geom_point(alpha = 0.1)+
  #geom_abline(slope = 0, intercept = 0)+
  ggtitle("Budding speciation")+
  ylim(-1, 7.5)+
  xlim(0,15)+
  xlab("Estimated age")+
  ylab("relative error")+
  theme_minimal()

est.vs.error + mynamestheme

dev.off()

# Estimated age vs ratio -------------------------------------------------------
png("Estimated_age.vs.ratio.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

est.vs.error <- ggplot(data = budding.df, aes(x = Estimate.age,
                                              y= ratio_error))+
  geom_point(alpha = 0.1)+
  #geom_abline(slope = 0, intercept = 0)+
  ggtitle("Budding speciation")+
  ylim(-3, 2)+
  xlim(0,20)+
  xlab("Estimated age")+
  ylab("log10(Estimated age/True age)")+
  theme_minimal()

est.vs.error + mynamestheme

dev.off()

# True age vs ratio -------------------------------------------------------
png("True_age.vs.ratio.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

true.vs.ratio <- ggplot(data = budding.df, aes(x = True.age,
                                              y= ratio_error))+
  geom_point(alpha = 0.1)+
  #geom_abline(slope = 0, intercept = 0)+
  ggtitle("Budding speciation")+
  ylim(-3, 2)+
  xlim(0,15)+
  xlab("True age")+
  ylab("log10(Estimated age/True age)")+
  theme_minimal()

true.vs.ratio + mynamestheme

dev.off()

# true age vs error faceted by extinction rates ---------------------------

png("true_age.vs.error.extinction.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.extinction <- ggplot(data = budding.df,
                           aes(x = True.age, y= relative_error))+
  geom_point(alpha = 0.1)+
  ylim(-1, 10)+
  xlim(0, 20)+
  facet_wrap(~mu.int)+
  ggtitle("Budding speciation; Extinction rates intervals")+
  xlab("True age")+
  ylab("Relative error")+
  theme_minimal()

facet.extinction + mynamestheme

dev.off()

# true age vs ratio faceted by extinction rates ---------------------------

png("true_age.vs.ratio.extinction.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.extinction <- ggplot(data = budding.df,
                           aes(x = True.age, y= ratio_error))+
  geom_point(alpha = 0.1)+
  ylim(-3, 2)+
  xlim(0,20)+
  facet_wrap(~mu.int)+
  ggtitle("Budding speciation; Extinction rates intervals")+
  xlab("True age")+
  ylab("log10(Estimated age/ True age)")+
  theme_minimal()

facet.extinction + mynamestheme

dev.off()
# estimated age vs ratio faceted by extinction rates ---------------------------

png("estimated_age.vs.ratio.extinction.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.extinction <- ggplot(data = budding.df,
                           aes(x = Estimate.age, y= ratio_error))+
  geom_point(alpha = 0.1)+
  xlim(0, 20)+
  ylim(-3, 2)+
  facet_wrap(~mu.int)+
  ggtitle("Budding speciation; Extinction rates intervals")+
  xlab("Estimated age")+
  ylab("log10(Estimated age/ True age)")+
  theme_minimal()

facet.extinction + mynamestheme

dev.off()

# extinction vs absolute error -----------------------

png("extinction.vs.abs_error.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

facet.div <- ggplot(data = budding.df,
                    aes(x = mu, y= absolute_error))+
  geom_point(alpha = 0.1)+
  ylim(0, 25)+
  ggtitle("Budding speciation")+
  xlab("Extinction rates")+
  ylab("Absolute error")+
  theme_minimal()

facet.div + mynamestheme

dev.off()


# mean absolute error vs extinction rate

png("mean.abs.error.vs.extinction.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.ext <- ggplot(budding.mean.data, aes(x = mu, y = mean_abs_error))+
  geom_point()+
  geom_smooth(method = "loess", se = F)+
  #ggtitle("Budding speciation")+
  xlab("extinction rates")+
  ylab("mean absolute error")+
  theme_minimal()

mean.abs.ext + mynamestheme

dev.off()




########absolute error along different extinctions rates
png("absolute.error.hist.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

ta<- ggplot(data = budding.df, aes(absolute_error))+
  geom_histogram(bins = 35)+
  xlim(0,10)+
  ylim(0,100)+
  facet_wrap(~mu.int)+
  xlab("absolute error")+
  ggtitle("Budding speciation; Faceted by extinction rates")+
  theme_minimal()

ta + mynamestheme

dev.off()

#########estimated age along different extinction rates

png("Estimate.age.hist.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

ea<- ggplot(data = budding.df, aes(Estimate.age))+
  geom_histogram(bins = 45)+
  xlim(0,10)+
  facet_wrap(~mu.int)+
  ggtitle("Budding speciation; Faceted by extinction rates")+
  xlab("Estimated age")+
  theme_minimal()

ea + mynamestheme

dev.off()

#########True age along different extinction rates

png("True.age.hist.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

ea<- ggplot(data = budding.df, aes(True.age))+
  geom_histogram(bins = 45)+
  xlim(0,10)+
  facet_wrap(~mu.int)+
  ggtitle("Budding speciation; Faceted by extinction rates")+
  xlab("True age")+
  theme_minimal()

ea + mynamestheme

dev.off()

# mean ratio vs extinction rates ----------------------------

png("mean.ratio.vs.ext.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.div <- ggplot(budding.mean.data, aes(x = mu, y = mean_ratio))+
  geom_point()+
  ggtitle("Budding speciation")+
  xlab("Extinction rates")+
  ylab("mean log10(estimated age/true age)")+
  theme_minimal()

mean.abs.div + mynamestheme

dev.off()






# mean absolute error vs turnover rates -----------------------------------
png("mean.abs.error.vs.turnover.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.tur <- ggplot(budding.mean.data, aes(x = turnover,
                                            y = mean_abs_error))+
  geom_point()+
  xlab("Turnover rates")+
  ylab("Mean absolute error")+
  theme_minimal()

mean.abs.tur + mynamestheme

dev.off()
