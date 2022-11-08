
# graphics for anagenetic-budding speciation ------------------------------

ana.bud <- read.xlsx("ana.bud.training.dataset.xlsx", sheetIndex = 1)

##defining fonts, faces and sizes for all the graphics
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (12)),
                      axis.title = element_text(family = "serif", size = (13),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (12)))


###true vs estimated age
png("true.vs.estimated.ae.png", width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)


##ggplot basic information

true.vs.phylo <- ggplot(data = ana.bud, aes(x = True.age,
                                                      y= Estimated.age))+
  
  geom_point(alpha = 0.1)+
  #geom_smooth(method = "loess", se = F)+
  #geom_abline(slope = 1, intercept = 0, size= 1, color = "red")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  ggtitle("Anagenetic-budding speciation")+
  xlab("True age")+
  xlim(0, 10)+
  ylim(0, 10)+
  ylab("Estimated age")




true.vs.phylo + mynamestheme

dev.off()


# mean absolute error vs extinction rate
#group_by and then summary
mean.ana.bud <- ana.bud %>% group_by(div) %>% summarise(mean_abs_error = 
                                           abs(mean(True.age - Estimated.age)),
                                           mu = mean(-turnover*div/(turnover-1)))


png("mean.abs.error.vs.extinction.png", 
    width = 16, height = 12, units = "cm", 
    pointsize = 8, res = 300)

mean.abs.ext <- ggplot(mean.ana.bud, aes(x = mu, y = mean_abs_error))+
  geom_point()+
  #geom_smooth(method = "loess", se = F)+
  ggtitle("Anagenetic-budding speciation")+
  xlim(0, 0.4)+
  ylim(0,2)+
  xlab("extinction rates")+
  ylab("mean absolute error")+
  theme_minimal()

mean.abs.ext + mynamestheme

dev.off()

