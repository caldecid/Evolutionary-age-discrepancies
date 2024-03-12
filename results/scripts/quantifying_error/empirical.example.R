
# Fossil vs phylogenetic vs corrected age -----------------------------------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


###calling individual species
ind_sp <- read_delim(file.path(getwd(), pro, q_error, "ind_species.csv"))


##establishing different speciation and extinction scenarios from the diversification rates
###changing the value of the Homo sapiens phylogeny (Rivas-González et al. 2023)


ind_sp <- ind_sp %>% distinct(species, .keep_all = TRUE) %>%  
                  mutate(Estimated.age = replace(Estimated.age,
                                        Estimated.age == 6.175880, 10)) 
##just fossil data
ind_sp_fossil <- ind_sp %>% select(species, fossil.age) %>%
                                 filter(species %in% c("Acinonyx_jubatus",
                                                       "Ursus_arctos",
                                                       "Homo_sapiens",
                                                       "Balaena_mysticetus"))

##Calling the phylogenies

##A random mammal phylogeny from Upham et al. (2019)

mammals_phy <- read.tree("C:/Users/carlo/OneDrive/Desktop/Phylogenies/Mammals/mammals_phy.tre")

mammals_phy <- force.ultrametric(mammals_phy)

##mammals species names
mammals_tips <- tip.label(mammals_phy)

##root age
root_age <- tree.max(mammals_phy)


##calculate species ages for all mammals
mammals_age <- calculate_tip_ages(mammals_phy) %>% rename(species = tip,
                                                   Estimated.age = tip.age)

##matching fossil data and species age (changing the Homo Sapiens according toRivas-González et al. 2023)
ind_sp_ages <- left_join(ind_sp_fossil, mammals_age, by = "species") %>% 
                                mutate(Estimated.age = replace(Estimated.age,
                                       Estimated.age == 0.933, 10)) 



# correcting species ages -------------------------------------------------

###estimating extinction and speciation rate

#high extinction fraction
mammals_rates_high <- estimate_bd(mammals_phy,
                                 epsilon = 0.95, rho = 0.95)

##creating a tibble for posterior merging
mammals_high <- tibble(lambda = mammals_rates_high[1],
                      mu = mammals_rates_high[2],
                      epsilon = "high")

##intermediate epsilon
mammals_rates_int <- estimate_bd(mammals_phy,
                                epsilon = 0.50, rho = 0.95)

##creating a tibble for posterior merging
mammals_int <- tibble(lambda = mammals_rates_int[1],
                     mu = mammals_rates_int[2],
                     epsilon = "int")


##low epsilon
mammals_rates_low <- estimate_bd(mammals_phy,
                                epsilon = 0.15, rho = 0.95)

##creating a tibble for posterior merging
mammals_low <- tibble(lambda = mammals_rates_low[1],
                     mu = mammals_rates_low[2],
                     epsilon = "low")

##repeating each row the number of rows in mammal_age
mammals_div <- rbind(mammals_high[rep(seq_len(nrow(mammals_high)), nrow(mammals_age)), ],
                     mammals_int[rep(seq_len(nrow(mammals_int)), nrow(mammals_age)), ], 
                     mammals_low[rep(seq_len(nrow(mammals_low)), nrow(mammals_age)), ])

##repeating three times the mammal_ages df (3 extinction fractions)
mammals_age_rep <- mammals_age[rep(seq_len(nrow(mammals_age)), 3), ]

##binding
mammals_df = cbind(mammals_age_rep, mammals_div)

##changing the Homo sapiens estimated age according to Rivas-González et al.(2023)
homo_df = mammals_df %>% filter(species == "Homo_sapiens")

homo_df$Estimated.age = c(10, 10, 10)

mammals_df = mammals_df %>% filter(species != "Homo_sapiens")

mammals_df = rbind(mammals_df, homo_df)
 

# Corrected mean age ------------------------------------------------------

mean.age <- vector(mode = "numeric",
                           length = nrow(mammals_df))

for(i in 1:nrow(mammals_df)){
  mean.age[[i]] <- meanAgeFunction(lam = mammals_df$lambda[i],
                                           mu = mammals_df$mu[i],
                                           rho = 0.95,
                                           v = mammals_df$Estimated.age[i])
}

##binding datasets
mammals_emp <- cbind(mammals_df, mean.age)

##root age
mammals_emp$root = root_age

##standarizing ages
mammals_emp = mammals_emp %>% mutate(rPhylo = Estimated.age/root_age,
                                     rMean = mean.age/root_age)

mammals_emp$epsilon = factor(mammals_emp$epsilon, levels = c("low",
                                                             "int",
                                                                "high"),
                             ordered = TRUE)

##pivoting
mammals_emp_pivot <- mammals_emp %>% 
                          pivot_longer(cols = ends_with(".age"),
                                       names_to = "type",
                                       values_to = "ages")
##factor




##########Figure
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (11),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      axis.text.y = element_text(family = "serif", size = (10),
                                                 face = "italic"),
                      legend.title = element_text(family = "serif", size = (11),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (10)),
                      legend.background = element_rect(fill= "grey",
                                                       size=.5, linetype="dotted"),
                      legend.position = "bottom")




# low extinction plots ----------------------------------------------------

mammals_low = mammals_emp %>% filter(epsilon == "low")

##df for the selected species
ind_low = left_join(ind_sp_fossil, mammals_low, by = "species") %>%  
                  mutate(Estimated.age = replace(Estimated.age,
                                 Estimated.age == 6.175880, 10))



##all mammals phylo age vs corrected age
low.phy = ggplot(mammals_low, aes(x = Estimated.age, y = mean.age))+
           geom_point(alpha = 0.1)+
           xlim(0,11)+
           ylim(0, 11)+
           geom_point(data = ind_low, aes(x = Estimated.age, y = mean.age), 
                      size = 3.5, color = "red", alpha = 0.9)+
          xlab("Phylogenetic age (myr)")+
          ylab("Corrected age (myr)")+
          #ggtitle("Low extinction fraction")+
          theme_bw()+
          mynamestheme

##comparison phylogenetic, fossil, and corrected age

##pivoting for plotting
ind_low_pivot = ind_low %>% 
                    pivot_longer(cols = ends_with(".age"), names_to = "type",
                                 values_to = "ages")

##organizing factors
ind_low_pivot$type <- factor(ind_low_pivot$type,
                                 levels = c("mean.age",
                                            "fossil.age",
                                            "Estimated.age"),
                                 ordered = TRUE)

ind_low_pivot$species = factor(ind_low_pivot$species, 
                               levels = c("Ursus_arctos",
                                          "Acinonyx_jubatus",
                                          "Homo_sapiens",
                                       "Balaena_mysticetus"),
                               ordered = TRUE)

low.comp <- ind_low_pivot %>% 
            ggplot(aes(x = species, y = ages, fill = type))+
          geom_col(position = "dodge")+
          theme_bw()+
          scale_x_discrete(limits = c(levels(ind_low_pivot$species)),
                           labels = c("Ursus arctos",
                                      "Acinonyx jubatus",
                                      "Homo sapiens",
                             "Balaena mysticetus"))+
          scale_fill_manual(name = "Ages",
                            labels = c("Mean Age", "Fossil Age",
                                       "Phylogenetic Age"),
                            values = c("#66c2a5",
                                       "#fc8d62",
                                       "#8da0cb"))+
          scale_y_reverse()+
          coord_flip()+
          guides(y = "none",
                 y.sec = "axis")+
          ylab("Time (myrs)")+
          xlab(NULL)+
          #scale_y_discrete(position = "right")+
          mynamestheme+
          theme(legend.position="none")

png("text/supplementary/Figure.low.empirical.ex.png", width = 25,
    height = 15, units = "cm", 
    pointsize = 8, res = 300)

##left
top <- text_grob("Low extinction fraction",
                  family = "serif", size = 13,  face = "bold")
grid.arrange(low.phy, low.comp,
             ncol = 2, top = top)

dev.off()



# Intermediate extinction fraction ----------------------------------------


mammals_int = mammals_emp %>% filter(epsilon == "int")

##df for the selected species
ind_int = left_join(ind_sp_fossil, mammals_int, by = "species")



##all mammals phylo age vs corrected age
int.phy = ggplot(mammals_int, aes(x = Estimated.age, y = mean.age))+
  geom_point(alpha = 0.1)+
  xlim(0,11)+
  ylim(0, 11)+
  geom_point(data = ind_int, aes(x = Estimated.age, y = mean.age), 
             size = 3.5, color = "red", alpha = 0.9)+
  geom_poitn(data = )
  xlab("Phylogenetic age (myr)")+
  ylab("Corrected age (myr)")+
  #ggtitle("int extinction fraction")+
  theme_bw()+
  mynamestheme

##comparison phylogenetic, fossil, and corrected age

##pivoting for plotting
ind_int_pivot = ind_int %>% 
  pivot_longer(cols = ends_with(".age"), names_to = "type",
               values_to = "ages")

##organizing factors
ind_int_pivot$type <- factor(ind_int_pivot$type,
                             levels = c("mean.age",
                                        "fossil.age",
                                        "Estimated.age"),
                             ordered = TRUE)

ind_int_pivot$species = factor(ind_int_pivot$species, 
                               levels = c("Ursus_arctos",
                                          "Acinonyx_jubatus",
                                          "Homo_sapiens",
                                          "Balaena_mysticetus"),
                               ordered = TRUE)

int.comp <- ind_int_pivot %>% 
  ggplot(aes(x = species, y = ages, fill = type))+
  geom_col(position = "dodge")+
  theme_bw()+
  scale_x_discrete(limits = c(levels(ind_int_pivot$species)),
                   labels = c("Ursus arctos",
                              "Acinonyx jubatus",
                              "Homo sapiens",
                              "Balaena mysticetus"))+
  scale_fill_manual(name = "Ages",
                    labels = c("Mean Age", "Fossil Age",
                               "Phylogenetic Age"),
                    values = c("#66c2a5",
                               "#fc8d62",
                               "#8da0cb"))+
  scale_y_reverse()+
  coord_flip()+
  guides(y = "none",
         y.sec = "axis")+
  ylab("Time (myrs)")+
  xlab(NULL)+
  #scale_y_discrete(position = "right")+
  mynamestheme+
  theme(legend.position="none")

png("text/figures/Figure.int.empirical.ex.png", width = 25,
    height = 15, units = "cm", 
    pointsize = 8, res = 300)

##left
top <- text_grob("Intermediate extinction fraction",
                 family = "serif", size = 13,  face = "bold")
grid.arrange(int.phy, int.comp,
             ncol = 2, top = top)

dev.off()



# High extinction fraction ------------------------------------------------



mammals_high = mammals_emp %>% filter(epsilon == "high")

##df for the selected species
ind_high = left_join(ind_sp_fossil, mammals_high, by = "species")



##all mammals phylo age vs corrected age
high.phy = ggplot(mammals_high, aes(x = Estimated.age, y = mean.age))+
  geom_point(alpha = 0.1)+
  xlim(0,11)+
  ylim(0, 11)+
  geom_point(data = ind_high, aes(x = Estimated.age, y = mean.age), 
             size = 3.5, color = "red", alpha = 0.9)+
  xlab("Phylogenetic age (myr)")+
  ylab("Corrected age (myr)")+
  #ggtitle("high extinction fraction")+
  theme_bw()+
  mynamestheme

##comparison phylogenetic, fossil, and corrected age

##pivoting for plotting
ind_high_pivot = ind_high %>% 
  pivot_longer(cols = ends_with(".age"), names_to = "type",
               values_to = "ages")

##organizing factors
ind_high_pivot$type <- factor(ind_high_pivot$type,
                             levels = c("mean.age",
                                        "fossil.age",
                                        "Estimated.age"),
                             ordered = TRUE)

ind_high_pivot$species = factor(ind_high_pivot$species, 
                               levels = c("Ursus_arctos",
                                          "Acinonyx_jubatus",
                                          "Homo_sapiens",
                                          "Balaena_mysticetus"),
                               ordered = TRUE)

high.comp <- ind_high_pivot %>% 
  ggplot(aes(x = species, y = ages, fill = type))+
  geom_col(position = "dodge")+
  theme_bw()+
  scale_x_discrete(limits = c(levels(ind_high_pivot$species)),
                   labels = c("Ursus arctos",
                              "Acinonyx jubatus",
                              "Homo sapiens",
                              "Balaena mysticetus"))+
  scale_fill_manual(name = "Ages",
                    labels = c( "Corrected Age",
                                "Fossil Age",
                                "Phylogenetic Age"),
                    values = c("#66c2a5",
                               "#fc8d62",
                               "#8da0cb"))+
  scale_y_reverse()+
  coord_flip()+
  guides(y = "none",
         y.sec = "axis")+
  ylab("Time (myrs)")+
  xlab(NULL)+
  #scale_y_discrete(position = "right")+
  mynamestheme+
  theme(legend.position="none")



##only the mean.ages changes for joining
ind_low_pivot = ind_low_pivot %>% filter(type == "mean.age")

ind_int_pivot = ind_int_pivot %>% filter(type == "mean.age")

#####joining pivot dataframes
ind_pivot = rbind(ind_low_pivot, ind_int_pivot, ind_high_pivot)


##organizing factors
ind_pivot$type <- factor(ind_pivot$type,
                              levels = c("mean.age",
                                         "fossil.age",
                                         "Estimated.age"),
                              ordered = TRUE)

ind_pivot$species = factor(ind_pivot$species, 
                                levels = c("Ursus_arctos",
                                           "Acinonyx_jubatus",
                                           "Homo_sapiens",
                                           "Balaena_mysticetus"),
                                ordered = TRUE)


ind_pivot$epsilon = factor(ind_pivot$epsilon, 
                           levels = c("low",
                                      "int",
                                      "high"),
                           ordered = TRUE)


png("text/supplementary/Figure.high.empirical.ex.png", width = 25,
    height = 15, units = "cm", 
    pointsize = 8, res = 300)

##left
top <- text_grob("High extinction fraction",
                 family = "serif", size = 13,  face = "bold")
grid.arrange(high.phy, high.comp,
             ncol = 2, top = top)

dev.off()


######all extinctions fractions grouped

##all mammals phylo age vs corrected age

######line plots
line_plot <-  ggplot(mammals_emp, aes(x = Estimated.age, y = mean.age, color = epsilon))+
  geom_point(alpha = 0.9)+
  xlim(0,11)+
  ylim(0, 11)+
  scale_color_manual(values = c("lightgrey", 
                                "#969696",
                                "#525252"),
                     name="Extinction fraction",
                     breaks=c("low", "int", "high"),
                     labels=c("Low",
                              "Intermediate",
                              "High"))+
  geom_point(data = ind_int, aes(x = Estimated.age, y = mean.age), 
             size = 3, color = "red", alpha = 0.9)+
  # geom_point(data = mammals_low, aes(x = Estimated.age, y = mean.age),
  # alpha = 0.1, color = "lightgrey")+
  geom_point(data = ind_low, aes(x = Estimated.age, y = mean.age), 
             size = 3, color = "red", alpha = 0.9)+
  # geom_point(data = mammals_high, aes(x = Estimated.age, y = mean.age),
  #alpha = 0.1, color = "#525252")+
  geom_point(data = ind_high, aes(x = Estimated.age, y = mean.age), 
             size = 3, color = "red", alpha = 0.9)+
  xlab("Phylogenetic age (myr)")+
  ylab("Corrected age (myr)")+
  #ggtitle("int extinction fraction")+
  theme_bw()+
  mynamestheme+ theme(
    legend.position = c(.95, .03),
    legend.justification = c("right", "bottom"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6)
  )

############bar plots

#####starting with high extinction

bar_empirical <-  ggplot(NULL, aes(x = species, y = ages, fill = type))+
  geom_col(data = ind_high_pivot, alpha = 0.9,
           position = "dodge")+
  geom_col(data = ind_int_pivot, alpha = 0.5,
           position = "dodge")+
  geom_col(data = ind_low_pivot, alpha = 0.3,
           position = "dodge")+
  theme_bw()+
  scale_x_discrete(limits = c(levels(ind_high_pivot$species)),
                   labels = c("Ursus arctos",
                              "Acinonyx jubatus",
                              "Homo sapiens",
                              "Balaena mysticetus"))+
  scale_fill_manual(name = "Ages",
                    labels = c( "Corrected Age",
                                "Fossil Age",
                                "Phylogenetic Age"),
                    values = c("#66c2a5",
                               "#fc8d62",
                               "#8da0cb"))+
  scale_y_reverse()+
  coord_flip()+
  guides(y = "none",
         y.sec = "axis")+
  ylab("Time (myrs)")+
  xlab(NULL)+
  #scale_y_discrete(position = "right")+
  mynamestheme+
  theme(legend.position="none")


##png object
png("text/figures/Figure.empirical.all.ex.png", width = 25,
    height = 15, units = "cm", 
    pointsize = 8, res = 300)


grid.arrange(line_plot, bar_empirical,
             ncol = 2)

dev.off()




##########mean age alone###############
png("text/supplementary/Figure.label.ex.png", width = 25,
    height = 15, units = "cm", 
    pointsize = 8, res = 300)

ind_pivot %>% filter(type == "mean.age") %>% 
  ggplot(aes(x = species, y = ages, alpha = epsilon))+
  geom_col(position = "identity", fill = "#66c2a5")+
  scale_x_discrete(limits = c(levels(ind_pivot$species)),
                   labels = c("Ursus arctos",
                              "Acinonyx jubatus",
                              "Homo sapiens",
                              "Balaena mysticetus"))+
  scale_alpha_manual(name = "Extinction fraction",
                     labels = c("Low", "Intermediate",
                                "High"),
                     values = c(0.3, 0.6, 0.9))+
  scale_y_reverse()+
  coord_flip()+
  guides(y = "none",
         y.sec = "axis")+
  ylab("Time (myrs)")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme
#theme(legend.position="none")
dev.off()



