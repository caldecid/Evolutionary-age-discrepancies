
# Fossil vs phylogenetic vs corrected age -----------------------------------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


###calling individual species
ind_sp <- read_delim(file.path(getwd(), pro, q_error, "ind_species.csv"))


##establishing different speciation and extinction scenarios from the diversification rates
###changing the value of the Homo sapiens phylogeny (Rivas-González et al. 2023)


ind_sp <- ind_sp %>% distinct(species, .keep_all = TRUE) %>%  
                  mutate(Estimated.age = replace(Estimated.age,
                                        Estimated.age == 6.175880, 10)) %>% 
                        mutate(high_ext = 0.45 * div,
                            int_ext = 0.3 * div,
                            low_ext = 0.15 * div,
                            f50_samp = 0.5,
                            f25_samp = 0.75,
                            full_samp = 1) %>%
              pivot_longer(cols = ends_with("_ext"), names_to = "ext_scenario",
                           values_to = "mu") %>% 
              pivot_longer(ends_with("samp"), names_to = "sampling",
                           values_to = "fraction_samp") %>% 
              mutate(lambda = 1 - mu)

####Calculating the corrected ages

mean.age <- vector(mode = "numeric",
                           length = nrow(ind_sp))

for(i in 1:nrow(ind_sp)){
  mean.age[[i]] <- meanAgeFunction(lam = ind_sp$lambda[i],
                                           mu = ind_sp$mu[i],
                                           rho = ind_sp$fraction_samp[i],
                                           v = ind_sp$Estimated.age[i])
}


ind_sp <- cbind(ind_sp, mean.age)

ind_sp_pivot <- ind_sp %>% select(-root.age) %>% 
                      pivot_longer(cols = ends_with(".age"), names_to = "type",
                                        values_to = "ages")

ind_sp_pivot$ext_scenario <- factor(ind_sp_pivot$ext_scenario, 
                                    levels = c("low_ext",
                                               "int_ext",
                                               "high_ext"),
                                    ordered = TRUE)

ind_sp_pivot$sampling <- factor(ind_sp_pivot$sampling, 
                                levels = c("full_samp","f25_samp",
                                           "f50_samp"),
                                ordered = TRUE)

ind_sp_pivot$type <- factor(ind_sp_pivot$type,
                            levels = c("mean.age",
                                       "fossil.age",
                                       "Estimated.age"),
                            ordered = TRUE)


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
                      legend.background = element_rect(fill="gray90",
                                                       size=.5, linetype="dotted"),
                      legend.position = "bottom")

##png object
png("text/figures/Figure1.phylo.fossil.mean.full.samp.png", width = 22,
    height = 15,
    units = "cm", 
    pointsize = 8, res = 300)

##remember to change the order of the labels

ind_sp_pivot %>% filter(sampling == "full_samp") %>% 
ggplot(aes(x = species, y = ages, fill = type))+
  geom_col(position = "dodge")+
  theme_bw()+
  scale_x_discrete(limits = c(unique(ind_sp_pivot$species)),
                   labels = c("Carcharhinus obscurus", "Trianenodon obesus",
                              "Vulpes velox", "Acinonyx jubatus", "Ursus arctos",
                              "Homo sapiens", "Balaena mysticetus"))+
  scale_fill_manual(name = "Ages",
                      labels = c("Mean Age", "Fossil Age", "Phylogenetic Age"),
                    values = c("#66c2a5",
                               "#fc8d62",
                               "#8da0cb"))+
  scale_y_reverse()+
  coord_flip()+
  guides(y = "none",
         y.sec = "axis")+
  facet_grid(~ ext_scenario,
             labeller = as_labeller(c("low_ext" = "Low extinction",
                                     "int_ext" = "Intermediate extinction",
                                     "high_ext" = "High extinction"
                                     )))+
  ylab("Time (myrs)")+
  xlab(NULL)+
  #scale_y_discrete(position = "right")+
  mynamestheme

dev.off()  
