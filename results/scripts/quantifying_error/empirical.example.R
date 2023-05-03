
# Fossil vs phylogenetic age ----------------------------------------------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


###calling individual species
ind_sp <- read_delim(file.path(getwd(), pro, q_error, "ind_species.csv"))

##removing duplicate rows
ind.sp <- ind_sp %>% distinct(species, .keep_all = TRUE) %>% 
          select(species, Estimated.age, fossil.age) %>% 
          pivot_longer(cols = ends_with(".age"), names_to = "type",
                       values_to = "ages")

##########Figure
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (11),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (11),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (10)),
                      legend.background = element_rect(fill="gray90",
                                                       size=.5, linetype="dotted"),
                      legend.position = "bottom")

##png object
png("text/figures/Figure1.sp.phylo.vs.fossil.png", width = 15,
    height = 12,
    units = "cm", 
    pointsize = 8, res = 300)

ggplot(ind.sp, aes(x = species, y = ages, fill = type))+
  geom_col(position = "dodge")+
  theme_bw()+
  scale_x_discrete(limits = c(unique(ind.sp$species)),
                   labels = c("Carcharhinus obscurus", "Trianenodon obesus",
                              "Vulpes velox", "Acinonyx jubatus", "Ursus arctos",
                              "Homo sapiens", "Balaena mysticetus"))+
  scale_fill_discrete(name = "Ages",
                      labels = c("Phylogenetic", "Fossil"))+
  scale_y_reverse()+
  coord_flip()+
  guides(y = "none",
         y.sec = "axis")+
  ylab("Time (myrs)")+
  xlab(NULL)+
  #scale_y_discrete(position = "right")+
  mynamestheme

dev.off()  
