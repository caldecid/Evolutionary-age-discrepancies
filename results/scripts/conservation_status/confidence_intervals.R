
# Evaluating confidence intervals -----------------------------------------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))

#####ages total from the conservation status section
load(file = file.path(pro, st_sp, "ages.total.join.RData"))


ages.total.join %>% mutate(CI_up_sc = CI_up/root.age,
                           CI_low_sc = CI_low/root.age,
                           max_sc = max.prob.age/root.age) %>% 
  select(CI_up_sc, CI_low_sc, max_sc, extinction) %>%
    pivot_longer(!extinction, names_to = "age.type", values_to = "age.value") %>% 
    ggplot(aes(x = age.value, fill = age.type))+
    geom_histogram(alpha = 0.2, position = "identity")+
    facet_wrap(~extinction)


ages.total.join %>% select(CI_up, CI_low, max.prob.age, mean.age, root.age, extinction) %>% 
  ggplot(aes(x = max.prob.age/root.age))+
  geom_histogram()+
  facet_wrap(~extinction)

ages.total.join %>% 
  ggplot(aes((CI_up - CI_low)/root.age))+
  geom_histogram()+
  facet_wrap(~extinction)+
  theme_bw()
