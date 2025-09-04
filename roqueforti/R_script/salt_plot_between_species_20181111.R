# setting the work folder
setwd("/Users/Ying-chu/Desktop/R/R_morphology_test/salt/")

library(ggplot2)
library(dplyr)
library(gcookbook)
library(grid)
library(gtable)
library(gridExtra)
library(scales)


salt_species_mixday<- read.csv("salt_plot_mixdays_all_apecies.csv",header = T, sep = ";")
salt_Pnal_Pchry <- read.csv("salt_plot_mixdays_Pnal_Pchry.csv",header = T, sep = ",")
salt_Psal_Pols <- read.csv("salt_plot_mixdays_Psal_Pols.csv", header = T, sep = ",")

# group depends on salt concentration
salt_species_mixday  %>% group_by(Salt.concentration) %>% tally()
salt_Pnal_Pchry  %>% group_by(Salt.concentration) %>% tally()
salt_Psal_Pols %>% group_by(Salt.concentration) %>% tally()


# group depends on days
salt_species_mixday  %>% group_by(day) %>% tally()
salt_Pnal_Pchry  %>% group_by(day) %>% tally()
salt_Psal_Pols %>% group_by(day) %>% tally()


# group by different species
salt_species_mixday  %>% group_by(Species) %>% tally()
salt_Pnal_Pchry %>% group_by(Species) %>% tally()
salt_Psal_Pols %>% group_by(Species) %>% tally()


# out put the plot

pd <- position_dodge(.3)

species_mixdays <-ggplot(salt_species_mixday , aes(x=Salt.concentration, y=Colony.area, colour=Species, shape=Species)) +
  geom_errorbar(data =salt_species_mixday , aes(ymin=Colony.area-se, ymax=Colony.area+se), 
                width=.2, size=0.25, colour="black", position = pd) +
  geom_point(data= salt_species_mixday, position=pd, size=3, show.legend=T) +
  scale_colour_manual(values = c("blue4","red4","blue4","red4")) + scale_shape_manual(values=c(1,2,16,17))  + 
  theme_bw()  + expand_limits(x=c(0,2,10,18), y=c(0, 4)) + facet_wrap(~ day)

species_mixdays + scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 4.5)) + scale_x_continuous(breaks=c(0, 2, 10, 18))



### Pnal_Pchry

salt_Pnal_Pchry_1 <-ggplot(salt_Pnal_Pchry , aes(x=as.factor(Salt.concentration), y=Colony.area, colour=Species, shape=Species)) +
  geom_errorbar(data =salt_Pnal_Pchry , aes(ymin=Colony.area-se, ymax=Colony.area+se), 
                width=.2, size=0.25, colour="black", position = pd) +
  geom_point(data= salt_Pnal_Pchry, position=pd, size=3, show.legend=T) +
  scale_colour_manual(values = c("blue4","red4","blue4","red4")) + scale_shape_manual(values=c(16,16,16,16))  + 
  theme_bw() + facet_wrap(~ day) + theme(aspect.ratio=1) 

salt_Pnal_Pchry_1 



### Psal_Pols


salt_Psal_Pols <-ggplot(salt_Psal_Pols , aes(x=as.factor(Salt.concentration), y=Colony.area, colour=Species, shape=Species)) +
  geom_errorbar(data =salt_Psal_Pols , aes(ymin=Colony.area-se, ymax=Colony.area+se), 
                width=.2, size=0.25, colour="black", position = pd) +
  geom_point(data= salt_Psal_Pols, position=pd, size=3, show.legend=T) +
  scale_colour_manual(values = c("blue4","red4","blue4","red4")) + scale_shape_manual(values=c(16,16,16,16))  + 
  theme_bw() + facet_wrap(~ day) + theme(aspect.ratio=1) 

salt_Psal_Pols



