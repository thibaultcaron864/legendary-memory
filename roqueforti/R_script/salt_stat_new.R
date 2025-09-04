# setting the work folder
setwd("/Users/Ying-chu/Desktop/R/R_morphology_test/salt/")

# read the library
library(readr)
library(tidyr)


# read input 
first_all <- read.csv("/home/thibault/Desktop/R/R_morphology_test/salt/first_all.csv", header=T)
second_all <- read.csv("/home/thibault/Desktop/R/R_morphology_test/salt/second_all.csv", header = T)
third_all <- read.csv("/home/thibault/Desktop/R/R_morphology_test/salt/third_all.csv", header = T)

# read input 
first_Psal_Pols_Pnal_Pchry <- read.csv("/home/thibault/Desktop/R/R_morphology_test/salt/first_Psal_Pols_Pnal_Pchry.csv", header=T)
second_Psal_Pols_Pnal_Pchry <- read.csv("/home/thibault/Desktop/R/R_morphology_test/salt/second_Psal_Pols_Pnal_Pchry.csv", header = T)
third_Psal_Pols_Pnal_Pchry <- read.csv("/home/thibault/Desktop/R/R_morphology_test/salt/third_Psal_Pols_Pnal_Pchry.csv", header = T)

# attached the data
attach(first_all)
detach(first_all)

attach(second_all)
detach(second_all)

attach(third_all)
detach(third_all)

# attached the data
attach(first_Psal_Pols_Pnal_Pchry)
detach(first_Psal_Pols_Pnal_Pchry)

attach(second_Psal_Pols_Pnal_Pchry)
detach(second_Psal_Pols_Pnal_Pchry)

attach(third_Psal_Pols_Pnal_Pchry)
detach(third_Psal_Pols_Pnal_Pchry)

# write the model

Modelanova_yl_1log <-lm(log(colonyarea)~Species+ordered(Medium)+Species*ordered(Medium)+ Sausage.NonSausage%in%Species)

Modelanova_yl_2log <-lm(log(colonyarea)~Species+ordered(Medium)+Species*ordered(Medium)+ Sausage.NonSausage%in%Species)

Modelanova_yl_3log <-lm(log(colonyarea)~Species+ordered(Medium)+Species*ordered(Medium)+ Sausage.NonSausage%in%Species)


hist(residuals(Modelanova_yl_1log))
plot(Modelanova_yl_1log)
     
     
# use ANOVA

anova(Modelanova_yl_1log)

anova(Modelanova_yl_2log)

anova(Modelanova_yl_3log)

# post hoc test 

TukeyHSD(aov(Modelanova_yl_1log))

TukeyHSD(aov(Modelanova_yl_3log))

# consider medium as number

# show "Medium" column
Medium

# setting Medium -> MediumNum
MediumNum<-Medium

# setting "0%" as "0", "2%" as "2" ...so on

MediumNum[MediumNum=="0%"]<-0
MediumNum[MediumNum=="2%"]<-20
MediumNum[MediumNum==20]<-2
MediumNum[MediumNum=="10%"]<-10
MediumNum[MediumNum=="18%"]<-18

# show the result of Medium number
MediumNum

# setting as numbers by using as.numeric function
MediumNum<-as.numeric(MediumNum)

MediumNumSq<-(MediumNum-2)^2

# re-write new model
Modelanova_yl_1log<-lm(log(colonyarea)~Species+MediumNumSq+Species*MediumNumSq+`Sausage/NonSausage`%in%Species)

# Do ANOVA test
Anova(Modelanova_yl_1log)


# calculate the 97.5% - 2.5% for drawing figures
cbind(confint(lm(log(colonyarea)~ 0+Medium:Sausag.NonSausage%in%Species)),summary(lm(log(colonyarea)~ 0+Medium:Sausage.NonSausage%in%Species))$coefficients[,1])


cbind(confint(lm(log(colonyarea)~ 0+Medium:Species)),summary(lm(log(colonyarea)~ 0+Medium:Species))$coefficients[,1])



write.table(cbind(confint(lm(log(colonyarea)~ 0+Medium:Species)),summary(lm(log(colonyarea)~ 0+Medium:Species))$coefficients[,1]),file="Salt_Medium_Species_av_confint.txt",sep="\t",quote=F)









?confint
?coefficients







