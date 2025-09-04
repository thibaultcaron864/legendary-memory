library(doBy)
library(car)
library(compositions)
library(MASS)
library(mda)

#sources
#https://rstudio-pubs-static.s3.amazonaws.com/133584_8726e36a0d4348d8918952505956e388.html
#https://www.statmethods.net/stats/anova.html

#working directory
setwd("/home/thibault/Documents/WORK/CROQ/METABARCODING")

#loading data
data<-read.csv("Synthese_totale_CROQ2017.csv", header=TRUE,na.strings="NA",dec=",",sep="\t")

##metabarcoding data
#datasets
metab9_20=acomp(data[37:108,20:31])
metab9=acomp(data[37:72,20:31])
metab20=acomp(data[73:108,20:31])
#linear models
wave=data$Wave[37:108]
fabrication=data$Fabrication[37:108]
population=data$Population[37:108]
stage=data$Stage[37:108]
lm9_20_wfps_sum <-lm(ilr(metab9_20)~ wave + fabrication + population + stage)
lm9_20_wfps_prod <-lm(ilr(metab9_20)~ wave*fabrication*population*stage)
wave=data$Wave[37:72]
fabrication=data$Fabrication[37:72]
population=data$Population[37:72]
lm9_wfp_sum <-lm(ilr(metab9)~ wave + fabrication + population)
lm9_wfp_prod <-lm(ilr(metab9)~ wave*fabrication*population)
wave=data$Wave[73:108]
fabrication=data$Fabrication[73:108]
population=data$Population[73:108]
lm20_wfp_sum <-lm(ilr(metab20)~ wave + fabrication + population)
lm20_wfp_prod <-lm(ilr(metab20)~ wave*fabrication*population)

#anovas
Anova(lm9_20_wfps_sum)
Anova(lm9_20_wfps_prod)
Anova(lm9_wfp_sum)
Anova(lm9_wfp_prod)
Anova(lm20_wfp_sum)
Anova(lm20_wfp_prod)

#la fabrication semble meilleure que la production car cette dernièremasque l'effet vague
#a priori les produits sont plus intéressants mais les valeurs diffèrent des sommes pour les facteurs simples ?
#ilr est choisi car valeurs infinie ou 0 présents avec clr et alr
#tester les différentes statistiques de test : test.statistic = "Wilks" / "Hotelling-Lawley" / "Roy" / "Pillai"
#vérifier la distribution normale et test de Levene

##microbial and physicochemical data
ubio <- data[,c(1,3,6,7,9:19)] #microbial counts columns
pc <- data[,c(1,3,6,7,32:51)] #physicochemical columns

#ANOVA for each variable at each stage
subsett <- ubio #choose your subset
attach(subsett)
for(i in 5:ncol(subsett)){ #for each column except wave, fab, pop and stage
  for(j in levels(as.factor(Stage))){ #for each stage
    tmp <- subsett[Stage==j,]
    x <- tmp[,i] #select x
    if (grepl("NA", paste(c(x), collapse='')) | max(x)==0){ #if there is NA or all values are zero
    cat(paste(colnames(subsett)[i], j, "NA", "NA", sep = "\t"), sep = "\n") #write variable, stage, p-val and hypothesis
    }
    else {
      fit <- aov(x ~ Wave + Fabrication + Population, data=tmp)
      if(drop1(fit,~.,test="F")[2,6]<0.05){sig="H1"}else{sig="H0"} #set H0 or H1 with alpha=5%
      cat(paste(colnames(subsett)[i], j, drop1(fit,~.,test="F")[2,6], sig, sep="\t"), sep="\n") #write variable, stage, p-val and hypothesis
    }
  }
}
#MG + MAS at 180 days found significantly different -> let's see how
mean(Data[Data$Stage=="180", Data$Population=="Roquefort"])