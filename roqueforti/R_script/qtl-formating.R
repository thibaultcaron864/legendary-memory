#Thibault CARON, 09/08/19, final qtl analysis of growth, lipolysis, proteolysis and color
#updated from the 01/07/2022
#last update on the 19/09/2022

library(nlme)
library(xlsx)
library(qtl)
library(qtl2)
library(qtlcharts)
library(FactoMineR)
library(factoextra)
library(corrplot)
#library(linkageMapView) # install me

setwd("/home/thibault/Documents/WORK/CROQ/QTL/qtl/")
load("qtl.RData")

INDELs <- read.csv("genotypes_INDELs.csv", header = T) #genotypes from the new assembly
phenotypes <- read.csv("all_phenotypes_outlierOut.csv", header = T) #all phenotypes without every value inf or sup to 8*(Q75-Q25)
phenotypes_C3 <- phenotypes[which(phenotypes$cross==3),] #only take the phentypes for the cross 3
phenotypes_C3 <- subset(phenotypes_C3, select = -cross) #remove the useless column cross
gbs_all <- read.csv("data_RxN_N.csv", header = T) #gbs
nameCorr <- read.csv("../name_correspondance.tsv", header = F, sep = "\t") #name correspondance
alldata <- merge(phenotypes_C3, INDELs, by = "strain", all = T) #merge phenotypes and genotypes
alldata <- rbind(alldata[390:393,], alldata[388:389,], alldata[1:387,])
#write.csv(alldata, "all_data.csv", row.names = F) #backup
header <- alldata[1:6,] # extract header
#change strain names
alldata$strain<-as.vector(alldata$strain)
nameCorr$V1<-as.vector(nameCorr$V1)
nameCorr$V2<-as.vector(nameCorr$V2)
strainHeader<-alldata$strain[1:6]
alldata$strain<-unlist(lapply(alldata$strain, function(x) nameCorr$V2[match(x, nameCorr$V1)]))
alldata$strain[1:6]<-strainHeader
gbs<-gbs_all[,c(1,60:7592)] #subset to only take genotypes
coln<-colnames(gbs)
gbs<-as.data.frame(t(apply(gbs,1,function(x) as.vector(x))),stringsAsFactors = F)
colnames(gbs)<-coln
gbs<-rbind(gbs[1:2,], c("LCP06133",rep("CC",7533)), c("LCP06136",rep("RR",7533)),gbs[3:389,])
dataQTL <- cbind(alldata[c(6,1,5,2,7:393),c(1,3:215)],gbs[,2:7534])
dataQTL[1:2,1:31] <- "" #no NA allowed
colnames(dataQTL)[1] <- "" #empty first colname
dataQTL[dataQTL=="-"]<-NA
dataQTL[dataQTL=="S"]<-"RR"
dataQTL[dataQTL=="NR"]<-"CC"
#dataQTL2<-dataQTL2[,order(dataQTL2[1,],dataQTL2[2,])] #reorder by increasing contig then position
#does not work, no time to find why, shamely cheating instead
write.csv(t(dataQTL), "data_qtl_GBS-temp.csv", row.names = T) #reorder with office
dataQTL<-read.csv("data_qtl_GBS-temp.csv") #read it back
colnames(dataQTL)[1] <- "" #empty first colname
dataQTL[1:2,1:31] <- "" #no NA allowed

write.csv(dataQTL, "data_qtl_GBS.csv", row.names = F) #backup

#dataQTL <- alldata[c(1:2,5:393),c(1,3:215)]
#dataQTL$strain <- as.vector(dataQTL$strain)
#write.csv(dataQTL, "data_qtl.csv", row.names = F) #backup
