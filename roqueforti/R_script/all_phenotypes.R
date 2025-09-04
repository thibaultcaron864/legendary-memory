#Thibault CARON, july 2022, concatenate all phenotypes together
#please titi, update this when you decide to change the value of a phenotype

setwd("/home/thibault/Documents/WORK/CROQ/QTL/phenotypes")

#load data (keep strain as independant vector)
growth <- read.csv("growth/reglin/reglin.csv")
lipo <- read.csv("lipoProteo/2_regression/lipoReg_all.csv")
proteo <- read.csv("lipoProteo/2_regression/proteoReg_all.csv")
lipoProteoMax <- read.csv("lipoProteo/2_regression/lipoProteo_max.csv")
#mycotoxinsArea <- read.csv("mycotoxins/mycotoxins_area_mean_all.csv")
#mycotoxinsConc <- read.csv("mycotoxins/mycotoxins_conc_mean_all.csv")
mycotoxins <- read.csv("mycotoxins/mycotoxins_conc_mean_all.csv")
colors <- read.csv("color/HSB_Lab_RGB/color_full_means.csv")

#just paste the max after the regression
lipo <- cbind(lipo, lipoProteoMax$lipoMax)
proteo <- cbind(proteo, lipoProteoMax$proteoMax)

#remove useless measures
lipo <- lipo[,c(1:5,7,10)]
proteo <- proteo[,c(1:5,7,10)]
growth <- growth[,c(1:3)]

#change colnames
colnames(lipo) <- c("strain","cross","lipoAsym","lipoXmid","lipoScal","lipoSlope","lipoMax")
colnames(proteo)<-c("strain","cross","proteoAsym","proteoXmid","proteoScal","proteoSlope","proteoMax")
colnames(growth)<-c("strain","cross","growthSlope")
#colnames(mycotoxinsArea)<-c("strain","cross","PRtoxinArea","MPAArea","AndrastinAArea","EremefortinBArea","IsofumigaclavineAArea","EremefortinAArea","RoquefortinCArea")
#colnames(mycotoxinsConc)<-c("strain","cross","PRtoxinConc","MPAConc","AndrastinAConc","EremefortinBConc","IsofumigaclavineAConc","EremefortinAConc","RoquefortinCConc")
colnames(mycotoxins)<-c("strain","cross","PRtoxin","MPA","AndrastinA","EremefortinB","IsofumigaclavineA","EremefortinA","RoquefortinC")
colnames(colors)[1]<-"strain"

## merge
#lipo & proteo
all <- merge.data.frame(lipo, proteo, by = "strain")
all <- subset(all, select = -cross.y)
colnames(all)[2] <- "cross"
# + colors
#remove duplicated measures in color per cross
colors <- colors[-c(which(colors$strain=="LCP06133-1"), which(colors$strain=="LCP06133-2"), which(colors$strain=="LCP06136-1"), 
                    which(colors$strain=="LCP06136-2"), which(colors$strain=="LCP06039")),]
all2 <- merge.data.frame(all, colors, by = "strain", all.x = T, all.y = T)
all2$cross <- pmin(all2$cross.x, all2$cross.y, na.rm = T)
all2 <- subset(all2, select = -c(cross.x, cross.y))
# + growth
all3 <- merge(all2, growth, by = "strain", all = T)
all3$cross <- pmin(all3$cross.x, all3$cross.y, na.rm = T)
all3 <- subset(all3, select = -c(cross.x, cross.y))
# + mycotoxins
all4 <- merge(all3, mycotoxins, by = "strain", all = T)
all4$cross <- pmin(all4$cross.x, all4$cross.y, na.rm = T)
all4 <- subset(all4, select = -c(cross.x, cross.y))
#all5 <- merge(all4, mycotoxinsConc, by = "strain", all = T)
#all5$cross <- pmin(all5$cross.x, all5$cross.y, na.rm = T)
#all5 <- subset(all5, select = -c(cross.x, cross.y))

#allPheno <- all5
allPheno <- all4
rm(all, all2, all3, all4)

#duplicate line for multicross parents
LCP06136_C5 <- c(list("strain" = "LCP06136"), allPheno[allPheno$strain=="LCP06136",3:length(allPheno[1,])-1],list("cross" = 5))
LCP06136_C8 <- c(list("strain" = "LCP06136"), allPheno[allPheno$strain=="LCP06136",3:length(allPheno[1,])-1],list("cross" = 8))
LCP06133_C13 <- c(list("strain" = "LCP06133"), allPheno[allPheno$strain=="LCP06133",3:length(allPheno[1,])-1],list("cross" = 13))
LCP06133_C16 <- c(list("strain" = "LCP06133"), allPheno[allPheno$strain=="LCP06133",3:length(allPheno[1,])-1],list("cross" = 16))
allPheno <- rbind(allPheno, LCP06136_C5, LCP06136_C8, LCP06133_C13, LCP06133_C16)

#create a population vector
allPheno$population <- allPheno$strain
allPheno$population[grepl("I|J.*", allPheno$strain)] <- "progeny" #progeny
#parents
allPheno$population[allPheno$strain=="LCP06133"] <- "non-Roquefort"
allPheno$population[allPheno$strain=="LCP06136"] <- "Roquefort"
allPheno$population[allPheno$strain=="LCP06043"] <- "silage"
allPheno$population[grepl("LCP06039", allPheno$strain)] <- "wood"
allPheno$population[allPheno$strain=="LCP06037"] <- "wood"
allPheno$population[allPheno$strain=="LCP06059"] <- "silage"
#other
allPheno$population[allPheno$strain=="BRFM3019"] <- "Roquefort"
allPheno$population[allPheno$strain=="BRFM3020"] <- "Roquefort"
allPheno$population[allPheno$strain=="BRFM3021"] <- "non-Roquefort"
allPheno$population[allPheno$strain=="BRFM3022"] <- "non-Roquefort"
allPheno$population[allPheno$strain=="BRFM3023"] <- "non-Roquefort"
allPheno$population[allPheno$strain=="BRFM3024"] <- "Roquefort"
allPheno$population[allPheno$strain=="LCP03914"] <- "silage"
allPheno$population[allPheno$strain=="LCP06040"] <- "wood"
allPheno$population[allPheno$strain=="LCP06064"] <- "wood"
allPheno$population[allPheno$strain=="LCP05885"] <- "silage"

#reorder
allPheno <- subset(allPheno, select = c(strain, cross, population, lipoAsym:RoquefortinC))
write.csv(allPheno, "all_phenotypes.csv", row.names = F)

## remove outliers
allPhenoProper <- allPheno
for (i in 4:length(allPhenoProper[1,])){
  dataMean <- mean(allPhenoProper[,i], na.rm = T)
  dataSd <- sd(allPhenoProper[,i], na.rm = T)
  q25 <- quantile(allPhenoProper[,i], probs = 0.25, na.rm = T)
  q75 <- quantile(allPhenoProper[,i], probs = 0.75, na.rm = T)
  iqr <- q75 - q25
  #cutOff <- dataSd*10
  cutOff <- iqr*8
  #allPhenoProper[,i][allPhenoProper>(dataMean+cutOff)]<-NA
  #allPhenoProper[,i][allPhenoProper<(dataMean-cutOff)]<-NA
  allPhenoProper[,i][allPhenoProper[,i]<(q25-cutOff)] <- NA
  allPhenoProper[,i][allPhenoProper[,i]>(q75+cutOff)] <- NA
}

write.csv(allPhenoProper, "all_phenotypes_outlierOut.csv", row.names = F)

