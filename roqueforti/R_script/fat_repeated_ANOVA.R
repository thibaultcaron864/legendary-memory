setwd("/Users/Ying-chu/Desktop/R/R_morphology_test/fat/")
library(readr)

fat <- read.table("/Users/Ying-chu/Desktop/R/R_morphology_test/fat/Lipolyse_2.csv", sep=";", header=T)


#這個方法叫做repeated measure ANOVA 如果有重複測驗的點，例如說：以我的資料而言我每個smaple各有四個點(依據天數不同，有重複檢測)
# this method called repeated ANOVA
#這種狀況下就可以使用repeated measure ANOVA
#Repeated Measure ANOVA
library(car)
fatSub<-fat[as.character(fat$Species)=="Psal" | as.character(fat$Species)=="Pnal",] #sub file of protein, for statistic analysis in species
DayFrame<-data.frame(c("seven","fourteen","twentyone","twentyeight"))# this is setting for idata, idata= dataframe
DayFactor<-ordered(as.factor(c("seven","fourteen","twentyone","twentyeight")))#Factors can be number or charaters; this setting is for idesign, idesign= datafactor
Model_fat<-lm(cbind(seven,fourteen,twentyone,twentyeight)~Species+Media%in%Species,data=fatSub)# write your stat model here
Anov_fat<-Anova(Model_fat,idata=DayFrame,idesign=~DayFactor,type="III") #Type III 的意思是 unbalance samplying (Psal/Pnal sample number are different)
summary(Anov_fat,multivariate=F)

#Post Hoc test considering independent tubes
library(lsmeans)
library(multcompView)

# 重新建立新的表單，因為這是重複的統計，而我重複的測驗點有四個，所以每個選項重複四次
# rebuild new data list, as I have four time points, each categories put them four times
# Day= 7, 14, 21, 28; population= sausage/nonsausage; Species="charater" Psal/Pnal; Lysis="number"
fat2<-data.frame(cbind(Day=c(rep("seven",length(fatSub$Media)),rep("fourteen",length(fatSub$Media)),rep("twentyone",length(fatSub$Media)),rep("twentyeight",length(fatSub$Media))),
                       Population=c(as.character(fatSub$Media),as.character(fatSub$Media),as.character(fatSub$Media),as.character(fatSub$Media)),
                       Species=c(as.character(fatSub$Species),as.character(fatSub$Species),as.character(fatSub$Species),as.character(fatSub$Species)),
                       Lysis=as.numeric(c(fatSub$seven,fatSub$fourteen,fatSub$twentyone,fatSub$twentyeight))))

fat2$Lysis<-as.numeric(as.character(fat2$Lysis))

# interaction with population nested in species depends on days
fat_Model2<-lm(Lysis~(Population%in%Species):Day,data=fat2)

# post hoc analysis
# If use TukeyHSD() -> 就會沒有個天數的比較，而是綜合的統計結果，所以才會寫這一段code 單純比較個天數的統計
# can not use TukeyHSD here, otherwise you don't have the comparison between different days, so the following code is for compare in certain time point
leastsquare_fat = lsmeans(fat_Model2, ~(Population%in%Species)|Day, data=fat2,adjust="tukey")
PH_Lipo_all<-cld(leastsquare_fat, alpha=.05, Letters=letters,sort=F,details=T)

PH_Lipo_all
