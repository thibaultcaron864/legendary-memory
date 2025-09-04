library(nlme)
library(xlsx)

setwd("/home/thibault/Documents/WORK/CROQ/+/Charlotte_PERRAULT")

Data <- read.csv("./proteolysis.csv" , sep = ";")
Data2 <- read.csv("./lipolysis.csv", sep = ";")
rownames(Data)<-Data[,1]
rownames(Data2)<-Data2[,1]
CoeffLogisLipo<-NULL
for(i in 1:length(Data2[,1])){
  Day<-c(7,14,21,28)
  Value<-Data2[i,2:5]
  Test<-as.data.frame(cbind(t(Value),Day))
  colnames(Test)[1]<-"Value"
  Test$Value[Test$Value==0]<-0.01
  TOTO<-try(gnls(Value ~ SSlogis(Day, Asym, xmid, scal),data=Test),silent=T)
  if(substr(TOTO,1,1)=="E"){
    CoeffLogisLipo<-rbind(CoeffLogisLipo,c(NA,NA,NA,NA))
  }else
  {
    fm1 <- gnls(Value ~ SSlogis(Day, Asym, xmid, scal),data=Test)
    CoeffLogisLipo<-rbind(CoeffLogisLipo,c(fm1$coefficients,fm1$logLik))
  }
}


pdf("Downloads/FitLipolyse.pdf")
par(mfrow=c(8,7),mar=c(1,1,1,1));for(i in 1:length(Data2[,1])){plot(c(7,14,21,28),c(0,0,0,40),col="white");points(c(7,14,21,28),as.numeric(Data2[i,2:5]),pch=19);lines(c(7,14,21,28),CoeffLogisLipo[i,1]/(1+exp((CoeffLogisLipo[i,2]-c(7,14,21,28))/CoeffLogisLipo[i,3])))}
dev.off()
write.xlsx(CoeffLogisLipo, "coeffLogisLipo.xlsx")

CoeffLogisProteo<-NULL
for(i in 1:length(Data[,1])){
  Day<-c(7,14,21,28)
  Value<-Data[i,2:5]
  Test<-as.data.frame(cbind(t(Value),Day))
  colnames(Test)[1]<-"Value"
  Test$Value[Test$Value==0]<-0.01
  if(sum(Data[i,2:5])==0){
    CoeffLogisProteo<-rbind(CoeffLogisProteo,c(NA,NA,NA,NA))
  } else{
    TOTO<-try(gnls(Value ~ SSlogis(Day, Asym, xmid, scal),data=Test),silent=T)
    if(substr(TOTO,1,1)=="E"){
      CoeffLogisProteo<-rbind(CoeffLogisProteo,c(NA,NA,NA,NA))
    }else
    {
      fm1 <- gnls(Value ~ SSlogis(Day, Asym, xmid, scal),data=Test)
      CoeffLogisProteo<-rbind(CoeffLogisProteo,c(fm1$coefficients,fm1$logLik))
    }
  }
}

pdf("Downloads/FitLipolyse.pdf")
par(mfrow=c(8,7),mar=c(1,1,1,1));for(i in 1:length(Data[,1])){plot(c(7,14,21,28),c(0,0,0,40),col="white");points(c(7,14,21,28),as.numeric(Data[i,2:5]),pch=19);lines(c(7,14,21,28),CoeffLogisProteo[i,1]/(1+exp((CoeffLogisProteo[i,2]-c(7,14,21,28))/CoeffLogisProteo[i,3])))}
dev.off()
write.xlsx(CoeffLogisProteo, "coeffLogisProteo.xlsx")