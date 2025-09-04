##ANOVA 

tab<-read.csv(file="C:/Users/Jade/Desktop/Charlotte Stage MAG  2019/INDEL_pour_R_remaniement3.csv", sep=";")

par(mfrow=c(2,2))
Geno<-tab[c(-1,-2),]
Anovaloopgrowth<-function(x){
  Anovgrowth<-lm(as.numeric(as.character(Geno$growth))~as.factor(x))
  return(anova(Anovgrowth)$"Pr(>F)"[1])
}
PValugrowth<-apply(Geno[,4:length(Geno[1,])],2,Anovaloopgrowth)
PValuadjgrowth <- p.adjust(PValugrowth,method = "fdr")
plot(-log10(PValuadjgrowth), main="Growth")
which(p.adjust(PValuadjgrowth,method = "fdr")<0.05)

Anovaloopwal<-function(x){
  Anovwal<-glm(as.numeric(as.character(Geno$WALLABY))~as.factor(x))
  return(anova(Anovwal,test="Chisq")$"Pr(>Chi)"[2])
}
PValuwal<-apply(Geno[,4:length(Geno[1,])],2,Anovaloopwal)
PValuadjwal <- p.adjust(PValuwal,method = "fdr")
plot(-log10(PValuadjwal), main="Wallaby")
which(p.adjust(PValuadjwal,method = "fdr")<0.05)

Anovaloopprot<-function(x){
  Anovprot<-glm(as.numeric(as.character(Geno$proteolysis))~as.factor(x))
  return(anova(Anovprot,test="Chisq")$"Pr(>Chi)"[2])
}
PValuprot<-apply(Geno[,4:length(Geno[1,])],2,Anovaloopprot)
PValuadjprot <- p.adjust(PValuprot,method = "fdr")
plot(-log10(PValuadjprot), main="¨protéolyse : oui ou non")
which(p.adjust(PValuadjprot,method = "fdr")<0.05)

tabprot<-read.csv(file="C:/Users/Jade/Desktop/Charlotte Stage MAG  2019/INDEL_pour_R_proteolysisafterremaniement.csv", sep=";")
Genoprot<-tabprot[c(-1,-2),]

par(mfrow=c(2,2))
Anovalooppasym<-function(x){
  Anovpasym<-lm(as.numeric(as.character(Genoprot$Asym))~as.factor(x))
  return(anova(Anovpasym)$"Pr(>F)"[1])
}
PValupasym<-apply(Genoprot[,4:length(Genoprot[1,])],2,Anovalooppasym)
PValuadjpasym <- p.adjust(PValupasym,method = "fdr")
plot(-log10(PValuadjpasym), main="Proteolysis Asym")
which(p.adjust(PValuadjpasym,method = "fdr")<0.05)

Anovalooppxmid<-function(x){
  Anovpxmid<-lm(as.numeric(as.character(Genoprot$xmid))~as.factor(x))
  return(anova(Anovpxmid)$"Pr(>F)"[1])
}
PValupxmid<-apply(Genoprot[,4:length(Genoprot[1,])],2,Anovalooppxmid)
PValuadjpxmid <- p.adjust(PValupxmid,method = "fdr")
plot(-log10(PValuadjpxmid), main="Proteolysis xmid")
which(p.adjust(PValuadjpxmid,method = "fdr")<0.05)

Anovalooppscal<-function(x){
  Anovpscal<-lm(as.numeric(as.character(Genoprot$scal))~as.factor(x))
  return(anova(Anovpscal)$"Pr(>F)"[1])
}
PValupscal<-apply(Genoprot[,4:length(Genoprot[1,])],2,Anovalooppscal)
PValuadjpscal <- p.adjust(PValupscal,method = "fdr")
plot(-log10(PValuadjpscal), main="Proteolysis scal")
which(p.adjust(PValuadjpscal,method = "fdr")<0.05)

tablip<-read.csv(file="C:/Users/Jade/Desktop/Charlotte Stage MAG  2019/INDEL_pour_R_lipolysisafterremaniement.csv", sep=";")
Genolip<-tablip[c(-1,-2),]

par(mfrow=c(2,2))
Anovalooplasym<-function(x){
  Anovlasym<-lm(as.numeric(as.character(Genolip$Asym))~as.factor(x))
  return(anova(Anovlasym)$"Pr(>F)"[1])
}
PValulasym<-apply(Genolip[,4:length(Genolip[1,])],2,Anovalooplasym)
PValuadjlasym <- p.adjust(PValulasym,method = "fdr")
plot(-log10(PValuadjlasym), main="Lipolysis Asym")
which(p.adjust(PValuadjlasym,method = "fdr")<0.05)

Anovalooplxmid<-function(x){
  Anovlxmid<-lm(as.numeric(as.character(Genolip$xmid))~as.factor(x))
  return(anova(Anovlxmid)$"Pr(>F)"[1])
}
PValulxmid<-apply(Genolip[,4:length(Genolip[1,])],2,Anovalooplxmid)
PValuadjlxmid <- p.adjust(PValulxmid,method = "fdr")
plot(-log10(PValuadjlxmid), main="Lipolysis xmid")
which(p.adjust(PValuadjlxmid,method = "fdr")<0.05)

Anovalooplscal<-function(x){
  Anovlscal<-lm(as.numeric(as.character(Genolip$scal))~as.factor(x))
  return(anova(Anovlscal)$"Pr(>F)"[1])
}
PValulscal<-apply(Genolip[,4:length(Genolip[1,])],2,Anovalooplscal)
PValuadjlscal <- p.adjust(PValulscal,method = "fdr")
plot(-log10(PValuadjlscal), main="Lipolysis scal")
which(p.adjust(PValuadjlscal,method = "fdr")<0.05)
