library(qtl)

#genetic map and QTL analysis script for Proteolysis

setwd("/home/thibault/Documents/WORK/CROQ/+/Charlotte_PERRAULT")

tab <- read.csv(file = "./INDEL_pour_R_proteolysis.csv", sep=";")
cross <- read.cross("csv", file = "./INDEL_pour_R_proteolysis.csv", sep = ";", na.strings = "-", genotypes = c("RR","RC","CC")) #cr?ation du croisement g?n?tique
length(cross$geno)<-6 #remove scaffold >6
summary(cross) #checking phenotypes and genotypes (display nind, nphe, nchr, totmar and nmar)
plot(cross) #g?notypes manquants, petite carte des marqueurs, distrib des ph?nos
plotMap(cross) #marker disposition on genome

#QTL mapping for the phenotype Asym using cross and not the real genetic map
asym.out.em <- scanone(cross, pheno.col = 1, model = "normal") #maximum de vraisemblance EM algorithm (Lander and Botstein 1989)
asym.out.hk <- scanone(cross, pheno.col = 1, method = "hk", model = "normal")#Haley-Knott regression (Haley and Knott 1992).
asym.out.imp <- scanone(cross, pheno.col = 1, method = "imp", model = "normal")# multiple imputation method of Sen and Churchill (2001), need to do sim.geno before ?
par(mfrow=c(3,2))
for(i in 1:6){
  plot(asym.out.hk, asym.out.imp, asym.out.em, chr=i, lty=1, col=c("blue", "red", "black"),ylim=c(0,300))
}
summary(growth.out.em) #maximum lod scores markers per chr
summary(growth.out.em, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(growth.out.hk, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(growth.out.imp, threshold=3)#maximum lod scores markers per chr with lod > 3

#is it significiant ? 
operm.em <- scanone(cross, method="em", model="normal", n.perm=100)
summary(operm.em, alpha=0.05)
summary(asym.out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)



#QTL mapping for the phenotype xmid using cross and not the real genetic map
xmid.out.em<-scanone(cross, pheno.col = 2, model = "normal") #maximum de vraisemblance EM algorithm (Lander and Botstein 1989)
xmid.out.hk<-scanone(cross, pheno.col = 2,method="hk", model = "normal")#Haley-Knott regression (Haley and Knott 1992).
xmid.out.imp<-scanone(cross, pheno.col = 2,method="imp", model = "normal")# multiple imputation method of Sen and Churchill (2001), need to do sim.geno before ?
par(mfrow=c(3,2))
for(i in 1:6){
  plot(xmid.out.hk, xmid.out.imp, xmid.out.em, chr=i, lty=1, col=c("blue", "red", "black"),ylim=c(0,100))
}
summary(xmid.out.em) #maximum lod scores markers per chr
summary(xmid.out.em, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(xmid.out.hk, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(xmid.out.imp, threshold=3)#maximum lod scores markers per chr with lod > 3

#is it significiant ?
operm.em <- scanone(cross, method="em", pheno.col = 2,model="normal", n.perm=100)
summary(operm.em, alpha=0.05)
summary(xmid.out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)


#QTL mapping for the phenotype scal using cross and not the real genetic map
scal.out.em<-scanone(cross, pheno.col = 3, model = "normal") #maximum de vraisemblance EM algorithm (Lander and Botstein 1989)
scal.out.hk<-scanone(cross, pheno.col = 3,method="hk", model = "normal")#Haley-Knott regression (Haley and Knott 1992).
scal.out.imp<-scanone(cross, pheno.col = 3,method="imp", model = "normal")# multiple imputation method of Sen and Churchill (2001), need to do sim.geno before ?
par(mfrow=c(3,2))
for(i in 1:6){
  plot(scal.out.hk, scal.out.imp, scal.out.em, chr=i, lty=1, col=c("blue", "red", "black"),ylim=c(0,15))
}
summary(scal.out.em) #maximum lod scores markers per chr
summary(scal.out.em, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(scal.out.hk, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(scal.out.imp, threshold=3)#maximum lod scores markers per chr with lod > 3

#is it significiant ?
operm.em <- scanone(cross, method="em", pheno.col = 3,model="normal", n.perm=100)
summary(operm.em, alpha=0.05)
summary(scal.out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)
