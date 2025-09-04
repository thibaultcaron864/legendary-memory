library(qtl)

#genetic map and QTL analysis script

setwd("C:/Users/Jade/Desktop/Charlotte Stage MAG  2019")

tab<-read.csv(file="C:/Users/Jade/Desktop/Charlotte Stage MAG  2019/INDEL_pour_R_remaniement3.csv", sep=";")
cross<-read.cross("csv",file="C:/Users/Jade/Desktop/Charlotte Stage MAG  2019/INDEL_pour_R_remaniement3.csv", sep=";",na.strings="-", genotypes=c("RR","RC","CC")) #création du croisement génétique
length(cross$geno)<-6
summary(cross) #checking phenotypes and genotypes (display nind, nphe, nchr, totmar and nmar)
plot(cross) #génotypes manquants, petite carte des marqueurs, distrib des phénos
plotMap(cross) #marker disposition on genome

geneticmap <- est.map(cross, error.prob=0.05, verbose=TRUE) #genetic map
plotMap(geneticmap) #huge distance in chrom 1 and 3 -> hotspots of recombination ! check what's in there
recomb<- est.rf(cross) #get recombination faction between pairs of adjecent markers -> LOD score
plotRF(recomb) #LOD score map
plotMissing(recomb) #missing genotypes

geneticmap<- replace.map(cross,geneticmap) #replace the genetic map of the cross
growth <- calc.errorlod(growth, error.prob=0.05) #calcul des erreurs LOD scores
top.errorlod(growth)
plotErrorlod(growth) #get LOD score errors



#QTL mapping for the phenotype "growth" using cross and not the real genetic map
#growth <- calc.genoprob(growth, step=100, err=0.05)#QTL mapping (on peut garder hyper)
growth.out.em<-scanone(geneticmap,pheno.col = 1, model = "normal") #maximum de vraisemblance EM algorithm (Lander and Botstein 1989)
growth.out.hk<-scanone(geneticmap,pheno.col = 1,method="hk", model = "normal")#Haley-Knott regression (Haley and Knott 1992).
growth.out.imp<-scanone(geneticmap,pheno.col = 1,method="imp", model = "normal")# multiple imputation method of Sen and Churchill (2001), need to do sim.geno before ?
par(mfrow=c(3,2))
for(i in 1:6){
  plot(growth.out.hk, growth.out.imp, growth.out.em, chr=i, lty=1, col=c("blue", "red", "black"),ylim=c(0,max(growth.out.hk$lod, growth.out.imp$lod, growth.out.em$lod)))
}
summary(growth.out.em) #maximum lod scores markers per chr
summary(growth.out.em, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(growth.out.hk, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(growth.out.imp, threshold=3)#maximum lod scores markers per chr with lod > 3

#is it significiant ? 
operm.em <- scanone(geneticmap, method="em", model="normal", n.perm=1000)
summary(operm.em, alpha=0.05)
summary(growth.out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)
operm.imp <- scanone(geneticmap, method="imp", model="normal", n.perm=1000)
summary(operm.em, alpha=0.05)
summary(growth.out.imp, perms=operm.imp, alpha=0.05, pvalues=TRUE)


#QTL mapping for the phenotype "WALLABY" using the real genetic map
#LOD score for WALLABY with cross map (not the real genetic map)
wal.out.em<-scanone(geneticmap, pheno.col = 3, model = "binary") #maximum de vraisemblance EM algorithm (Lander and Botstein 1989)
wal.out.hk<-scanone(geneticmap, pheno.col = 3,method="hk", model = "binary")#Haley-Knott regression (Haley and Knott 1992).
wal.out.imp<-scanone(geneticmap, pheno.col = 3,method="imp", model = "binary")# multiple imputation method of Sen and Churchill (2001), need to do sim.geno before ?
par(mfrow=c(3,2))
for(i in 1:6){
  plot(wal.out.hk, wal.out.imp, wal.out.em, chr=i, lty=1, col=c("blue", "red", "black"),ylim=c(0,max(wal.out.hk$lod, wal.out.imp$lod, wal.out.em$lod)))
}
summary(wal.out.em) #maximum lod scores markers per chr
summary(wal.out.em, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(wal.out.hk, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(wal.out.imp, threshold=3)#maximum lod scores markers per chr with lod > 3

#is it significiant ?
operm.em <- scanone(geneticmap, method="em", pheno.col = 2,model="binary", n.perm=1000)
summary(operm.em, alpha=0.05)
summary(wal.out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)



#QTL mapping for the phenotype "proteolysis" (proteolyse or not) using cross and not the real genetic map
prot.out.em<-scanone(cross, pheno.col = 3, model = "binary") #maximum de vraisemblance EM algorithm (Lander and Botstein 1989)
prot.out.hk<-scanone(cross, pheno.col = 3,method="hk", model = "binary")#Haley-Knott regression (Haley and Knott 1992).
prot.out.imp<-scanone(cross, pheno.col = 3,method="imp", model = "binary")# multiple imputation method of Sen and Churchill (2001), need to do sim.geno before ?
par(mfrow=c(2,2))
for(i in 1:4){
  plot(prot.out.hk, prot.out.imp, prot.out.em, chr=i, lty=1, col=c("blue", "red", "black"),ylim=c(0,max(prot.out.hk$lod, prot.out.imp$lod, prot.out.em$lod)))
}
summary(prot.out.em) #maximum lod scores markers per chr
summary(prot.out.em, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(prot.out.hk, threshold=3)#maximum lod scores markers per chr with lod > 3
summary(prot.out.imp, threshold=3)#maximum lod scores markers per chr with lod > 3

#is it significiant ?
operm.em <- scanone(cross, method="em", pheno.col = 3,model="binary", n.perm=100)
summary(operm.em, alpha=0.05)
summary(prot.out.em, perms=operm.em, alpha=0.05, pvalues=TRUE)


cross2 <- sim.geno(cross2, step=10, n.draws=64, err=0.01)

top.errorlod(cross2)
plotGeno(hyper, chr=1) #what is the blue cross ?
help(ploplotGeno(cross2, chr=2, ind=c(24:34, 71:81))
plotGeno(cross2, chr=3, ind=c(24:34, 71:81))
plotGeno(cross2, chr=4, ind=c(24:34, 71:81))
plotGeno(cross2, chr=5, ind=c(24:34, 71:81))
plotGeno(cross2, chr=6, ind=c(24:34, 71:81))
