#Thibault CARON, 09/08/19, final qtl analysis of growth, lipolysis, proteolysis and color
#updated from the 01/07/2022
#last update on the 19/09/2022

library(nlme)
library(xlsx)
library(qtl)
library(qtl2)
library(qtlcharts)
#library(linkageMapView) # install me

setwd("/home/thibault/Documents/WORK/CROQ/QTL/qtl/")
load("qtl.RData")

#load datasets
phenotypes <- read.csv("all_phenotypes_outlierOut.csv", header = T) #all phenotypes without every value inf or sup to 8*(Q75-Q25)
phenotypes_C3 <- phenotypes[which(phenotypes$cross==3),] #only take the phentypes for the cross 3
dataQTL<-read.csv("data_qtl_GBS.csv", header = T) 

#### make a first genetic map ####
# keep the rearranged markers duplicated at the two parental position
#cross <- read.cross("csv", file = "./data_carto_diplo_dup.csv", sep = ",", na.string = "NA", genotypes = c("RR", "CC"), estimate.map = FALSE, crosstype = "bc") #with re-arranged region duplicated
cross <- read.cross("csv", file = "./data_carto_diplo_dup.csv", sep = ",", na.string = "NA", genotypes = c("RR", "CC"), estimate.map = FALSE, crosstype = "bc") #with re-arranged region duplicated


## looking for missing values
plotMissing(cross)
par(mfrow = c(1,2), las = 1)
plot(ntyped(cross), ylab = "No. types markers", main = "No. genotypes by individual") #4 individuals with a lot of NA
plot(ntyped(cross, "mar"), ylab = "No. typed individuals", main = "No. genotypes by marker") # 2 markers with a lot of NA
#drop individuals (four here)
which(nmissing(cross, what = "ind")>50) #185=I0217, #232=I1453, 256=I1477, 282=I1504
cross <- subset(cross, ind = ntyped(cross)>50) #remove them
nmissing(cross, what = "mar") #Wallaby = 12
#check
dev.off()
plot(ntyped(cross), ylab = "No. types markers", main = "No. genotypes by individual") #ok now
#drop markers (none here)
nt.bymar <- ntyped(cross, "mar")
todrop <- names(nt.bymar[nt.bymar <200])
cross <- drop.markers(cross, todrop)

## identify duplicate individuals
cg <- comparegeno(cross) #matrix of comparison
hist(cg[lower.tri(cg)], breaks = seq(0, 1, len = 101), xlab = "No. matching genotypes", 
     main = "Markers with pairs of individuals matching genotypes")
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.9, arr = TRUE) #find pairs of individuals with 90% of matching genotypes
wh <- wh[wh[,1] < wh[,2],]
g <- pull.geno(cross)
cross <- subset(cross, ind=-wh[,2]) #remove
#check-it
cg <- comparegeno(cross) 
hist(cg[lower.tri(cg)], breaks = seq(0, 1, len = 101), xlab = "No. matching genotypes", 
     main = "Markers with pairs of individuals matching genotypes")
rug(cg[lower.tri(cg)])
print(dup <- findDupMarkers(cross, exact.only=FALSE)) 
#INDEL154R=31R, INDEL58C=31C > to keep
#INDEL34=37, INDEL84=112, INDEL191=INDEL45, INDEL171=INDEL64, INDEL134=INDEL185, INDEL103=INDEL135=INDEL110 > to remove
cross <- drop.markers(cross, unlist(unname(dup))[3:8]) #only remove the last 6 ones

## markers with distorted segregation patterns
gt <- geno.table(cross)
print(distort <- gt[gt$P.value < 0.05/totmar(cross),]) # 95 markers with distorted segregation pattern...
todrop <- rownames(gt[gt$P.value < 1e-10,]) #44 with 1e-10
#cross <- drop.markers(cross, todrop) #keep them as we can expect such pattern from this cross

## individuals’ genotype frequencies
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq) )
par(mfrow = c(1,3), las = 1)
for(i in 1:3){
  plot(gfreq[i,], ylab = "Genotype frequency", main = c("RR", "CC", "RC")[i], ylim = c(0,1))
}
dev.off()

## recombination fraction
cross <- est.rf(cross)
rf <- pull.rf(cross)
lod <- pull.rf(cross, what = "lod")
plot(as.numeric(rf), as.numeric(lod), xlab = "Recombination fraction", ylab = "LOD score")

## form linkage groups
# max.rf = 0.20, min.lod = 6 -> 15 linkage groups -> can do less
# max.rf = 0.28, min.lod = 6 -> 8  linkage groups -> quite ok -> choose that one 
# max.rf = 0.28, min.lod = 6 -> 6  linkage groups -> bullshit

lg <- formLinkageGroups(cross, max.rf = 0.28, min.lod = 6) 
table(lg[,2]) # eight linkage group with cheesyTer alone (absent in one parent at located at the end of chromosome 3)
cross <- formLinkageGroups(cross, max.rf = 0.28, min.lod = 6, reorgMarkers = TRUE)
plotRF(cross, alternate.chrid = TRUE)

## print marker names for each lg
for (i in 1:length(cross$geno)){
  map <- as.data.frame(cross$geno[[i]][1])
  print(paste("lg",i))
  if(length(names(map))>1){
  print(sapply(strsplit(names(map), "\\."), "[[", 2))}
  else{
    print(names(map))}
}

crossProper <- cross #backup

# and finally estimate a genetic map from the physical position
MLE <- 0.0025 # afterward estimated genotyping error rate
map <- est.map(cross, error.prob = MLE, map.function = "kosambi", verbose = TRUE) #estimate map from physical position
cross <- replace.map(cross, map) #replace it

crossProperMap<- cross

##### checking to reorder markers ####
## check for each linkage group except the 8 (only HGT cheesyTer alone)
# # first method with the number of crossover (skipped)
# rip <- list()
# for (i in rev(1:(length(cross$geno)-1))){rip[[i]]<-summary(ripple(cross,chr=i,window=7))}
# save.image("qtl.RData") # save environment
# # have a look at the rip to see for each linkage group if there is a better alternative (with less obligate crossing-over)
# # replace better combinations
# crossXO <- cross
# for(i in 1:length(rip)){if(rip[[i]][2,dim(rip[[i]])[2]]<rip[[i]][1,dim(rip[[i]])[2]]){crossXO<-switch.order(crossXO,chr=i,rip[[i]][2,1:(dim(rip[[i]])[2]-1)],error.prob=0.0001)}}

# second method more meticulous with likelihood (window = 4 takes some time)
riplik <- list()
for (i in rev(1:(length(cross$geno)-1))){print(paste("linkage group",i));riplik[[i]]<-summary(ripple(cross,chr=i,window=4,method="likelihood",error.prob=MLE))}
save.image("qtl.RData") # save environment
# have a look at the riplik to see for each linkage group if there is a better alternative (with >= +3 LOD  and chrlen <= initial combination)
# replace better combinations
LOD_threshold <- 3
for (i in 1:length(riplik)){if(riplik[[i]][2,dim(riplik[[i]])[2]-1]>LOD_threshold && riplik[[i]][2,dim(riplik[[i]])[2]]<riplik[[i]][1,dim(riplik[[i]])[2]]){cross<-switch.order(cross,chr=i,riplik[[i]][2,1:(dim(riplik[[i]])[2]-2)],error.prob=MLE)}}
# plot and compare genetic maps
svg("genetic_map_LIK.svg", width = 20, height = 30)
plotMap(cross, show.marker.names = TRUE, xlab = "linkage group")
dev.off()

crossProperMapLIK <- cross #backup

# important distances on a linkage group may be due to problematic markers (either high error rate or large number of individuals)
# drop one marker at a time to detect the faulty one
dropone <- droponemarker(cross, error.prob = MLE)
plot(dropone, lod = 1, ylim = c(-100,50));abline(h=0) # one bad marker on lg4
plot(dropone, lod = 2, ylab = "decrease in chromosome length") # on marker 
summary(dropone, lod.column=2) # droping INDEL89 + INDEL INDEL7

# dropping bad markers and re-estimate the genetic map
badmar <- rownames(summary(dropone, lod.column=2))[3:4]
cross <- drop.markers(cross, badmar)
newmap <- est.map(cross, error.prob = MLE)
cross <- replace.map(cross, newmap)
summaryMap(cross)

crossProperMapLIKdrop <- cross # backup

# look for problem individuals
plot(countXO(cross), ylab = "Number of crossovers");abline(h=40) # five individuals showing more than 40 XO
cross <- subset(cross, ind=(countXO(cross) < 40)) # remove them
plot(countXO(cross), ylab = "Number of crossovers") # check it

crossProperMapLIKdropProper <- cross # backup

# reorderlinkage group
crosslg1<-cross$geno$`1`
crosslg2<-cross$geno$`2`
crosslg3<-cross$geno$`3`
crosslg4<-cross$geno$`4`
crosslg5<-cross$geno$`5`
crosslg6<-cross$geno$`6`
crosslg7<-cross$geno$`7`
cross$geno$`1`<-crosslg5
cross$geno$`2`<-crosslg3
cross$geno$`3`<-crosslg6
cross$geno$`4`<-crosslg2
cross$geno$`5`<-crosslg1
cross$geno$`6`<-crosslg7
cross$geno$`7`<-crosslg4
rm(crosslg1,crosslg2,crosslg3,crosslg4,crosslg5,crosslg6,crosslg7)

crossProperMapLIKdropProperReordered <- cross # backup

# check marker order again (with window = 3 this time)
riplik2 <- list()
for (i in rev(1:(length(cross$geno)-1))){print(paste("linkage group",i));riplik2[[i]]<-summary(ripple(cross,chr=i,window=3,method="likelihood",error.prob=MLE))}
save.image("qtl.RData") # save environment
for (i in 1:length(riplik2)){if(riplik2[[i]][2,dim(riplik2[[i]])[2]-1]>LOD_threshold && riplik2[[i]][2,dim(riplik2[[i]])[2]]<riplik2[[i]][1,dim(riplik2[[i]])[2]]){cross<-switch.order(cross,chr=i,riplik2[[i]][2,1:(dim(riplik2[[i]])[2]-2)],error.prob=MLE)}}
svg("genetic_map_LIK2.svg", width = 20, height = 30)
plotMap(cross, show.marker.names = TRUE, xlab = "linkage group")
dev.off()

crossProperMapLIKdropProperReorderedLIK <- cross # backup

## re-estimate the genetic map
newmap <- est.map(cross, error.prob=MLE)
cross <- replace.map(cross, newmap)
summaryMap(cross)

crossProperMapLIKdropProperReorderedLIKMap <- cross # backup

# estimate genotyping error rate
loglik <- err <- c(0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02)
for(i in seq(along=err)) {
  cat(i, "of", length(err), "\n")
  tempmap <- est.map(cross, error.prob=err[i])
  loglik[i] <- sum(unlist(sapply(tempmap, attr, "loglik")))
}
lod <- (loglik - max(loglik))/log(10)
plot(err, lod, xlab = "Genotyping error rate", xlim = c(0,0.02), ylab = expression(paste(log[10], " likelihood"))) 
#MLE = 0.0025 -> write it 
cross <- calc.errorlod(cross, error.prob = MLE)
print(toperr <- top.errorlod(cross, cutoff = 6))
plotGeno(cross, chr = 7, ind = toperrD$id[toperrD$chr==7], cutoff = 6, include.xo = FALSE)
# remove them

for(i in 1:nrow(toperr)) {
  chr <- toperr$chr[i]
  id <- toperr$id[i]
  mar <- toperr$marker[i]
  cross$geno[[chr]]$data[cross$pheno$id==id, mar] <- NA
}

crossProperMapLIKdropProperReorderedLIKMapClean <- cross.clean # backup

# revisit segregation distortion
gt <- geno.table(cross, scanone.output = TRUE)
par(mfrow = c(2,1))
plot(gt, ylab = expression(paste(-log[10], " P-value"))) #very high at the end of lg7
plot(gt, lodcolumn = 3, ylab = "Genotype frequency")
abline(h = c(0.25, 0.5), lty = 2, col = "gray")

svg("genetic_map_final.svg", width = 20, height = 30)
plotMap(cross, show.marker.names=TRUE) #still a hole in lg2 corresponding to re-arranged region
dev.off()

save.image("qtl.RData")

#### QTL analysis ####
summary(cross)
cross$pheno$NA. <- NULL #remove strain names taken as phenotype
par(mar=c(2,2,2,2))
plot(cross)

# single-QTL analysis
cross <- calc.genoprob(cross, step = 1)
cross <- sim.geno(cross, step = 1, n.draws = 100, error.prob = 0.0025)

# all phenotypes with 50 permutations (2 minutes)
out.em <- out.em.perm <- list()
pdf("cross_scanone_em_allPheno_50perm_chro.pdf")
for (i in 1:length(cross$pheno)){
  out.em[[i]] <- scanone(cross, pheno.col = i, method = "em", model = "normal")
  out.em.perm[[i]] <- scanone(cross, pheno.col = i, method = "em", model = "normal", n.perm = 50)
  threshold <- summary(out.em.perm[[i]], alpha=c(0.05))
  par(oma=c(0,0,2,0), mar=c(4,4,1,1), mfrow = c(5,1))
  plot(out.em[[i]][out.em[[i]]$chr==c(1,2),], col = "black", ylim = c(0,6))
  abline(h = threshold[1], col = "red")
  mtext(text = names(cross$pheno)[i], side = 3, line = 1, font = 2)
  plot(out.em[[i]][out.em[[i]]$chr==c(3,4),], col = "black", ylim = c(0,6))
  abline(h = threshold[1], col = "red")
  plot(out.em[[i]][out.em[[i]]$chr==c(5,6),], col = "black", ylim = c(0,6))
  abline(h = threshold[1], col = "red")
  plot(out.em[[i]][out.em[[i]]$chr==7,], col = "black", ylim = c(0,6))
  abline(h = threshold[1], col = "red")
  plot(out.em[[i]][out.em[[i]]$chr==8,], col = "black", ylim = c(0,6))
  abline(h = threshold[1], col = "red")
}
dev.off()

#test
out.em.test <- scanone(cross, pheno.col = 1, method = "em", model = "normal")
out.em.perm.test <- scanone(cross, pheno.col = 1, method = "em", model = "normal", n.perm = 50)
threshold.test <- summary(out.em.perm.test, alpha=c(0.05))
par(oma=c(0,0,2,0), mar=c(4,4,1,1), mfrow = c(5,1))
plot(out.em[[i]][out.em[[i]]$chr==c(1,2),], col = "black", ylim = c(0,6))
mtext(text = names(cross$pheno)[i], side = 3, line = 1, font = 2)
abline(h = threshold[1], col = "red")
plot(out.em[[i]][out.em[[i]]$chr==c(3,4),], col = "black", ylim = c(0,6))
abline(h = threshold[1], col = "red")
plot(out.em[[i]][out.em[[i]]$chr==c(5,6),], col = "black", ylim = c(0,6))
abline(h = threshold[1], col = "red")
plot(out.em[[i]][out.em[[i]]$chr==7,], col = "black", ylim = c(0,6))
abline(h = threshold[1], col = "red")
plot(out.em[[i]][out.em[[i]]$chr==8,], col = "black", ylim = c(0,6))
abline(h = threshold[1], col = "red")


# try interval mapping (last out.em and out.em.perm is Roquefortin C)
lodint(out.em[[1]], chr = 5, expandtomarkers = TRUE)
bayesint(out.em, chr = 5, expandtomarkers = TRUE)
merge(lodint(out.em, chr = 5, expandtomarkers = TRUE), bayesint(out.em, chr = 5, expandtomarkers = TRUE))

# try composite interval mapping on eremefortin A
testCIM <- cim(cross, pheno.col = 29, n.marcovar = 2, window = 10, method = "em", imp.method = "imp", error.prob = 0.0025,
               map.function = "kosambi", addcovar = NULL, n.perm = 0)
testCIM.perm <- cim(cross, pheno.col = 29, n.marcovar = 2, window = 10, method = "em", imp.method = "imp", error.prob = 0.0025,
                    map.function = "kosambi", addcovar = NULL, n.perm = 10)
summary(testCIM) #
summary(testCIM.perm, alpha=0.5)

#### old QTL analysis ####
## see package example
#data(hyper)
#newmap <- est.map(hyper, error.prob=0.01)
#plotMap(hyper, newmap)
#hyper <- replace.map(hyper, newmap)
#plotMap(hyper, newmap)

#load data
cross <- read.cross("csv", file = "./data_qtl.csv", sep = ",", na.string = "NA", genotypes = c("RR", "CR", "CC"))
cross$pheno$NA. <- NULL #remove strain names taken as phenotype

#summary & plot
info <- summary(cross) #checking phenotypes and genotypes (display nind, nphe, nchr, totmar and nmar)
#sink("summary.txt");print(info);sink()
#pdf("summary.pdf")
#par(mar=c(4,4,1,1))
#plot(cross2) #missing genotypes, marker disposition on genome and phenotype frequencies
#dev.off()
#some other useless plots
#plotMap(cross, chr = c(1:6), show.marker.names = F) #marker disposition on genome
#plotMissing(cross) #missing genotypes
#plotPheno(cross, pheno.col = 2) #phenotype frequencies
#recomb <-  drop.nullmarkers(recomb) #remove markers that have no genotype

#copy cross for genetic map
cross2 <- cross

#recombination fraction
cross2 <- est.rf(cross2) #get recombination faction between pairs of adjacent markers -> LOD score
plotRF(cross2) #LOD score map

#genetic map
map <- est.map(cross2, error.prob = 0.05, verbose = TRUE) #genetic map
#pdf("map.pdf")
#plotMap(map) #huge distance in chrom 1 and 3 -> hotspots of recombination ! check what's in there
#dev.off()
cross2 <- replace.map(cross2, map) #replace the initial physical map of the genetic estimated one
cross2 <- calc.errorlod(cross2, error.prob = 0.05) #get LOD score errors
#top.errorlod(cross2) #top error-lod markers
#pdf("errors.pdf")
#plotErrorlod(cross2) #plot them
#dev.off()

#have a look at some genotypes, meaning of the legend please?
#plotGeno(cross2, chr=1, ind=c(1:20))
#plotGeno(cross2, chr=2, ind=c(1:20))
#plotGeno(cross2, chr=3, ind=c(1:20))
#plotGeno(cross2, chr=4, ind=c(1:20))
#plotGeno(cross2, chr=5, ind=c(1:20))
#plotGeno(cross2, chr=6, ind=c(1:20))

#proportion of missing genotype for a chromosome
#pdf("missing.pdf")
#plotInfo(cross2, method = "entropy")
#plotInfo(cross2, method = "variance")
#dev.off()

#test all phenotypes in normal one-QTL model with 6 different methods

#preliminary steps
cross2 <- calc.genoprob(cross2, step = 10, error.prob = 0.05) # 10cM is probably not a good estimate?
cross2 <- sim.geno(cross2, step = 10, n.draws = 64, error.prob = 0.05) # 10cM is probably not a good estimate?

##create the three method output dataframes
#first make a first phenotype
lod.em <- scanone(cross2, pheno.col = 1, method = "em", model = "normal")
lod.mr <- scanone(cross2, pheno.col = 1, method = "mr", model = "normal")
#then only keep the two first columns
lod.em <- lod.em[,1:2]
lod.hk <- lod.em[,1:2]
lod.imp<- lod.em[,1:2]
lod.imp<- lod.em[,1:2]
lod.ehk<- lod.em[,1:2]
lod.mr <- lod.mr[,1:2]
lod.mri<- lod.mr[,1:2]
lod.mra<- lod.mr[,1:2]
pval.em <- lod.em[,1:2]
pval.hk <- lod.em[,1:2]
pval.imp <-lod.em[,1:2]
pval.ehk <-lod.em[,1:2]
pval.mr <- lod.mr[,1:2]
pval.mri<- lod.mr[,1:2]
pval.mra<- lod.mr[,1:2]

for (i in 1:31){
  model <- "normal"
  perm <- 1000 #number of permutation, caution, 100 permutations takes a while for the six methods, you may want to put to 2 for tests
  alpha <- 0.05 #significativity threshold
  #scan (without perm)
  out.em <- scanone(cross2, pheno.col = i, method = "em", model = model)
  out.hk <- scanone(cross2, pheno.col = i, method = "hk", model = model)
  out.imp<- scanone(cross2, pheno.col = i, method = "imp",model = model)
  out.ehk<- scanone(cross2, pheno.col = i, method = "ehk",model = model)
  out.mr <- scanone(cross2, pheno.col = i, method = "mr", model = model)
  out.mri<- scanone(cross2, pheno.col = i, method = "mr-imp",model=model)
  out.mra<- scanone(cross2, pheno.col = i, method = "mr-argmax",model=model)
  #scan with permutation
  out.em.perm <- scanone(cross2, pheno.col = i, method = "em", model = model, n.perm = perm)
  out.hk.perm <- scanone(cross2, pheno.col = i, method = "hk", model = model, n.perm = perm)
  out.imp.perm<- scanone(cross2, pheno.col = i, method = "imp",model = model, n.perm = perm)
  out.ehk.perm<- scanone(cross2, pheno.col = i, method = "ehk",model = model, n.perm = perm)
  out.mr.perm <- scanone(cross2, pheno.col = i, method = "mr", model = model, n.perm = perm)
  out.mri.perm<- scanone(cross2, pheno.col = i, method = "mr-imp",model = model, n.perm = perm)
  out.mra.perm<- scanone(cross2, pheno.col = i, method = "mr-argmax",model = model, n.perm = perm)
  #get significance with scan and permutation
  sum.em <- summary(out.em, perms = out.em.perm, alpha = alpha, pvalues = TRUE)
  sum.hk <- summary(out.hk, perms = out.hk.perm, alpha = alpha, pvalues = TRUE)
  sum.imp<- summary(out.imp,perms = out.imp.perm,alpha = alpha, pvalues = TRUE)
  sum.ehk<- summary(out.ehk,perms = out.ehk.perm,alpha = alpha, pvalues = TRUE)
  sum.mr <- summary(out.mr, perms = out.mr.perm, alpha = alpha, pvalues = TRUE)
  sum.mri<- summary(out.mri,perms = out.mri.perm,alpha = alpha, pvalues = TRUE)
  sum.mra<- summary(out.mra,perms = out.mra.perm,alpha = alpha, pvalues = TRUE)
  #rename to pull into data.frame
  colnames(out.em)[3] <- colnames(cross2$pheno)[i]
  colnames(out.hk)[3] <- colnames(cross2$pheno)[i]
  colnames(out.imp)[3]<- colnames(cross2$pheno)[i]
  colnames(out.ehk)[3]<- colnames(cross2$pheno)[i]
  colnames(out.mr)[3] <- colnames(cross2$pheno)[i]
  colnames(out.mri)[3]<- colnames(cross2$pheno)[i]
  colnames(out.mra)[3]<- colnames(cross2$pheno)[i]
  #cbind the lod scores
  lod.em <- cbind(lod.em, out.em)
  lod.hk <- cbind(lod.hk, out.hk)
  lod.imp<- cbind(lod.imp,out.imp)
  lod.ehk<- cbind(lod.ehk,out.ehk)
  lod.mr <- cbind(lod.mr, out.mr)
  lod.mri<- cbind(lod.mri,out.mri)
  lod.mra<- cbind(lod.mra,out.mra)
  #for pval
  if(substr(sum.em,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.em))){ #for every result/row of "sum", fill the pval
      pval.em$pval[row.names(pval.em)==row.names(sum.em)[j]] <- sum.em$pval[j]
    }
  }
  if(substr(sum.em,1,1)[1]=="i"){
    pval.em$pval <- NA #if not, create an NA column
  }
  if(substr(sum.hk,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.hk))){ #for every result/row of "sum", fill the pval
      pval.hk$pval[row.names(pval.hk)==row.names(sum.hk)[j]] <- sum.hk$pval[j]
    }
  }
  if(substr(sum.hk,1,1)[1]=="i"){
    pval.hk$pval <- NA #if not, create an NA column
  }
  if(substr(sum.imp,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.imp))){ #for every result/row of "sum", fill the pval
      pval.imp$pval[row.names(pval.imp)==row.names(sum.imp)[j]] <- sum.imp$pval[j]
    }
  }
  if(substr(sum.imp,1,1)[1]=="i"){
    pval.imp$pval <- NA #if not, create an NA column
  }
  if(substr(sum.ehk,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.ehk))){ #for every result/row of "sum", fill the pval
      pval.ehk$pval[row.names(pval.ehk)==row.names(sum.ehk)[j]] <- sum.ehk$pval[j]
    }
  }
  if(substr(sum.ehk,1,1)[1]=="i"){
    pval.ehk$pval <- NA #if not, create an NA column
  }
  if(substr(sum.mr,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.mr))){ #for every result/row of "sum", fill the pval
      pval.mr$pval[row.names(pval.mr)==row.names(sum.mr)[j]] <- sum.mr$pval[j]
    }
  }
  if(substr(sum.mr,1,1)[1]=="i"){
    pval.mr$pval <- NA #if not, create an NA column
  }
  if(substr(sum.mri,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.mri))){ #for every result/row of "sum", fill the pval
      pval.mri$pval[row.names(pval.mri)==row.names(sum.mri)[j]] <- sum.mri$pval[j]
    }
  }
  if(substr(sum.mri,1,1)[1]=="i"){
    pval.mri$pval <- NA #if not, create an NA column
  }
  if(substr(sum.mra,1,1)[1]!="i"){ #if there are significant markers
    for (j in 1:length(row.names(sum.mra))){ #for every result/row of "sum", fill the pval
      pval.mra$pval[row.names(pval.mra)==row.names(sum.mra)[j]] <- sum.mra$pval[j]
    }
  }
  if(substr(sum.mra,1,1)[1]=="i"){
    pval.mra$pval <- NA #if not, create an NA column
  }
  #in any case and any method, rename the pval column by the phenotype
  colnames(pval.em)[length(pval.em[1,])] <- colnames(cross2$pheno)[i]
  colnames(pval.hk)[length(pval.hk[1,])] <- colnames(cross2$pheno)[i]
  colnames(pval.imp)[length(pval.imp[1,])]<- colnames(cross2$pheno)[i]
  colnames(pval.ehk)[length(pval.ehk[1,])]<- colnames(cross2$pheno)[i]
  colnames(pval.mr)[length(pval.mr[1,])] <- colnames(cross2$pheno)[i]
  colnames(pval.mri)[length(pval.mri[1,])]<- colnames(cross2$pheno)[i]
  colnames(pval.mra)[length(pval.mra[1,])]<- colnames(cross2$pheno)[i]
}
#remove useless temporary objects
rm(out.em, out.em.perm, out.ehk, out.ehk.perm, out.hk, out.hk.perm, out.imp, out.imp.perm, out.mr, out.mr.perm, out.mra, out.mra.perm,
   out.mri, out.mri.perm, sum.em, sum.ehk, sum.hk, sum.imp, sum.mr, sum.mra, sum.mri)
save.image("qtl.RData")

###chose the best method
##get the number of significative markers for each phenotype for each method
res.em <- colSums(!is.na(pval.em[,3:33])) #total 4 (lipoAsym, saturation, a, andrastin)
res.hk <- colSums(!is.na(pval.hk[,3:33])) #total 19 (lipoAsym, lipoXmid, lipoScal, lipoSlopex2, proteoSlope, Hue, saturation, a, blue, H, MPA, andrastinx3, Isofumx3, RoqC)
res.imp <- colSums(!is.na(pval.imp[,3:33])) #total 18 (lipoAsym, lipoXmid, lipoAsym, lipoSlopex3, proteoSlope, Hue, saturation, a, blue, H, andrastinx3, Isofumx3, roqC)
res.ehk <- colSums(!is.na(pval.ehk[,3:33])) #total 1 (growthSlope)
res.mr <- colSums(!is.na(pval.mr[,3:33])) #total 14 (lipoSlopex2, proteoSlope, Hue, saturation, a, H, MPA, Andrastinx3, Isofumx2, RoqC)
res.mri <- colSums(!is.na(pval.mri[,3:33])) #total 17 (lipoAsym, lipoSlopex3, proteoSlope, Hue, Saturation, a, blue, H, Andrastinx4, Isofumx2, RoqC)
res.mra <- colSums(!is.na(pval.mra[,3:33])) #total 15 (lipoSlopex2, proteoSlope, Hue, Saturation, a, Blue, H, MPA, Andrastinx3, Isofumix2, RoqC)
res <- rbind(res.em, res.hk, res.imp, res.ehk, res.mr, res.mri, res.mra)

rownames(pval.em)[which(!is.na(pval.em[,3]))]
##compare them
for (i in 3:33){
  print(paste(colnames(pval.imp)[i], rownames(pval.imp)[which(!is.na(pval.imp[,i]))], sep = " "))
}

##filtering
#taking the common pvalues
pvalAll <- pval.em[,1:2]
for (i in 3:33){
  pvalAll[,i] <- pmin(pval.em[,i], pval.hk[,i], pval.imp[,i], na.rm = T)
}
colnames(pvalAll) <- colnames(pval.em)
pvalAll <- pvalAll[rowSums(is.na(pvalAll)) != ncol(pvalAll)-2,]

##plot
lod.em


