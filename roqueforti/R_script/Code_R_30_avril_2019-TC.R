library(qtl)

#Charlotte Perrault, 30/04/2019, script used to build a genetic map and QTL analysis

#see package example
data(hyper)
newmap <- est.map(hyper, error.prob=0.01)
plotMap(hyper, newmap)
hyper <- replace.map(hyper, newmap)
plotMap(hyper, newmap)

setwd("/home/thibault/Documents/WORK/CROQ/+/Charlotte_PERRAULT/")

cross <- read.cross("csv", file = "./INDEL_pour_R3_growth.csv", sep = "\t", na.strings = "-", genotypes = c("RR","RC","CC")) #data importation, see QTL example format
length(cross$geno) <- 6 #remove scaffold >6
cross$pheno$NA. <- NULL #remove strains names taken as phenotype
#trying to remove INDEL manually (doesnt work)
#cross$geno$`1`$map$INDEL149 <- NULL
#cross$geno$`1`$data <- cross$geno$`1`$data[,-1]

summary(cross) #checking phenotypes and genotypes (display nind, nphe, nchr, totmar and nmar)
plot(cross) #missing genotypes, marker disposition on genome and phenotype frequencies (this last one doesnt work)
plotMap(cross, chr = c(1:6), show.marker.names = F) #marker disposition on genome
plotMissing(cross) #missing genotypes
plotPheno(cross, pheno.col = 1) #phenotype frequencies (doesnt work)
recomb <-  drop.nullmarkers(recomb) #remove markers that have no genotype

map <- est.map(cross, error.prob = 0.05, verbose = TRUE) #genetic map
plotMap(map) #huge distance in chrom 1 and 3 -> hotspots of recombination ! check what's in there
recomb <- est.rf(cross) #get recombination faction between pairs of adjecent markers -> LOD score
plotRF(recomb) #LOD score map
plotMissing(recomb) #missing genotypes

plotMap(cross, map)
hyper <- replace.map(cross, map) #replace the genetic map of the cross -> why ?
plotMap(hyper, map)
hyper <- calc.errorlod(hyper, error.prob = 0.05) #get LOD score errors
top.errorlod(hyper) #top errorlod markers
plotErrorlod(hyper) #plot them -> how ?

#dont run this !
cross2 <- calc.genoprob(cross, step = 10, err = 0.01)
cross2 <- sim.geno(cross2, step=10, n.draws=64, err=0.01)
out2.hk <- scantwo(cross2, method="hk")
out2.em <- scantwo(cross2)
out2.imp <- scantwo(cross2, method="imp")

cross2 <- calc.errorlod(cross2, error.prob=0.05)
top.errorlod(cross2)
plotGeno(cross2, chr=1, ind=c(24:34, 71:81))
plotGeno(cross2, chr=2, ind=c(24:34, 71:81))
plotGeno(cross2, chr=3, ind=c(24:34, 71:81))
plotGeno(cross2, chr=4, ind=c(24:34, 71:81))
plotGeno(cross2, chr=5, ind=c(24:34, 71:81))
plotGeno(cross2, chr=6, ind=c(24:34, 71:81))
