library(FactoMineR)
library(factoextra)
library(missMDA)
library(corrplot)
library(ggpubr)
library(ade4)
library(adegenet)
library(car)
library(xlsx)

##working directory
setwd("/home/thibault/Documents/WORK/CROQ/QTL/stability_duplication")

##loading data
data <- as.data.frame(read.csv("growth_indel.csv", header = TRUE, dec = ",", sep = "\t"))
data <- data[-c(1,2),] #removes 2 first lines
strains <- sapply(as.vector(data[,1]), function(x){strsplit(strsplit(x, "_")[[1]][2], " ")[[1]][1]}) #get proper strain names
data <- data[,-1] #remove strain names column

##transform CC (fake diploid non-Roquefort), RC (fake heterozygous, duplicated), RR (fake diploid Roquefort) and "-" into 0, 1, 2, NA
data[,2:length(data[1,])] <- apply(data[,2:length(data[1,])], 2, function(x){gsub("CC",0,x)})
data[,2:length(data[1,])] <- apply(data[,2:length(data[1,])], 2, function(x){gsub("RC",1,x)})
data[,2:length(data[1,])] <- apply(data[,2:length(data[1,])], 2, function(x){gsub("RR",2,x)})
data <- apply(data, 2, function(x){gsub("-",NA,x)})
data <- as.data.frame(apply(data, 2, function(x){as.numeric(gsub(",",".",x), na.rm=TRUE)})) #get dot decimal for growth and numeric for all
rownames(data) <- strains

##duplications and growth
data$duplications <- apply(data[,2:length(data[1,])], 1, function(x){return(length(x[x==1]))}) #get number of duplications 
data <- data[,c(1,length(data[1,]),3:length(data[1,])-1)] #put it left side in df

##correlation + fisher tests and scatterplot between growth and duplications
shapiro.test(data$growth); shapiro.test(data$duplications) #p-value<2.2e-16 -> not normally distributed
ggqqplot(data, x = "growth", conf.int = TRUE, title= "Growth slope percentage qqplot") #right skewed ! neither log10, sqrt or ^(1/3) can correct it
ggqqplot(data, x = "duplications", title = "Duplicated marker frequencies qqplot") #left and righ skewed + flatbed (same)
cor.test(data$growth, data$duplications, method = "spearman") #p-val = 0.2578 (rho=0.066)
cor.test(data$growth, data$duplications, method = "kendall") #p-val = 0.2327 (tau=0.048)
tbl <- table(data$growth, data$duplications) #frequency matrix
tbl <- tbl[apply(tbl,1,function(x){!all(x==0)}),apply(tbl,2,function(x){!all(x==0)})] #remove rows and columns full of 0 (only col26)
fisher.test(tbl, simulate.p.value = T, conf.int = T, B = 5000) #estimated p-value 1
chisq.test(tbl) #p-val = 0.4211
ggscatter(data, x = "growth", y = "duplications", add = "reg.line", conf.int = T, cor.coef = F, cor.method = "spearman",
          xlab = "Growth slope", ylab = "Duplicated marker number") #cannot calculate cor.coef for unclear reason (probably because no normality?)

##1.pca
data.imp <- imputePCA(data[,1:length(data[1,])], ncp = 40, scale = T)
pca <- PCA(data.imp$completeObs, ncp = 2, scale.unit = T)

dup_classes <- gsub("\\,","/", gsub("\\]","",gsub("\\(","", cut(data$duplications, breaks = hist(data$duplications, plot=FALSE)$breaks))))
dup_classes[is.na(dup_classes)] <- 0
gwth_classes <- gsub("\\,","/", gsub("\\]","",gsub("\\(","", cut(data$growth, breaks = hist(data$growth, plot=FALSE)$breaks))))
gwth_classes[is.na(gwth_classes)] <- 0

#plot
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 80))
#fviz_pca_var(pca, col.var = "black")
#corrplot(pca$var$contrib, is.corr=FALSE, title = "PCA contribution")
fviz_pca_biplot(pca)
#fviz_pca_var(pca, col.var = "cos2", gradient.cols = rainbow(5), repel = TRUE)#, labelsize= 0.8)
fviz_pca_ind(pca, geom.ind = "text", col.ind = as.factor(gwth_classes), palette = rainbow(length(levels(as.factor(gwth_classes)))),
             addEllipses = TRUE, ellipse.type = "confidence", legend.title = "growth", title = "strain names with growth classes", labelsize=2)
fviz_pca_ind(pca, geom.ind = "text", col.ind = as.factor(dup_classes), palette = rainbow(length(levels(as.factor(dup_classes)))),
             addEllipses = TRUE, ellipse.type = "confidence", legend.title = "duplications", title = "strain names with duplications marker number",
             labelsize=2)
selection <- c("6133","6136","I0037","I0217","I0001","I0186","I1564","I1469","I1581","I0069","I1570","I1566","I1555","I1582","I1616","I1571","I1577",
                             "I0146","I0133","I0011","I1502","I0117","I1321","I0192")
selection <- cbind(selection, unlist(sapply(selection, function(x){x<-round(data.imp$completeObs[which(rownames(data)%in% x)])})), 
                   unlist(sapply(selection, function(x){x<-data[which(rownames(data)%in% x),2]})))
colnames(selection)<- c("strain", "growth", "duplications")
#make plot with selected strains ?

##2.dapc
data2 <- data.imp$completeObs
grp <- find.clusters(data2, n.pca= 50, n.clust = 22, center = F, scale = T)
dapc <- dapc(data2, grp = grp$grp, n.pca = 50, n.da = 5, center = F, scale = T, var.contrib = T)
scatter(dapc, scree.pca = T, cstar=T, mstree=T, pch = 50) #does not work without any clear reason
contrib <- loadingplot(dapc$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
biplot(dapc$ind.coord[,1:2], dapc$var.contr)
selection2<-NULL
for (i in 1:22){
  selection2 <- c(selection2, names(grp$grp[grp$grp==i])[7])
}
selection2 <- c("6133", "6136", selection2)
selection2 <- cbind(unlist(sapply(selection2, function(x){x<-data$growth[which(rownames(data)%in% x)]})), 
                   unlist(sapply(selection2, function(x){x<-data[which(rownames(data)%in% x),2]})))
colnames(selection2)<- c("growth", "duplications")
print(selection2) #remplacer le I1565 'NA en croissance) opar un autre du cluster n°9 -> I1497, 0.3634260, 7
#plot
colSelect <- replace(as.logical(match(rownames(data),selection2)), is.na(as.logical(match(rownames(data),selection2))), FALSE)
dataEllipse(x = dapc$ind.coord[,1], y = dapc$ind.coord[,2], groups = as.factor(grp$grp), level = 0.2, center.pch = NULL, xlab = "DACP 1", ylab = "DACP 2",
            col = rainbow(22))
coeff <- 50
arrows(x0 = rep(0,169), y0 = rep(0,169), x1 = dapc$var.contr[,1]*coeff, y1 = dapc$var.contr[,2]*coeff, length = 0.15)
text(x = dapc$var.contr[,1]*coeff+0.5, y = dapc$var.contr[,2]*coeff+0.5, labels = colnames(data),cex = 0.4)
#legend("bottomright", xjust = 1, legend = levels(as.factor(DataSample$Wave)), pch = 19, cex = 0.8, bg = "transparent", col = rainbow(3))
title(main="DAPC segregating 22 groups in CROQ crossing n°3\nfor recombinant selection", sub="2 disciminant components whithin 50")

#3.tatiana & antoine
data3 <- apply(data[,-c(1:2)],2,function(x){as.numeric(x)})
rownames(data3) <- rownames(data)
#dataDUP <- apply(data3,1,function(x){sum(as.numeric(x)==1,na.rm=T)>1})
dataNDUP <- apply(data3,1,function(x){sum(as.numeric(x)!=1,na.rm=T)==167})
#heatmap(data3[dataDUP,], scale = "none", Colv = NA, Rowv = NA, cexRow = 0.4, cexCol = 0.8, labRow = "")
#axis(4, 1:nrow(data3), rownames(data3), las = 2, cex = 0.2, lwd = 0.5)
DUP.imp <- imputePCA(data3[dataDUP,], scale = T)
#NDUP.imp <- imputePCA(data3[dataNDUP,], scale = T)
#grpDUP <- find.clusters(DUP.imp$completeObs, n.pca= 50, n.clust = 20, center = F, scale = T)
grpNDUP <- find.clusters(data3[dataNDUP,], n.pca= 50, n.clust = 5, center = F, scale = T)

plot(hclust(dist(apply(data3==1,2,function(x){as.numeric(if(sum(x,na.rm=T)>19){as.numeric(x)}else{x*0})}))),cex=0.4)
Choice <- c(3,336,240,339,348,338,253,80,112,62,235,86,267,251,386,388,306,259,107)
DataHM <- heatmap(apply(data3==1 & is.na(data3)==F,2,function(x){as.numeric(if(sum(x)>19){as.numeric(x)}else{x*0})})[Choice,],scale="none")
DataHM2 <- heatmap(apply(data3==1 & is.na(data3)==F,2,function(x){as.numeric(if(sum(x)>19){as.numeric(x)}else{x*0})})[Choice,rev(DataHM$colInd)[1:9]],labRow = rownames(data3)[Choice],scale="none")
selection3 <- rownames(data3)[Choice]
selection3 <- c("I0029","I0039","I0048", "I0107", "I0123", selection3) #add 5 strains with no duplications
selection3 <- cbind(unlist(sapply(selection3, function(x){x<-data$growth[which(rownames(data)%in% x)]})), 
                    unlist(sapply(selection3, function(x){x<-data[which(rownames(data)%in% x),2]})))
colnames(selection3) <- c("growth", "duplications")
#write.table(selection3, file = "sequencing_selection.txt", sep = "\t", quote = F)

#get list of duplicated markers for selection
duplicatd <- NULL
for (i in selection4){
  duplicatd[[i]] <- colnames(data)[which(data[which(rownames(data)==i),]==1)]
    if(length(duplicatd[[i]])<21){
    duplicatd[[i]] <- c(duplicatd[[i]],rep(NA,(21-length(duplicatd[[i]]))))
  }
}
duplicatd <- as.data.frame(duplicatd)
occur <- table(as.vector(apply(duplicatd,2,function(x){return(x)}))[!is.na(as.vector(apply(duplicatd,2,function(x){return(x)})))])
write.xlsx(duplicatd, file = "duplicated_markers2.xlsx", sheet = "duplicated")
write.xlsx(occur, file = "duplicated_markers2.xlsx", sheet = "occurence", append = T)

#number of duplications for selection3
apply(duplicatd,2,function(x){return(sum(!is.na(x)))})
selection4<-rownames(data)[data$duplications>=min(sort(data$duplication,decreasing=T)[1:24])]
