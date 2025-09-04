library(FactoMineR)
library(factoextra)
library(car)
library(corrplot)

##working directory
setwd("/home/thibault/Documents/WORK/CROQ/METABARCODING/FROGS/2017/1er_essai")

##loading data
Data<-read.csv("first_results_metab17.csv",
               header=TRUE,dec=",",sep="\t")

##sub dataframe
cheese <- Data #all
j9 <- Data[Data$Stage=="9",] #J9
j20 <- Data[Data$Stage=="20",] #J20
j90 <- Data[Data$Stage=="90",] #J90
j180 <- Data[Data$Stage=="180",] #180

##variables to adjust
DataSample <- Data #sample to select, see above

##PCA
pca <- PCA(DataSample, ncp = 4, scale.unit = TRUE, quali.sup = c(1:5))

##plots
pdf("ACP.pdf")
fviz_eig(pca, addlabels = TRUE, ylim = c(0, 50))
#fviz_pca_var(pca, col.var = "black")
corrplot(pca$var$contrib, is.corr=FALSE)#, title = "PCA contribution"
fviz_pca_var(pca, col.var = "cos2", gradient.cols = rainbow(3), repel = TRUE)
fviz_pca_ind(pca, geom.ind = "point", col.ind = DataSample$Population, 
             palette = rainbow(4), addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Population", title = "PCA on 9-days cheeses : population effect")
fviz_pca_ind(pca, geom.ind = "point", col.ind = as.factor(DataSample$Fabrication), 
             palette = rainbow(9), addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Fabrication", title = "PCA on 9-days cheeses : fabrication effect")
fviz_pca_ind(pca, geom.ind = "point", col.ind = DataSample$Wave, 
             palette = rainbow(3), addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Wave", title = "PCA on 9-days cheeses : wave effect")
fviz_pca_ind(pca, geom.ind = "point", col.ind = as.factor(DataSample$Stage), 
             palette = rainbow(4), addEllipses = TRUE, ellipse.type = "confidence",
             legend.title = "Stage", title = "PCA on all cheeses : stage effect")
dev.off()
####
dataEllipse(x = pca$ind$coord[,1], y = pca$ind$coord[,2],
            groups = as.factor(DataSample$Wave), #points to project
            level = 0.2, center.pch = NULL, #proportion of points in ellipsis
            xlab = paste("PC1 : ", round(pca$eig[1,2], 2), "%", sep = ""),
            ylab = paste("PC2 : ", round(pca$eig[2,2], 2), "%", sep = ""),
            col = rainbow(3))
arrows(x0 = rep(0,12), y0 = rep(0,12),
       x1 = (pca$var$coord[,1]), y1 = (pca$var$coord[,2]), length = 0.15) #Draw the arrows corresponding to the variables and their influences on the axes
text(x = (pca$var$coord[,1]*1.5),
     y = (pca$var$coord[,2]*1.5),labels = colnames(DataSample[,-c(1:5)]),cex = 0.8) #Display variable names
legend("bottomright", xjust = 1, legend = levels(as.factor(DataSample$Wave)),
       pch = 19, cex = 0.8, bg = "transparent", col = rainbow(3))
title("J9 wave")
