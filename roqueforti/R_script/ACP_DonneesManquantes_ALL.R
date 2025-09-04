library(car)
library(pcaMethods)
library(FactoMineR)
library(factoextra)

#chargement des données
Data<-read.csv("C:/Users/tfccaron/Documents/CROQ/MICROBIOLOGIE/Synthese_microbio_croq_2017.csv",header=TRUE,na.strings="NA",dec=",",sep=";")

#sous-tableaux
DataPCA_cheese<-Data[Data$Stage>=9,] #tout sauf laits
DataPCA_j9<-Data[Data$Stage=="9",] #J9
DataPCA_j20<-Data[Data$Stage=="20",] #J20
DataPCA_j90<-Data[Data$Stage=="90",] #J90
DataPCA_j180<-Data[Data$Stage=="180",] #180

DataSample<-DataPCA_
test1<-"LogPPCAAllCenter"
test2<-"LogPPCAAllCenterWave"

#centrage+normalisation, attention pour les laits
md <- prep(DataSample[,-c(1:8)], center=TRUE)#, scale="uv"

#PCAs
resPPCA <- pca(md, method="ppca", center=FALSE, nPcs=4) #probabiliste permettant les NA
resBPCA <- pca(md, method="bpca", center=FALSE, nPcs=4) #bayesien
resSVDI <- pca(md, method="svdImpute", center=FALSE, nPcs=4) #algorithme de reconstruction des NA
resNipals <- pca(md, method="nipals", center=FALSE, nPcs=4) #non-linear iterative partial least squares
resNLPCA <- pca(md, method="nlpca", center=FALSE, nPcs=4, maxSteps=300) #non-linear pca permettant les NA

#r cumulés
pdf(paste(test1,"_r2cum.pdf",sep=""))
plot(1:4,R2cum(resPPCA),xlab="PCA axis",ylab=expression(paste("Cumulative",R^2,sep=" ")),col="magenta",type="l",ylim=c(0,1))
lines(1:4,R2cum(resBPCA),col="red")
lines(1:4,R2cum(resSVDI),col="green")
lines(1:4,R2cum(resNipals),col="black")
lines(1:4,R2cum(resNLPCA),col="blue")
legend("bottomright",legend=c("PPCA","BPCA","SVDI","Nipals","NLPCA"),pch=19,col=c("magenta","red","green","black","blue"))
dev.off()

#groups = le groupement que tu veux associer à tes points en couleur et forme, level= la proportion de point de chaque groupe que tu veux dans ton ellipse)
pdf(paste(test2,"_dataEllipse.pdf",sep=""))
dataEllipse(x=scores(resPPCA)[,1],y=scores(resPPCA)[,2],
            groups=as.factor(DataSample$Wave),
            level=0.2,center.pch=NULL,
            col=c("green","red","blue"))
#Dessiner les flèche correspondant aux variables et leurs influences sur les axes
zoom<-2
arrows(x0=rep(0,8),y0=rep(0,8),x1=zoom*loadings(resPPCA)[,1],y1=zoom*loadings(resPPCA)[,2],length=0.15)
#Afficher les nombs des variable
text(zoom*loadings(resPPCA)[,1],zoom*loadings(resPPCA)[,2],labels=rownames(loadings(resPPCA)),cex=0.8,pos=1)
#idem pour axes 3 et 4
dataEllipse(x=scores(resPPCA)[,3],y=scores(resPPCA)[,4],
            groups=as.factor(DataSample$Wave),
            level=0.2,center.pch=NULL,
            col=c("green","red","blue","purple"))
arrows(x0=rep(0,8),y0=rep(0,8),x1=zoom*loadings(resPPCA)[,1],y1=zoom*loadings(resPPCA)[,4],length=0.15)
text(zoom*loadings(resPPCA)[,1],zoom*loadings(resPPCA)[,2],labels=rownames(loadings(resPPCA)),cex=0.8,pos=1)
dev.off()
#plot avec factomineR (ne fnoctionne pas avec pcaMethods)
#fviz_pca_ind(resPPCA,
#             geom.ind = "point",
#             col.ind = "Cos2",
#             addEllipses = TRUE,
#             ellipse.level=0.95,
#             legend.title = "Population"
#)

#plot basique
#slplot(resPPCA,sl=DataSample$Strain,hotelling=1)


