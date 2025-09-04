setwd("/home/thibault/ownCloud/Thèse/gcandidum/gcandidum_snp")
library(scales)
library(FactoClass)
library("phangorn")
library("ggtree")
library("adegenet")
library("ggplot2")
library("grid")
library("gridExtra")
library("ggrepel")
library("svglite")
library("plyr")
library("ape")
library(scatterplot3d)
library(ggsci)
library(dplyr)
options(ggrepel.max.overlaps = Inf)
Pcam <- read.dna(file = "snp_wo_polyploid.fasta", format = "fasta") #snp.fasta

full_Pcam_trait <- read.csv2(file ="../../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t") # dataframe for information about each strains (species, origin)
full_Pcam_trait=full_Pcam_trait%>%group_by(Clonal.group)%>%mutate(clonal=ifelse(test = n()>1,yes = "Clonal",no = "Unique"))
MyColorcluster2 <- c("cheese/dairy"="blue","wild"="red","food"="purple")
MyColorcluster2<-MyColorcluster2[levels(full_Pcam_trait$Environment_simple)]
MyColorcluster4<- c("MAT_A"="blue","MAT_B"="red","MAT_A/MAT_B"="green")
MyColorcluster4<-MyColorcluster4[levels(full_Pcam_trait$Mating_type)]
palette_geo=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

####### PCA #####
Pcam_genind<-DNAbin2genind(x = Pcam)#Convert the DNAbin object into a genind to perform pca
x.cows <- tab(Pcam_genind, freq=TRUE, NA.method="mean")
pca.cows <- dudi.pca(df = x.cows, center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)

df_out<-pca.cows$l1
colnames(df_out)<-c("PC1","PC2","PC3")
df_out$strains <- rownames(df_out)
df_out<-merge.data.frame(df_out, full_Pcam_trait, by = "strains")#Add origin and species to pca data/strains

######### Plotting PCA ############
#Preparation
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))#Prepare nice and beautiful theme
axis_name<-c("PC1","PC2")
explained_variance<-c(pca.cows$eig[1]/sum(pca.cows$eig),pca.cows$eig[2]/sum(pca.cows$eig))
percentage <- round(explained_variance * 100, 2)#Get percentage of explanation of PCA
percentage <- paste( axis_name, "(", paste( as.character(percentage), "%", ")", sep="") )# Get x and y axis label with percentage
#Real plot
p<-ggplot(df_out,aes(x=PC1,y=PC2))+
   geom_point(aes(fill=population,shape=Environment_simple,color=clonal),size=10,stroke=2,alpha=0.65)+
   theme +
   scale_color_manual(values = c("red","black"))+
   xlab(percentage[1]) +
   ylab(percentage[2])+
   scale_fill_manual(values =palette_geo,na.value="white")+
   scale_shape_manual(values=c(21,22,24,23,25))+
   guides(fill=guide_legend(override.aes=list(shape=22)))+
   geom_point(data = subset(df_out,Polyploid==1),aes(x=PC1,y=PC2,size=as.character(Polyploid)),shape=8)+
   scale_size_manual(values=1,labels="Polyploid")+
   theme(legend.key=element_blank())+labs(fill="Population",shape="Origin",size="")
   ggsave(file="SVG/PCA.svg", plot=p, width=12, height=9.6)

p<-ggplot(df_out,aes(x=PC1,y=PC2, label=strains))+
   geom_point(aes(fill=population,shape=Environment_simple),size=10,alpha=0.65,colour="black")+
   theme +
   geom_text_repel(aes(color=population),show.legend = F,segment.alpha = 0.5,segment.colour = "grey60") +
   xlab(percentage[1]) + ylab(percentage[2])+
   scale_fill_manual(values = palette_geo,na.value="white")+
   scale_color_manual(values = palette_geo,na.value="black")+
   scale_shape_manual(values=c(21,22,24,23,25))+
   guides(fill=guide_legend(override.aes=list(shape=22)))+
   geom_point(data = subset(df_out,Polyploid==1),aes(x=PC1,y=PC2,size=as.character(Polyploid)),shape=8)+
   scale_size_manual(values=1,labels="Polyploid")+
   theme(legend.key=element_blank())+labs(fill="Population",shape="Origin",size="")
ggsave(file="SVG/PCA_annotated.svg", plot=p, width=15, height=15)

# 3D ####
axis_name<-c("PC1","PC2","PC3")
explained_variance<-c(pca.cows$eig[1]/sum(pca.cows$eig),pca.cows$eig[2]/sum(pca.cows$eig),pca.cows$eig[3]/sum(pca.cows$eig))
percentage <- round(explained_variance * 100, 2)#Get percentage of explanation of PCA
percentage <- paste( axis_name, "(", paste( as.character(percentage), "%", ")", sep="") )# Get x and y axis label with percentage

# 3D PCA ------------------------------------------------------------------
shapes<-c(21,22,24,23,25)
shapes <- shapes[as.numeric(as.factor(df_out$Environment_simple))]
color3d<-palette_geo
color3d<- color3d[as.numeric(as.factor(df_out$population))]


png(filename = "SVG/PCA_3D.png", width = 9, height = 9, units = "in", res = 400)
# 2. Empty 3D scatter plot using pch=""
s3d <- scatterplot3d(df_out[,2:4], pch = "", grid=FALSE, box=FALSE,y.margin.add = 2.2,xlab = percentage[1],ylab = percentage[2],zlab = percentage[3])
# 3. Add grids
addgrids3d(df_out[,2:4], grid = c("xy", "xz", "yz"))
# 4. Add points
s3d$points3d(df_out[,2:4], pch = shapes, bg=alpha(color3d,0.4),col="black", cex=3)
legend("bottomright", legend = levels(as.factor(df_out$Environment_simple)),
       pch = c(21,22,24,23), xpd = TRUE)
legend("topright", legend = levels(as.factor(df_out$population)),
       col=palette_geo, pch=15)
dev.off()
