library("ape") 
library("phangorn")
library("ggtree")
library("adegenet")
library("ggplot2")
library("grid")
library("gridExtra")
library("ggrepel")
library("svglite")
library("plyr")
library("ggnetworx")
library(stringr)
setwd("/home/thibault/ownCloud/Thèse/gcandidum/gcandidum_snp")

options(ggrepel.max.overlaps = Inf) 

###### All strains ######
##### Import SNP from fasta ####
DATA <- read.dna(file = "snp.fasta", format = "fasta")
#Convert to the format type that we need :
DATA_phyDat <- phyDat(DATA, type = "DNA", levels = NULL)
DATA_dnabin <- fasta2DNAbin(file = "snp.fasta")

dist_DATA <- dist.dna(DATA, model = "JC69")
DATA_nj <- njs(dist_DATA)

##### Annotated ggplot tree#################
full_DATA_trait <- read.csv2(file ="../../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t") # dataframe for information about each strains (species, origin)
#full_DATA_trait<-droplevels(full_DATA_trait[full_DATA_trait$Polyploid==0,])
full_DATA_trait$strains=str_replace(full_DATA_trait$strains)
full_DATA_trait$newcluster<-as.factor(full_DATA_trait$newcluster)

#From Paul Tol: https://personal.sron.nl/~pault/
# Tol_muted <- c('#88CCEE', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD','#E41A1C','black','cyan','darkred','darkgoldenrod1')
# Tol_muted= pal_d3("category10")(10)
palette_geo=c("Orange", "Orange Red", "blue", "Maroon 1", "darkorchid4", "gray", "gray30","black")

MyColorcluster2 <- c("cheese"="blue","dairy"="lightblue","wild"="red","food"="purple","industrial Standa"="black")
MyColorcluster2<-MyColorcluster2[levels(full_DATA_trait$Environment_simple)]
MyColorcluster3 <- c("GeoA"="blue","GeoB"="lightblue","GeoC"="red","GeoD"="pink")
MyColorcluster3<-MyColorcluster3[levels(full_DATA_trait$Cluster)]
MyColorcluster4<- c("MAT_A"="blue","MAT_B"="red","MAT_A/MAT_B"="green")
MyColorcluster4<-MyColorcluster4[levels(full_DATA_trait$Mating_type)]

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_text_repel(aes(label=label),segment.alpha = 0.3,segment.colour = "grey60")
ggsave(file="SVG/annotated_name_tree.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Environment_simple),colour="black",pch=21,size=5)+
  scale_fill_manual(values = MyColorcluster2)+#Add origin of strains
  guides(fill=guide_legend(title = "Origin"))
ggsave(file="SVG/annotated_tree_env.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Environment_simple),size=5,colour="black",pch=21)+
  scale_fill_manual(values = MyColorcluster2)+#Add origin of strains
  guides(fill=guide_legend(title = "Origin"))+
  geom_text_repel(aes(label=label),segment.alpha = 0.3,segment.colour = "grey60")
ggsave(file="SVG/annotated_name_tree_env.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Cluster),colour="black",size=5,shape=21)+
  scale_fill_manual(values = MyColorcluster3)+#Add origin of strains
  guides(fill = guide_legend(title = "Population",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_tree_cluster.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Mating_type),colour="black",size=5,shape=21)+
  scale_fill_manual(values = MyColorcluster4)+#Add origin of strains
  guides(fill = guide_legend(title = "Mating type",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_tree_mat.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Mating_type),size=5,colour="black",pch=21)+
  geom_text_repel(aes(label=label,color=Mating_type),show.legend=FALSE,segment.alpha = 0.3,segment.colour = "grey60")+
  scale_fill_manual(values = MyColorcluster4)+
  scale_color_manual(values = MyColorcluster4)+
  guides(fill = guide_legend(title = "Mating type",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_name_tree_mat.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=newcluster,shape=Environment_simple),colour="black",size=5)+
  scale_shape_manual(values = c(21:25),na.value=25)+
  scale_fill_manual(values = palette_geo)+#Add origin of strains
guides(shape=guide_legend(title ="Origin"),fill = guide_legend(title ="Population" ,override.aes=list(shape=22)))
ggsave(file="SVG/annotated_tree_env_new_cluster.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01) +#tree scale added
  geom_tippoint(aes(fill=newcluster),colour="black",size=5,shape=21)+
  scale_fill_manual(values = palette_geo)+#Add origin of strains
  guides(fill = guide_legend(title = "Population",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_tree_new_cluster.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=newcluster),colour="black",size=5,shape=21)+
  scale_fill_manual(values = palette_geo)+#Add origin of strains
  guides(fill = guide_legend(title = "Population",override.aes=list(shape=21)))+
  geom_text_repel(aes(label=label,color=newcluster),show.legend=FALSE,segment.alpha = 0.3,segment.colour = "grey60")+
  scale_color_manual(values = palette_geo,na.value ="grey10")
ggsave(file="SVG/annotated_tree_name_new_cluster.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file