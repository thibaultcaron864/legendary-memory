library(ape)
library(phangorn)
library(seqinr)

###### All strains ######
##### Import SNP from fasta ####
DATA <- read.dna(file = "snp_wo_polyploid.fasta", format = "fasta")
#Convert to the format type that we need :
DATA_phyDat <- phyDat(DATA, type = "DNA", levels = NULL)
# random_sample=sample(1:99, 10)
# subDATA_phydat <- subset(DATA_phyDat, random_sample)

#mt <- modelTest(subDATA_phydat)
# Model df   logLik     AIC          AICw    AICc         AICcw     BIC
# 23   GTR+G 26 -3088712 6177476  1.000000e+00 6177476  1.000000e+00 6177784
# 24 GTR+G+I 27 -3088956 6177966 5.761958e-107 6177966 5.761652e-107 6178285
# 22   GTR+I 26 -3089168 6178389 7.613307e-199 6178389 7.613307e-199 6178696
# 15   HKY+G 22 -3097702 6195449  0.000000e+00 6195449  0.000000e+00 6195709
# 14   HKY+I 22 -3097971 6195987  0.000000e+00 6195987  0.000000e+00 6196247
# 16 HKY+G+I 23 -3098169 6196385  0.000000e+00 6196385  0.000000e+00 6196657
# 19   SYM+G 23 -3100317 6200680  0.000000e+00 6200680  0.000000e+00 6200953
# 18   SYM+I 23 -3100612 6201270  0.000000e+00 6201270  0.000000e+00 6201542
# 20 SYM+G+I 24 -3100712 6201473  0.000000e+00 6201473  0.000000e+00 6201757
# 11   K80+G 19 -3106028 6212094  0.000000e+00 6212094  0.000000e+00 6212319
# 10   K80+I 19 -3106295 6212628  0.000000e+00 6212628  0.000000e+00 6212852
# 12 K80+G+I 20 -3106490 6213021  0.000000e+00 6213021  0.000000e+00 6213257
# 21     GTR 25 -3171571 6343192  0.000000e+00 6343192  0.000000e+00 6343488
# 13     HKY 21 -3181833 6363709  0.000000e+00 6363709  0.000000e+00 6363957
# 17     SYM 22 -3184836 6369717  0.000000e+00 6369717  0.000000e+00 6369977
# 9      K80 18 -3191940 6383915  0.000000e+00 6383915  0.000000e+00 6384128
# 6    F81+I 21 -3234804 6469651  0.000000e+00 6469651  0.000000e+00 6469899
# 7    F81+G 21 -3235357 6470756  0.000000e+00 6470756  0.000000e+00 6471005
# 8  F81+G+I 22 -3236304 6472651  0.000000e+00 6472651  0.000000e+00 6472912
# 2     JC+I 18 -3241503 6483043  0.000000e+00 6483043  0.000000e+00 6483256
# 3     JC+G 18 -3242048 6484132  0.000000e+00 6484132  0.000000e+00 6484345
# 4   JC+G+I 19 -3243008 6486054  0.000000e+00 6486054  0.000000e+00 6486278
# 5      F81 20 -3312079 6624197  0.000000e+00 6624197  0.000000e+00 6624434
# 1       JC 17 -3319477 6638988  0.000000e+00 6638988  0.000000e+00 6639189

#Best model is GTR+G
# choose best model from the table according to AICc
#bestmodel <- mt$Model[which.min(mt$AICc)]

dna_dist <- dist.ml(DATA_phyDat, model="JC69")
NJ_tree<- NJ(dna_dist)

fitStart <- pml(tree = NJ_tree,data = DATA_phyDat,model = "GTR")
fit <- optim.pml(fitStart, rearrangement = "stochastic",
                  optInv=TRUE, model="GTR")
bs <- bootstrap.pml(fit, bs=100, optNni=TRUE, multicore=TRUE,mc.cores = 50)


plotBS(midpoint(fit$tree), bs, p = 50, type="p")

write.tree(bs, file="bootstrap_mltree.tre")
write.nexus(bs,file = "bootstrap_mltree.nex")
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





####################Cheese only #######################

##### Import SNP from fasta ####
DATA <- read.dna(file = "snp_cheese_wo_polyploid.fasta", format = "fasta")
#Convert to the format type that we need :
DATA_phyDat <- phyDat(DATA, type = "DNA", levels = NULL)
DATA_dnabin <- fasta2DNAbin(file = "snp_cheese.fasta")

dist_DATA <- dist.dna(DATA, model = "JC69")
DATA_nj <- njs(dist_DATA)

##### Annotated ggplot tree#################
MyColorcluster2 <- c("cheese"="blue","dairy"="lightblue","wild"="red","food"="purple","industrial Standa"="black")
MyColorcluster2<-MyColorcluster2[levels(full_DATA_trait$Environment_simple)]
MyColorcluster3 <- c("GeoA"="blue","GeoB"="lightblue","GeoC"="red","GeoD"="pink")
MyColorcluster3<-MyColorcluster3[levels(full_DATA_trait$Cluster)]
MyColorcluster4<- c("MAT_A"="blue","MAT_B"="red","MAT_A/MAT_B"="green")
MyColorcluster4<-MyColorcluster4[levels(full_DATA_trait$Mating_type)]

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_text_repel(aes(label=label),segment.alpha = 0.3,segment.colour = "grey60")
ggsave(file="SVG/annotated_name_tree_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Environment_simple),colour="black",pch=21,size=5)+
  scale_fill_manual(values = MyColorcluster2)+#Add origin of strains
  guides(fill=guide_legend(title = "Origin"))
ggsave(file="SVG/annotated_tree_env_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Environment_simple),size=5,colour="black",pch=21)+
  scale_fill_manual(values = MyColorcluster2)+#Add origin of strains
  guides(fill=guide_legend(title = "Origin"))+
  geom_text_repel(aes(label=label),segment.alpha = 0.3,segment.colour = "grey60")
ggsave(file="SVG/annotated_name_tree_env_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Cluster),colour="black",size=5,shape=21)+
  scale_fill_manual(values = MyColorcluster3)+#Add origin of strains
  guides(fill = guide_legend(title = "Population",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_tree_cluster_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Mating_type),colour="black",size=5,shape=21)+
  scale_fill_manual(values = MyColorcluster4)+#Add origin of strains
  guides(fill = guide_legend(title = "Mating type",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_tree_mat_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file

annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=Mating_type),size=5,colour="black",pch=21)+
  geom_text_repel(aes(label=label,color=Mating_type),show.legend=FALSE,segment.alpha = 0.3,segment.colour = "grey60")+
  scale_fill_manual(values = MyColorcluster4)+
  scale_color_manual(values = MyColorcluster4)+
  guides(fill = guide_legend(title = "Mating type",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_name_tree_mat_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=newcluster,shape=Environment_simple),colour="black",size=5)+
  scale_shape_manual(values = c(21:25),na.value=25)+
  scale_fill_manual(values = palette_geo)+#Add origin of strains
  guides(shape=guide_legend(title ="Origin"),fill = guide_legend(title ="Population" ,override.aes=list(shape=22)))
ggsave(file="SVG/annotated_tree_env_new_cluster_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01) +#tree scale added
  geom_tippoint(aes(fill=newcluster),colour="black",size=5,shape=21)+
  scale_fill_manual(values = palette_geo)+#Add origin of strains
  guides(fill = guide_legend(title = "Population",override.aes=list(shape=21)))
ggsave(file="SVG/annotated_tree_new_cluster_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


annotated_tree<-ggtree(DATA_nj,layout = "daylight")%<+% full_DATA_trait+ #Adding dataframe to the ggtree object
  geom_treescale(fontsize=6, linesize=2, offset=0.01)+ #tree scale added
  geom_tippoint(aes(fill=newcluster),colour="black",size=5,shape=21)+
  scale_fill_manual(values = palette_geo)+#Add origin of strains
  guides(fill = guide_legend(title = "Population",override.aes=list(shape=21)))+
  geom_text_repel(aes(label=label,color=newcluster),show.legend=FALSE,segment.alpha = 0.3,segment.colour = "grey60")+
  scale_color_manual(values = palette_geo,na.value ="grey10")
ggsave(file="SVG/annotated_tree_name_new_cluster_cheese.svg", plot=annotated_tree, width=15, height=10)# Save the plot into a svg file


