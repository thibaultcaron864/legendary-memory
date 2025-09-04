#Thibault CARON, 2024/11/18, plot histograms from a bunch of busco short summary texts
library(rjson)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(dplyr)
library(viridis)
library(RColorBrewer)
library(ggtext)

setwd("/home/thibault/Documents/WORK/MARASMIUS/busco")

old <-        paste(getwd(), "all-old",    sep = "/")
assemblies <- paste(getwd(), "assemblies", sep = "/")
genomes <-    paste(getwd(), "genomes",    sep = "/")
all_phylo <-  paste(getwd(), "all_phylo",  sep = "/")

jsonretrieve<-function(folder,splitdot, splitunder,format){
  data<-data.frame(strain=character(0), matrix(ncol=(format),nrow=0))
  file_list<-list.files(path=folder, pattern = ".json")
  for (i in 1:length(file_list)){
    #print(file_list[i])
    strain <- strsplit((strsplit(file_list[i], "\\.")[[1]][splitdot]),"_")[[1]][splitunder]
    tmp <- fromJSON(file=paste(folder,file_list[i],sep='/'))
    res <- as.data.frame(t(as.data.frame(unlist(tmp$results))))
    res<-cbind(strain,res)
    if(length(res)<(format+1)){
    avg_identity<-NA
    res<-cbind(res[1:13],avg_identity,res[14:format])
    }
    data<-rbind(data, res)
  }
  rownames(data)<-NULL
  data<-data[,-which(colnames(data)=="one_line_summary")]
  return(data)
}
i<-1 #debug
#rm(data)
dataOld <- jsonretrieve(old,5,1,22)
dataPhylo <- jsonretrieve(all_phylo, 4, 2, 22)
dataAss<-jsonretrieve(assemblies,4,2,21)
dataGen<-jsonretrieve(genomes,4,2,21)
data<-rbind(dataAss, dataGen)

write.csv(data,"busco.merged.csv",row.names=F) #print full table

#definitive set
online<-c(baemy1=F, colco1=T,colpe1=F,flave1=T,gyman1=T,gymaq1=F,gymer1=F,gymlu1=T,lened1=T,marcr1=T,marfi1=T,maror1=T,marro1=F,marsc1=T,marsi1=F,
          marsp1=T,marte1=T,marto1=F,monro1=T,mycal1=F,mycal2=F,mycsc1=F,mycsc2=F,rhodo1=F,schco1=T)
phylo<-c("flave1","schco1","baemy1","monro1","maror1","marfi1","marsi1","marro1","marcr1","marte1","marsp1","marsc1","mycsc1","mycal1","mycal2",
         "gyman1","gymer1","gymaq1","lened1","rhodo1","colpe1","colco1","gymlu1")
species<-c("marro1"="Marasmius rotula 1","marsi1"="Marasmius siccus 1","marfi1"="Marasmius fiardii 1","maror1"="Marasmius oreades 1",
           "lened1"="Lentinula edodes 1","rhodo1"="Rhodocollybia sp. 1","colco1"="Collybiopsis confluens 1",
           "gymlu1"="Gymnopus luxurians 1","colpe1"="Collybiopsis peronata 1","gymaq1"="Gymnopus aquosus 1","gymer1"="Gymnopus erythropus 1",
           "gyman1"="Gymnopus androsacues 1","mycal2"="Mycetinis alliaceus 2","mycal1"="Mycetinis alliaceus 1","mycsc1"="Mycetinis scorodonius 1",
           "marsc1"="Marasmiellus scandens 1","baemy1"="Baeorospora myosura 1","flave1"="Flammulina velutipes 1","schco1"="Schizophillum commune 1",
           "monro1"="Moniliophtora roreri 1","marte1"="Marasmius tenuissimus 1","marsp1"="Marasmius sp. 1","marcr1"="Marasmius crinis-equi 1", 
           "marto1"="Marasmius torquescens 1", "mycsc2"="Mycetinis scorodonius 2")
#dataplotV<-data[,c(1,5,7,9,11)] #select raw figures
dataplotV<-data[,c(1,4,6,8,10)] #take percentages
dataplotV<-cbind(dataplotV$strain,online[dataplotV$strain],species[dataplotV$strain], dataplotV[2:length(dataplotV)])
colnames(dataplotV)<-c("strain","online","species","single-copy","multi-copy","fragmented","missing")
dataplotV$`single-copy`<-as.numeric(dataplotV$`single-copy`)
dataplotV$`multi-copy`<-as.numeric(dataplotV$`multi-copy`)
dataplotV$fragmented<-as.numeric(dataplotV$fragmented)
dataplotV$missing<-as.numeric(dataplotV$missing)
dataplotV$strain<-as.factor(dataplotV$strain)
dataplotV$species<-as.factor(dataplotV$species)
dataplotV$strain<-factor(dataplotV$strain, levels=rev(phylo))
dataplotV$species<-factor(dataplotV$species, levels=species[levels(dataplotV$strain)])

write.csv(dataplotV,"busco.merged.filtered.csv",row.names=F) #print full table

#plot definitive
dataplotV <- melt(dataplotV, id.vars = c("species","online"), measure.vars = c(4:7), variable.name = "gene")
#dataplotV$online<-ifelse(dataplotV$online,"red","black")
dataplotV$label <- ifelse(dataplotV$online, paste0("<span style='color:black;'>", dataplotV$species, "</span>"), 
                          paste0("<span style='color:red;'>", dataplotV$species, "</span>"))

buscoplotV<-ggplot(dataplotV, aes(x=species, y=value, fill=gene))+
  geom_col(position="stack")+
  theme_classic()+ 
  theme(axis.text.x =element_markdown(size=12, angle=90, vjust=.5, hjust=1, face="italic"),
        axis.title.x=element_text(vjust=3, margin=margin(8,0,0,0), size=15), 
        axis.title.y=element_text(size=15),
        plot.title =element_text(vjust=0, margin=margin(5,0,5,0)))+
  scale_x_discrete(labels = setNames(dataplotV$label, dataplotV$species))+
  #geom_label_repel(aes(label=value,color=gene), fill="white",position=position_stack(vjust=0.5), size=3, force=0.01, force_pull=2, 
  #                direction="y", alpha=0.75, box.padding = 0, label.size=0, show.legend = F)+
  #labs(title="BUSCO analysis with agaricales odb10 lineage (3870 genes))
  #scale_x_discrete(labels = function(x) {sub("^(\\S+)\\s+", "\\1\n", x)}) +
  labs(fill="Gene classification",y="gene percentage",x="Species")+
  #scale_fill_viridis_d(option = "E")
  scale_fill_brewer(palette = "Pastel1")
buscoplotV

png("busco-gene.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
buscoplotV
dev.off()

svg("busco-gene.svg", width = 13, height = 11)
buscoplotV
dev.off()

#horizontal next to phylogeny
dataplotH<-dataplotV
dataplotH$species <- factor(dataplotH$species, levels=rev(levels(dataplotH$species)))
buscoplotH<-ggplot(dataplotH, aes(y=species, x=value, fill=gene))+
  geom_col(position="stack")+
  theme_classic()+ 
  theme(axis.text.x = element_text(size=12, angle=70, vjust=.5),
        axis.title.x=element_text(vjust=3, margin=margin(15,0,0,0), size=15), 
        axis.title.y=element_text(size=15),
        plot.title =element_text(vjust=0, margin=margin(5,0,5,0)))+
  #geom_label_repel(aes(label=value,color=gene), fill="white",position=position_stack(vjust=0.5), size=3, force=0.01, force_pull=2, 
  #                direction="y", alpha=0.75, box.padding = 0, label.size=0, show.legend = F)+
  #labs(title="BUSCO analysis with agaricales odb10 lineage (3870 genes)")+
  labs(fill="gene percentage",y="",x="")
buscoplotH

png("busco-geneH.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
buscoplotH
dev.off()

svg("busco-geneH.svg", width = 13, height = 11)
buscoplotH
dev.off()


#plot haplotype comparison
dataplothap<-dataOld[,c(1,4,6,8,10)] #select percentages
dataplothap$step<-sapply(strsplit(dataplothap$strain, "-"),"[[",2)
dataplothap$strain<-sapply(strsplit(dataplothap$strain, "-"),"[[",1)
dataplothap<-cbind(dataplothap$strain,dataplothap$step,online[dataplothap$strain],species[dataplothap$strain], dataplothap[3:length(dataplothap)-1])
colnames(dataplothap)<-c("strain","step","online","species","single-copy","multi-copy","fragmented","missing")
dataplothap$`single-copy`<-as.numeric(dataplothap$`single-copy`)
dataplothap$`multi-copy`<-as.numeric(dataplothap$`multi-copy`)
dataplothap$fragmented<-as.numeric(dataplothap$fragmented)
dataplothap$missing<-as.numeric(dataplothap$missing)
dataplothap$strain<-as.factor(dataplothap$strain)
dataplothap$step[which(dataplothap$step=="hifiasm")]<-"primary"
dataplothap$step<-as.factor(dataplothap$step)
dataplothap$species<-as.factor(dataplothap$species)
dataplothap<-dataplothap[-c(which(dataplothap$strain=="marto1"),which(dataplothap$strain=="mycsc2")),]
dataplothap$strain<-factor(dataplothap$strain, levels=rev(phylo))
dataplothap$species<-factor(dataplothap$species, levels=species[levels(dataplothap$strain)])

dataplothap<-dataplothap[dataplothap$step=="primary" | dataplothap$step=="hap1" | dataplothap$step=="hap2",]
dataplothap <- melt(dataplothap, id.vars = c("strain","step"), measure.vars = c(5:8), variable.name = "gene")

plothap<-ggplot(dataplothap, aes(x=paste(strain, step), y=value, fill=gene))+
  geom_col(position="stack")+
  theme_classic()+ 
  theme(axis.text.x =element_text(size=12, angle=90, vjust=.5, hjust=1, face="italic"),
        axis.title.x=element_text(vjust=3, margin=margin(8,0,0,0), size=15), 
        axis.title.y=element_text(size=15),
        plot.title =element_text(vjust=0, margin=margin(5,0,5,0)))+
  scale_fill_brewer(palette = "Pastel1")+
  labs(y="gene percentage",x="Species")
plothap

png("busco-hap.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plothap
dev.off()

svg("busco-hap.svg", width = 13, height = 11)
plothap
dev.off()

#############################################

#assembly data
dataplotassembly <- melt(data[-c(which(data$strain=="marto1"), which(data$strain=="mycsc2")),c(1,9:14)], id = "strain")
dataplotassembly$value<-gsub("%","",dataplotassembly$value)
dataplotassembly$value<-as.numeric(dataplotassembly$value)
#without online strains
dataplotassembly <- melt(data[which(data$online),c(1,9:14)], id = "strain")
#remove purged ones
dataplotassembly <- melt(data[-c(2,5,9,18,23,25,27,30),c(1,9:14)], id = "strain")

assemblyplot<-ggplot(dataplotassembly, aes(x=strain, y=value, fill=variable))+
  geom_bar(stat="identity")+
  theme_classic()+ 
  theme(axis.text.x = element_text(size=9, angle=90, vjust=.7))+
  facet_wrap(~variable, axes="all", scales = "free_y")+
  labs(title="Assemblies general features")
assemblyplot

png("busco-assembly.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
assemblyplot
dev.off()

svg("busco-assembly.svg", width = 13, height = 11)
assemblyplot
dev.off()
