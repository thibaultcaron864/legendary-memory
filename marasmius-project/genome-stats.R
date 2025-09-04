#Thibault CARON, 2025/05/18, plots for genome features

library("TreeTools")
library("ggtree")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("ggrepel")
library("dplyr")
library("RColorBrewer")
library("cowplot")
library("patchwork")
library("ggh4x")

#compiling the data
setwd("/home/thibault/Documents/WORK/MARASMIUS/iqtree/")
run90<-read.tree("SUPERMATRIX90.phylip.treefile.cf.tree.cf.tree")
setwd("/home/thibault/Documents/WORK/MARASMIUS/busco")
busco<-read.csv("busco.merged.filtered.csv")
setwd("/home/thibault/Documents/WORK/MARASMIUS/annotation")
general<-read.csv("general.csv")
te<-read.csv("earlgrey/coverage_superfamily_matrix.tsv", sep ="\t")

#merge
genomic<-merge.data.frame(general,te,by="strain")
data<-merge.data.frame(busco,genomic,by="strain")

phylo<-c("flave1","schco1","baemy1","monro1","maror1","marsi1","marfi1","marro1","marcr1","marte1","marsp1","marsc1","mycsc1","mycal2", "mycal1","gyman1",
         "gymaq1","gymer1","lened1","rhodo1","colpe1","colco1","gymlu1")

#filter
dataf<-data[,-which(colnames(data)=="mRNA_number")] #remove mRNA
dataf$strain<-as.factor(dataf$strain)
dataf$strain<-factor(dataf$strain, levels=phylo)
dataf<-cbind(dataf[,1:12],dataf[,46:47],dataf[,13:45]) #re-arrange total repeat and TE pb columns

#PLOTS
#phylo
run90$tip.label <- sapply(strsplit(gsub("BUSCO_","",run90$tip.label),"\\."),"[[",1)
run90<-root(run90, outgroup="schco1")
tip.color <- sapply(run90$tip.label,function(x) ifelse(dataf$online[which(dataf$strain==x)],"black","red"))
run90$tip.label<-sapply(run90$tip.label, function(x) dataf$species[which(dataf$strain==x)])
run90$tip.label <- gsub("^([A-Za-z])[a-z]+\\s+([a-z]+\\s*\\d*)", "\\1. \\2", run90$tip.label)
run90$tip.label <- gsub("Rhodocollybia / ","",run90$tip.label)
#run90$tip.label <- sub(" ", "\n", run90$tip.label)
node_labels <- strsplit(run90$node.label, "/")
alrt      <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[1] else NA))
bootstrap <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[2] else NA))
genes     <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 3) x[3] else NA))
sites     <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[4] else NA))

tree<-ggtree(run90, branch.length = "none", ladderize = T)
tree$data$tip.color<-c(unname(tip.color),rep("NA",21))
tree<-tree+
  geom_tiplab(aes(color=tip.color),align=TRUE, fontface="italic")+
  scale_color_identity()+
  xlim(0,16)+
  theme(plot.margin=margin(t=25,r=0, b=38,l=0))
tree

#busco
colnames(dataf)[11]<-"total_gene_number"
databusco<-melt(dataf, id.vars = c("strain","species","online"), measure.vars = c(4:7), variable.name = "gene")
databusco$species<-factor(databusco$species, levels=unique(as.vector(sapply(levels(databusco$strain), function(x) databusco$species[which(databusco$strain==x)]))))

buscoH<-ggplot(databusco, aes(y=species, x=value, fill=gene))+
  geom_col(position="stack")+
  theme_classic()+
  guides(x = guide_axis(angle=45), fill=guide_legend(title="BUSCO gene category"))+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), axis.title =element_blank())+
  scale_fill_brewer(palette = "Set1")+
  facet_wrap(~ "BUSCO completeness (%)")+
  theme(plot.margin=margin(t=5,r=0, b=22,l=0))
buscolegend<-get_legend(buscoH)
buscoH<-buscoH+theme(legend.position="none")
buscoH

#genomic features
colnames(dataf)[c(9,10,12,13,14)]<-c("N50 (pb)","genome_length (pb)","GC.content (%)","total_repeat_pb (%)","total_TE_pb (%)")
datagen <- melt(dataf, id.vars = c("strain","species","online"), measure.vars = c(8:14), variable.name = "feature")
datagen$species<-factor(datagen$species, levels=unique(as.vector(sapply(levels(datagen$strain), function(x) datagen$species[which(datagen$strain==x)]))))

scaffTH<-N50TH<-10000
datagen <- datagen %>% mutate(extreme = ifelse(feature == "scaffold_number" & value > 10000 | feature == "N50" & value < 10000, value, ""), 
                              value = ifelse(extreme==value, NA, value))

plotgen<-ggplot(datagen, aes(y=strain, x=value,fill=feature))+
  geom_col()+
  facet_grid(~feature, scales = "free_x")+
  guides(x = guide_axis(angle=45))+
  geom_text(data=datagen, aes(x=1,label=extreme), hjust=0, color = "black", size =3)+
  theme_classic()+
  theme(axis.text.y=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), plot.title =element_blank())+
  scale_fill_brewer(palette = "Set2")
genlegend<-get_legend(plotgen)
plotgen<-plotgen+theme(legend.position="none")
plotgen

#TE
te2$TE<-paste(te2$subclass,te2$superfamily,sep=":")
datate<-melt(te2[1:33,],id.vars=c("subclass","TE"), measure.vars = c(3:25), variable.name = "strain")
datate$strain<-factor(datate$strain, levels=phylo)
plottesf<-ggplot(datate, aes(y=strain, x=value,fill=TE))+
  geom_col(position="stack")+
  facet_grid(~subclass,  scales = "free_x")+
  guides(x = guide_axis(angle=45), fill=guide_legend(ncol=7, title="TE subclass:superfamily"))+
  scale_x_continuous(labels = scales::label_number())+
  theme_classic()+
  theme(axis.text.y=element_blank(), axis.title=element_blank(), axis.ticks.y=element_blank(), axis.line.y=element_blank(), plot.title =element_blank())+
  theme(plot.margin=margin(t=5,r=0, b=16,l=0))
telegend<-get_legend(plottesf)
plottesf<-plottesf+theme(legend.position="none")
plottesf

#ggarrange(tree,buscoH,plotgen,plottesf,ncol=4, widths = c(2,1,3,4), heights = c(0.8,0.7,1,0.8), common.legend = T, legend = "right")
legends <- plot_grid(buscolegend, telegend, nrow = 1)
plots<-plot_grid(tree,buscoH,plotgen,plottesf,ncol=4, rel_widths = c(1.5,1,5,4))
plotleg<-plot_grid(plots,legends, ncol=1, rel_heights = c(6,1))
plotleg

png("full-annot.png", width = 2200, height = 1200, units = "px", pointsize = 12, bg = "white", res = 110)
plotleg
dev.off()

svg("full-annot.svg", width = 22, height = 11)
plotleg
dev.off()

##SEPARATE DATA
#genplot
general<-read.csv("general.csv")
general$strain<-factor(general$strain, levels=phylo)

general<-general[,-c(6)] #remove mRNA
general[7,2]<-NA
general[13,2]<-NA

genplot<-melt(general,id.vars="strain",measure.vars = c(2:6), variable.name = "feature",na.rm=F)
genplot$strain <- factor(genplot$strain, levels=rev(levels(genplot$strain)))

plotgen<-ggplot(genplot, aes(y=strain, x=value,fill=feature))+
  geom_col()+
  facet_grid(~feature, scales = "free_x")+
  guides(x = guide_axis(angle=45), fill="none")+
  #geom_text(data = df_labels, aes(x=1,label = value), hjust=0,color = "black", size = 3) +
  #geom_text(data = genplot, aes(x=1,label = value), hjust=0,color = "black", size = 3) +
  theme_classic()
plotgen

png("genome-features.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plotgen
dev.off()

svg("genome-feature.svg", width = 13, height = 11)
plotgen
dev.off()

#TEplot
te<-read.csv("TE.csv")
te<-rbind(te,c("genome_length","genome_length",general$genome_length))
te <- te %>% mutate(across(c(3:25), as.numeric))
teplotsf<-melt(te,id.vars=c("subclass","superfamily"), measure.vars = c(3:25), variable.name = "strain")
teplotsf$strain<-factor(teplotsf$strain, levels=phylo)
teplotsf$strain <- factor(teplotsf$strain, levels=rev(levels(teplotsf$strain)))
plottesf<-ggplot(teplotsf, aes(y=strain, x=value,fill=superfamily))+
  geom_col(position="stack")+
  facet_grid(~subclass, scales = "free_x")+
  guides(x = guide_axis(angle=45))+
  scale_x_continuous(labels = scales::label_number())+
  #geom_text(data = teplot, aes(label = value), position=position_stack(vjust=0.5),color = "black", size = 3) +
  theme_classic()
plottesf

png("genome-te-sf.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plottesf
dev.off()

svg("genome-te-sf.svg", width = 13, height = 11)
plottesf
dev.off()

teplotsc<-melt(te,id.vars="subclass", measure.vars = c(3:25), variable.name = "strain")
teplotsc$subclass<-as.factor(teplotsc$subclass)
teplotsc<-aggregate(value ~ subclass + strain, data = teplotsc, sum)
teplotsc$strain<-factor(teplotsc$strain, levels=c("schco1","flave1","baemy1","monro1","maror1","marfi1","marsi1","marro1","marcr1",
                                                  "marte1","marsp1","mycsc1","mycal1","mycal2","gyman1","gymer1","gymaq1","lened1",
                                                  "rhodo1","colpe1","gymlu1","colco1","marsc1"))
teplotsc$strain <- factor(teplotsc$strain, levels=rev(levels(teplotsc$strain)))
plottesc<-ggplot(teplotsc, aes(y=strain, x=value,fill=subclass))+
  geom_col()+
  facet_grid(~subclass, scales = "free_x")+
  guides(x = guide_axis(angle=45))+
  #geom_text(data = teplot, aes(label = value), position=position_stack(vjust=0.5),color = "black", size = 3) +
  theme_classic()
plottesc

png("genome-te-sc.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plottesc
dev.off()

svg("genome-te-sc.svg", width = 13, height = 11)
plottesc
dev.off()

#te and genome length
telength<-as.data.frame(t(te[c(34:36),c(2:25)]))
colnames(telength)<-telength[1,]
telength<-telength[-c(1),]
telength$strain<-rownames(telength)
rownames(telength)<-NULL
telength$total_pb_repeat <- as.numeric(telength$total_pb_repeat)
telength$total_pb_TE <- as.numeric(telength$total_pb_TE)
telength$genome_length <- as.numeric(telength$genome_length)
ggplot(telength, aes(x=genome_length,y=total_pb_repeat))+
  geom_point(aes(color=strain))+
  theme_classic()+
  geom_smooth(method = "lm", se = FALSE, color = "black")
ggplot(telength, aes(x=genome_length,y=total_pb_TE))+
  geom_point(aes(color=strain))+
  theme_classic()+
  geom_smooth(method = "lm", se = FALSE, color = "black")

#cazyme plot
cazyme<-read.csv("cazyme.csv")
cazyme<-cazyme[-c(5),]
cazyme<-rbind(cazyme,c("genome_length",general$genome_length))
cazyme<-rbind(cazyme,c("gene_number",general$gene_number))
cazyme <- cazyme %>% mutate(across(c(2:24), as.numeric))
cazplot<-melt(cazyme,id.vars="enzyme",measure.vars=c(2:24),variable.name="strain")
cazplot$strain<-factor(cazplot$strain, levels=c("schco1","flave1","baemy1","monro1","maror1","marfi1","marsi1","marro1","marcr1",
                                                  "marte1","marsp1","mycsc1","mycal1","mycal2","gyman1","gymer1","gymaq1","lened1",
                                                  "rhodo1","colpe1","gymlu1","colco1","marsc1"))
cazplot$strain <- factor(cazplot$strain, levels=rev(levels(cazplot$strain)))
plotcaz<-ggplot(cazplot, aes(y=strain, x=value,fill=enzyme))+
  geom_col()+
  facet_grid(~enzyme, scales = "free_x")+
  guides(x = guide_axis(angle=45))+
  scale_fill_discrete(name = "CAZyme", labels = c("Auxilary activities","Carbohydrate Esterases","Gene number","Genome size","Glycoside Hydrolases",
                                                  "Polysaccharide Lyases"))+
  theme_classic()
plotcaz

png("cazyme.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plotcaz
dev.off()

svg("cazyme.svg", width = 13, height = 11)
plotcaz
dev.off()

#enzymes (CAZ + proteases)
enzyme<-read.csv("enzymes.csv")
enzyme<-rbind(enzyme,c("genome_length",general$genome_length))
enzyme<-rbind(enzyme,c("gene_number",general$gene_number))
enzyme <- enzyme %>% mutate(across(c(2:24), as.numeric))
enzplot<-melt(enzyme,id.vars="enzyme",measure.vars=c(2:24),variable.name="strain")
enzplot$strain<-factor(enzplot$strain, levels=c("schco1","flave1","baemy1","monro1","maror1","marfi1","marsi1","marro1","marcr1",
                                                "marte1","marsp1","mycsc1","mycal1","mycal2","gyman1","gymer1","gymaq1","lened1",
                                                "rhodo1","colpe1","gymlu1","colco1","marsc1"))
enzplot$strain <- factor(enzplot$strain, levels=rev(levels(enzplot$strain)))
plotenz<-ggplot(enzplot, aes(y=strain, x=value,fill=enzyme))+
  geom_col()+
  facet_grid(~enzyme, scales = "free_x")+
  guides(x = guide_axis(angle=45))+
  scale_fill_discrete(name = "enzyme", labels = c("Auxilary activities","Carbohydrate Esterases","Glycoside Hydrolases","Polysaccharide Lyases",
                                                  "Gene number","Genome size","Aspactic protease","Metalloprotease","Serine protease"))+
  theme_classic()
plotenz

png("enzymes.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plotenz
dev.off()

svg("enzymes.svg", width = 13, height = 11)
plotenz
dev.off()

