#Thibault CARON, 2025-01-23, editing an iqtree

#library(iqtree)
library(treeio)
library(TreeTools)
library(ggplot2)
library(ggtree)

setwd("/home/thibault/Documents/WORK/MARASMIUS/iqtree/")

run90<-read.tree("SUPERMATRIX90.phylip.treefile.cf.tree.cf.tree")
run90$tip.label <- sapply(strsplit(gsub("BUSCO_","",run90$tip.label),"\\."),"[[",1)
run90<-root(run90, outgroup="schco1", resolve.root = T)
new<-c("marro1"=T,"marsi1"=T,"marfi1"=F,"maror1"=F,"lened1"=F,"rhodo1"=T,"colco1"=F,"gymlu1"=F,"colpe1"=T,"gymaq1"=T,"gymer1"=T,"gyman1"=F,
       "mycal2"=T,"mycal1"=T,"mycsc1"=T,"marsc1"=F,"baemy1"=T,"flave1"=F,"schco1"=F,"monro1"=F,"marte1"=F,"marsp1"=F,"marcr1"=F)
tip.color <- ifelse(run90$tip.label %in% names(new) & new[run90$tip.label], "red", "black")
species<-c("marro1"="Marasmius rotula 1","marsi1"="Marasmius siccus 1","marfi1"="Marasmius fiardii 1","maror1"="Marasmius oreades 1",
           "lened1"="Lentinula edodes 1","rhodo1"="Rhodocollybia sp. 1","colco1"="Gymnopus peronatus 1",
           "gymlu1"="Gymnopus luxurians 1","colpe1"="Collybiopsis peronata 1","gymaq1"="Gymnopus aquosus 1","gymer1"="Gymnopus erythropus 1",
           "gyman1"="Gymnopus androsacues 1","mycal2"="Mycetinis alliaceus 2","mycal1"="Mycetinis alliaceus 1","mycsc1"="Mycetinis scorodonius 1",
           "marsc1"="Marasmiellus scandens 1","baemy1"="Baeorospora myosura 1","flave1"="Flammulina velutipes 1","schco1"="Schizophillum commune 1",
           "monro1"="Moniliophtora roreri 1","marte1"="Marasmius tenuissimus 1","marsp1"="Marasmius sp. 1","marcr1"="Marasmius crinis-equi 1",
           "marto1"="Marasmius torquescens 1", "mycsc2"="Mycetinis scorodonius 2")
run90$tip.label<-species[run90$tip.label]
node_labels <- strsplit(run90$node.label, "/")
alrt      <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[1] else NA))
bootstrap <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[2] else NA))
gcf     <- as.numeric(sapply(node_labels,   function(x) if (length(x) >= 3) x[3] else NA))
scf     <- as.numeric(sapply(node_labels,   function(x) if (length(x) >= 2) x[4] else NA))

plot<-ggtree(run90, layout = "rectangular", branch.length = "branch.length", ladderize = TRUE)
plot$data$tip.color<-c(unname(tip.color),rep("NA",22))
plot$data$node.alrt <- c(rep("NA",23),alrt)
plot$data$node.bootstrap <- c(rep("NA",23),bootstrap)
plot$data$node.genes <- c(rep("NA",23),gcf)
plot$data$node.sites <- c(rep("NA",23),scf)

plot<-plot+
  #geom_treescale(x=0, y=-1, width=0.1, label="subst./site") +
  geom_tiplab(aes(color=tip.color),align=TRUE, fontface="italic")+
  scale_color_identity()+
  theme_tree2()+
  #geom_text2(aes(subset = !isTip, label = paste(node.bootstrap,node.sites,node.genes, sep="/")), hjust = 0, vjust = .5, size = 3) +
  geom_text2(aes(subset = !isTip, label = node.sites),     hjust = 1.1, vjust = -0.2, size = 3) +
  geom_text2(aes(subset = !isTip, label = node.genes),     hjust = 1.1, vjust = 1.2,  size = 3) +
  geom_cladelabel(node=29, label = "Marasmiaceae", color="gold", offset=0.28, align=TRUE)+
  geom_cladelabel(node=31, label = "Omphalotaceae", color="purple", offset=0.28, align=TRUE)+
  xlab("Substitutions per site")+
  xlim(0,1.15)+
  geom_hilight(node=29,fill="gold")+
  geom_hilight(node=31,fill="purple")
plot

################
run90<-read.iqtree("SUPERMATRIX90.phylip.treefile.cf.tree.cf.tree")
run90@phylo$tip.label <- sapply(strsplit(gsub("BUSCO_","",run90@phylo$tip.label),"\\."),"[[",1)
tip.color <- ifelse(run90@phylo$tip.label %in% names(new) & new[run90@phylo$tip.label], "red", "black")
run90@data$tip.color <- ifelse(rownames(run90@data) %in% names(new) & new[rownames(run90@data)], "red", "black")
run90@phylo <- root(run90@phylo, outgroup="schco1")
run90@phylo$tip.label <- species[run90@phylo$tip.label]
node_labels<-strsplit(run90@phylo$node.label, "/")
alrt      <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[1] else NA))
bootstrap <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[2] else NA))
genes     <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 3) x[3] else NA))
sites     <- as.numeric(sapply(node_labels, function(x) if (length(x) >= 2) x[4] else NA))

plot <- ggtree(run90, layout="rectangular", branch.length="branch.length", ladderize=TRUE)

plot <- plot +
  geom_tiplab(aes(color=tip.color), align=TRUE, fontface="italic") +
  scale_color_identity() +
  theme_tree2() +
  geom_text2(aes(subset = !isTip, label = paste(bootstrap, gCF, sCF, sep="/")), hjust = 0, vjust = .5, size = 3) +
  geom_cladelabel(node=29, label="Marasmiaceae", color="gold", offset=0.28, align=TRUE) +
  geom_cladelabel(node=31, label="Omphalotaceae", color="purple", offset=0.28, align=TRUE) +
  xlab("Substitutions per site") +
  xlim(0, 1.15) +
  geom_hilight(node=29, fill="gold") +
  geom_hilight(node=31, fill="purple")
plot

png("tree90.png", width = 1300, height = 1100, units = "px", pointsize = 12, bg = "white", res = 110)
plot
dev.off()

svg("tree90.svg", width = 13, height = 11)
plot
dev.off()
