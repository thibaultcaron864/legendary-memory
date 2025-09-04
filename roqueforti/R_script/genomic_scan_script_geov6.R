# Library import ----------------------------------------------------------
library(PopGenome)
library(tidyverse)
library(stats)
library(ggplot2)
library(stringr)
library(svglite)
library(reshape2)
library(gdata)
library(xtable)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(ggsci)
library(MASS)
library(scales)
setwd("/home/thibault/ownCloud/Thèse/gcandidum/genomic_scan")
# Definition of each population -------------------
#Strains info import
strains_info<- read.csv(file = "../../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t")
#Population and individual names
cheese1=strains_info[strains_info$population=="Cheese_1","strains"]
cheese2=strains_info[strains_info$population=="Cheese_2","strains"]
cheese3=strains_info[strains_info$population=="Cheese_3","strains"]
GeoB=strains_info[strains_info$population=="GeoB","strains"]
GeoC=strains_info[strains_info$population=="GeoC","strains"]

#WARNING : YOU HAVE TO DO THIS THE FIRST TIME
#Need to split vcf into scaffold-vcf file in scaffoldtest folder
#YOU HAVE TO LAUNCH THOSE THREE COMMAND ONCE. It will create a subfolder with a vcf for each scaffold. 
#By doing get.sum.data, you will be able to see which scaffolds have no biallelic.sites.
#They have to be removed for the analysis to work because whole.data = FALSE in sliding.window.transform function (see manual)
#VCF_split_into_scaffolds(VCF.file = "snps_pass.vcf","scaffolclib_918")
#df_genome<-readData(path = "scaffolclib_918",format = "VCF")
#get.sum.data(df_genome)

# Summary statistics mean over the genome ---------------------------------
#Once vcf file per scaffold is done using VCF_split_into_scaffold function you can read all scaffold vcf file in the right folder
df_genome<-readData(path = "scaffolclib_918",format = "VCF",include.unknown = T)
#set populations 
df_genome <- set.populations(df_genome,list(cheese1,cheese2,cheese3,GeoB,GeoC))

##Compute population statistics------------------
df_genome<-diversity.stats(df_genome,pi=T) 
#Divide by number of sites to have a statistics per site and make a weighted mean over scaffolds to have whole genome
pi=apply(df_genome@nuc.diversity.within/df_genome@n.sites, 2, weighted.mean,df_genome@n.sites)
#Same with neutrality stat
df_genome<-neutrality.stats(df_genome,FAST = T)
theta_wat=apply(df_genome@theta_Watterson/df_genome@n.sites, 2, weighted.mean,df_genome@n.sites)
##Preparing the table of Nei's Pi and Watterson theta
diversity=rbind(pi,theta_wat)
name=colnames(diversity)#Get colnames to modify them
name=str_replace_all(string = name,pattern = "pop 1","Cheese_1")
name=str_replace_all(string = name,pattern = "pop 2","Cheese_2")
name=str_replace_all(string = name,pattern = "pop 3","Cheese_3")
name=str_replace_all(string = name,pattern = "pop 4","GeoB")
name=str_replace_all(string = name,pattern = "pop 5","GeoC")
#Replace colnames by the modified vector
colnames(diversity)=name
diversity=formatC(diversity,digits = 3,format = "e")#Use a format with less digits and power of ten like "1.00e10"
write.table(diversity,"diversity.csv",sep="\t")#Save diversity stat^

df_genome<-F_ST.stats(df_genome,mode = "nucleotide")
fst=apply(t(df_genome@nuc.F_ST.pairwise), 2, weighted.mean,df_genome@n.sites)
dxy=apply(t(df_genome@nuc.diversity.between)/df_genome@n.sites, 2, weighted.mean,df_genome@n.sites)

pairwise=rbind(fst,dxy)
name=colnames(pairwise)
name=str_replace_all(string = name,pattern = "pop1","Cheese_1")
name=str_replace_all(string = name,pattern = "pop2","Cheese_2")
name=str_replace_all(string = name,pattern = "pop3","Cheese_3")
name=str_replace_all(string = name,pattern = "pop4","GeoB")
name=str_replace_all(string = name,pattern = "pop5","GeoC")
colnames(pairwise)=name
pairwise=as.data.frame(pairwise)
pairwise$stat=rownames(pairwise)
data_long <- gather(pairwise, population,value,colnames(pairwise)[1:10], factor_key=TRUE)
data_long=data_long%>%separate(col = population,into = c("pop1","pop2"),sep = "/")

pfst=ggplot(subset(data_long,stat=="fst"),aes(pop1,pop2,fill=value,label=formatC(value,digits = 3 )))+
  geom_tile()+geom_text(color="white")+theme_bw()+ylab("")+xlab("Fst")
pdxy=ggplot(subset(data_long,stat=="dxy"),aes(pop1,pop2,fill=value,label=formatC(value,digits = 2,format ="e" )))+
  geom_tile()+geom_text(color="white")+theme_bw()+ylab("")+xlab("Nucleotive diversity between (Dxy)")

ggsave(filename ="pairwise.png" ,plot = ggpubr::ggarrange(pdxy,pfst),width = 10,height = 2)

  # 
  # df_genome=calc.fixed.shared(df_genome)
  # fixed_sites=colSums(df_genome@n.fixed.sites)
  # shared_sites=colSums(df_genome@n.shared.sites)
  # monomorphic_sites=colSums(df_genome@n.monomorphic.sites)
  # 
  # snp_stat=rbind(fixed_sites,shared_sites,monomorphic_sites)
  # name=colnames(snp_stat)
  # name=str_replace_all(string = name,pattern = "pop1","Cheese_1")
  # name=str_replace_all(string = name,pattern = "pop2","Cheese_2")
  # name=str_replace_all(string = name,pattern = "pop3","Cheese_3")
  # name=str_replace_all(string = name,pattern = "pop4","GeoB")
  # name=str_replace_all(string = name,pattern = "pop5","GeoC")
  # colnames(snp_stat)=name
  # write.table(snp_stat,"snp_stat.csv",sep="\t")

df_genome<-readData(path = "scaffolclib_918",format = "VCF",include.unknown = T)
#set populations 
df_genome <- set.populations(df_genome,list(cheese1,cheese2,cheese3,GeoB,GeoC))

#get.sum.data(df_genome)
#Using sliding window
win<-7500 #Window size.
df_genome <- sliding.window.transform(df_genome,    width = win, 
                                    jump = 5000,
                                    type = 2,
                                    whole.data = FALSE)
#Number of window in the analysis
length(df_genome@region.names)

##### Computing statistical#####
df_genome<-neutrality.stats(df_genome,FAST = T)
df_genome<-F_ST.stats(df_genome,mode = "nucleotide")
df_genome<-diversity.stats(df_genome,pi=T) 
df_genome<-diversity.stats.between(df_genome)


######  Preparing data #####
#https://evolutionarygenetics.github.io/Chapter8.html
# extract nucleotide diversity and correct for window size
nd <- df_genome@nuc.diversity.within/win
#Thta watterson
theta<-df_genome@theta_Watterson/win
# make population name vector
pops<-c("cheese1","cheese2","cheese3","geob","geoc")
# set population names
colnames(nd) <- paste0(pops, "_pi")
colnames(theta) <- paste0(pops, "_wat")

# extract fst values
fst <- t(df_genome@nuc.F_ST.pairwise)
# extract dxy - pairwise absolute nucleotide diversity
dxy <- df_genome@nuc.diversity.between/win
#Get column names to replace it
x <- colnames(fst)
# replace all occurrences of each population with right name 
x <- sub("pop1", pops[1], x)
x <- sub("pop2", pops[2], x)
x <- sub("pop3", pops[3], x)
x <- sub("pop4", pops[4], x)
x <- sub("pop5", pops[5], x)
# replace forward slash
x <- sub("/", "_", x)
#Replace old colnames with new one 
colnames(fst) <- paste0(x, "_fst")
colnames(dxy) <- paste0(x, "_dxy")

#Create a dataframe with all statistics
cam_data<-cbind.data.frame(nd,theta,fst,dxy)

#Each window is named using scaffold and position, so it's easy to retrieve these informations
#Get scaffold name
cam_data$pos.scaffold<-word(df_genome@region.names, 1)
cam_data$pos.scaffold<-as.numeric(str_replace_all(cam_data$pos.scaffold,pattern = "CCBN0100|\\.1",replacement = ""))
#Get scaffold beginning
cam_data$pos.begin<-as.numeric(word(df_genome@region.names, 2))
#Get scaffold end
cam_data$pos.end<-as.numeric(word(df_genome@region.names, 4))
#Get all position in kb
cam_data$pos.begin<-cam_data$pos.begin/1000 #in kb
cam_data$pos.end<-cam_data$pos.end/1000 #in kb
#Calculate window center for plotting
cam_data$pos.mean<- (cam_data$pos.begin+cam_data$pos.end)/2
#Get a number identifier for each window
cam_data$window<-rownames(cam_data)
#Get a statistic summary of this dataframe
print.xtable(xtable(summary(cam_data),digits = -2),"html",file="summary.html")#Get summary of this dataframe
#__________________________________________________________________________________________________________________________
# Get dataframe ready for plot
df_plot<- melt(cam_data,id.vars = c("pos.mean","pos.begin","pos.end","pos.scaffold","window"),value.name = "value")
df_plot$analysis<-str_sub(df_plot$variable,start = -3)
df_plot$analysis<-str_replace(df_plot$analysis,pattern = "_",replacement = "")
df_plot$analysis<-str_replace(df_plot$analysis,pattern = "wat",replacement = "watterson_theta")
df_plot$variable<-str_replace(df_plot$variable,pattern = "_pi|_dxy|_fst|_wat",replacement = "")
df_plot$variable<-str_replace_all(df_plot$variable,pattern = "camemberti_var_",replacement = "cam.")
df_plot$pos.scaffold<-as.numeric(str_sub(df_plot$pos.scaffold,start = -2))
df_plot<-df_plot[order(df_plot$variable),]
df_plot$variable<-as.factor(df_plot$variable)

# Importing statistics output of vcftools for pi
#pi_cam<-read.table("whole_data_pi_cam_cas",header = T,sep = '\t', stringsAsFactors=FALSE)
#pi_cam$Population<-str_replace_all(string = pi_cam$Population,pattern = "pi_|.txt",replacement = "")
#pi_cam$Population<-str_replace_all(string = pi_cam$Population,pattern = "allcamemberti",replacement = "camemberti.sl")
#pi_cam$Population<-str_replace_all(string = pi_cam$Population,pattern = "camemberti",replacement = "cam.camemberti")
# pi_cam$Population<-str_replace_all(string = pi_cam$Population,pattern = "caseifulvum",replacement = "cam.caseifulvum")
# colnames(pi_cam)<- c("variable","pos.scaffold","pos.begin","pos.end", "N_VARIANTS","PI")
# pi_cam$pos.begin<-pi_cam$pos.begin/1000
# pi_cam$pos.end<-pi_cam$pos.end/1000
# pi_cam$pos.mean<- (pi_cam$pos.begin+pi_cam$pos.end)/2
# pi_cam$pos.scaffold<-as.numeric(str_sub(pi_cam$pos.scaffold,start = -2))
# pi_cam$variable<-as.factor(pi_cam$variable)

# Importing statistics output of vcftools for fst
# fst_cam<-read.table("whole_data_cam_cas",header = T,sep = '\t', stringsAsFactors=FALSE)
# fst_cam$PAIRWISE<-str_replace_all(string = fst_cam$PAIRWISE,pattern = "weir_fst_|.txt",replacement = "")
# fst_cam$PAIRWISE<-str_replace_all(string = fst_cam$PAIRWISE,pattern = "allcamemberti",replacement = "camemberti.sl")
# colnames(fst_cam)<- c("variable","pos.scaffold","pos.begin","pos.end", "N_VARIANTS","WEIGHTED_FST","MEAN_FST")
# fst_cam$pos.begin<-fst_cam$pos.begin/1000
# fst_cam$pos.end<-fst_cam$pos.end/1000
# fst_cam$pos.mean<- (fst_cam$pos.begin+fst_cam$pos.end)/2
# fst_cam$pos.scaffold<-as.numeric(str_sub(fst_cam$pos.scaffold,start = -2))
# fst_cam<-fst_cam[order(fst_cam$variable),]
# fst_cam$variable<-as.factor(fst_cam$variable)
# fst_cam$variable<-fct_inorder(fst_cam$variable)
# fst_var<- c("biforme_fuscoglaucum"="biforme_fuscoglaucum","biforme_camemberti"="cam.camemberti_biforme","camemberti_fuscoglaucum"="cam.camemberti_fuscoglaucum","biforme_caseifulvum"="cam.caseifulvum_biforme","camemberti_caseifulvum"="cam.caseifulvum_cam.camemberti","caseifulvum_fuscoglaucum"="cam.caseifulvum_fuscoglaucum")
# levels(fst_cam$variable)<-fst_var[levels(fst_cam$variable)]
# fst_cam$variable<-factor(fst_cam$variable,levels=sort(levels(fst_cam$variable)))


#MyColorcluster <- c("biforme"="#CCFFFF","cam.camemberti"="#00CC99","cam.caseifulvum"="#006666","fuscoglaucum"="#6666FF")
#MyColorcluster<-MyColorcluster[levels(df_plot_scaffold$pi$variable)]
#colpal<-brewer.pal(length(unique(df_plot_scaffold$dxy$variable)), "Dark2")
#colpal<-c("#e7298a","#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00")
colpal=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

outlier_pi<-df_plot[0,]
outlier_dxy<-df_plot[0,]

######## Overall plot by scaffold ####
for (i in unique(df_plot$pos.scaffold)){#setdiff(1:26,18)
# pi_cam_scaffold<-pi_cam[pi_cam$pos.scaffold==i,]  
# fst_cam_scaffold<-fst_cam[fst_cam$pos.scaffold==i,]
df_plot_scaffold<-df_plot[df_plot$pos.scaffold==i,]
df_plot_scaffold<-split(df_plot_scaffold, df_plot_scaffold$analysis,drop = T)
df_plot_scaffold$fst[df_plot_scaffold$fst$value<0&(!is.na(df_plot_scaffold$fst$value)),"value"]=0
df_plot_scaffold$pi<-droplevels.data.frame(df_plot_scaffold$pi)

df_pi_outlier<-df_plot_scaffold$pi %>% group_by(variable) %>% filter(quantile(value, 0.05,na.rm = T)>value)
df_pi_outlier=df_pi_outlier %>% group_by(pos.scaffold,window)%>%
            mutate(tokeep=any(variable%in%c("cheese1", "cheese2", "cheese3"))&!any(variable%in%"geoc")  )%>%
  filter(tokeep==T)
df_pi_outlier$tokeep=NULL

df_dxy_outlier<-df_plot_scaffold$dxy %>% group_by(variable) %>% filter(quantile(value, 0.99,na.rm = T)<value)



####pi
p1<-ggplot(data =df_plot_scaffold$pi,aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1.5,colour="black",aes(group=variable)) +
  geom_line(size=1)+
  scale_colour_manual(values=colpal)+
  geom_point(data = df_pi_outlier ,mapping = aes(y = value,x=pos.mean,color=variable),size=2,color="black")+
  theme_bw()+scale_x_continuous(expand = c(0,0))+
  ggtitle(paste0("Scaffold N°", i))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),na.value=0,
                                                  labels = trans_format("log10", math_format(10^.x))) +annotation_logticks(sides = "l",outside = T)  +coord_cartesian(clip = "off")+
  xlab("Position (kb)")+ylab("Nucleotide diversity within per site")+ theme(plot.title = element_text(hjust = 0.5),legend.key = element_rect(fill = "grey80"),legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),legend.title=element_blank())
###THETHA WATTERSON
p2<-ggplot(data =df_plot_scaffold$watterson_theta,aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1.5,colour="black",aes(group=variable)) +
  geom_line(size=1)+
  scale_colour_manual(values = colpal)+
  theme_bw()+scale_x_continuous(expand = c(0,0))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),na.value=0,
                            labels = trans_format("log10", math_format(10^.x))) +annotation_logticks(sides = "l",outside = T)  +coord_cartesian(clip = "off")+
  xlab("Position (kb)")+ylab("Watterson theta per site")+ theme(plot.title = element_text(hjust = 0.5),legend.key = element_rect(fill = "grey80"),legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),legend.title=element_blank())


#Dxy
p3<-ggplot(data =df_plot_scaffold$dxy,aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1)+
  geom_point(data = df_dxy_outlier ,mapping = aes(y = value,x=pos.mean,color=variable),size=2,color="black")+
 scale_color_npg()+#xlim(0,max(df_plot_scaffold$dxy$pos.end)+1)+
theme_bw()+scale_x_continuous(expand = c(0,0))+
  xlab("Position (kb)")+ ylab("Dxy")+ theme(legend.key = element_rect(fill = "grey80"),legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),legend.title=element_blank())

p4<-ggplot(data =df_plot_scaffold$fst,aes(y = value,x=pos.mean,color=variable))+
  geom_line(size=1)+scale_x_continuous(expand = c(0,0))+theme_bw()+
  scale_color_npg()+theme_bw()+#xlim(0,max(df_plot_scaffold$dxy$pos.end)+1)+
  xlab("Position (kb)")+ ylab("Fst")+ theme(legend.key = element_rect(fill = "grey80"),legend.background = element_rect(colour = 'black', fill = 'grey90', linetype='solid'),legend.title=element_blank())


ggsave(plot =egg::ggarrange(p1,p2,p3,p4,ncol =  1),filename = paste0("scaffold",i,".jpg"),width=30, height=15,dpi = 300)

#Save outlier windows for pi and fst
df_pi_outlier=df_pi_outlier%>%group_by(pos.begin,pos.scaffold)%>%mutate(variable=paste0(variable, collapse = ";"),value=max(value))%>%distinct(pos.begin,pos.scaffold,variable, .keep_all = T)
outlier_pi<-rbind.data.frame(outlier_pi,df_pi_outlier)


outlier_dxy<-rbind.data.frame(outlier_dxy,df_dxy_outlier)
}
#__________________________________________________________________________________________________________________________
outlier_genomic_scan=rbind.data.frame(outlier_pi,outlier_dxy)
outlier_genomic_scan$scaffold_name<- "CCBN0100000"
outlier_genomic_scan[outlier_genomic_scan$pos.scaffold<10,]$scaffold_name<-paste(outlier_genomic_scan[outlier_genomic_scan$pos.scaffold<10,]$scaffold_name,"0",sep = "")
outlier_genomic_scan$scaffold_name<-paste(sep = "",outlier_genomic_scan$scaffold_name,outlier_genomic_scan$pos.scaffold )
outlier_genomic_scan$scaffold_name=paste0(outlier_genomic_scan$scaffold_name,".1")
outlier_genomic_scan$begin=outlier_genomic_scan$pos.begin*1000
outlier_genomic_scan$end=outlier_genomic_scan$pos.end*1000
write_csv2(x = outlier_genomic_scan,file =   "outlier_genomic_scan.csv")


window_position<-data.frame(tot=df_genome@region.names)
window_position$scaffold<-word(window_position$tot,1)
window_position$begin<-word(window_position$tot,2)
window_position$end<-word(window_position$tot,4)
window_position$scaffold<-as.numeric(str_sub(window_position$scaffold,start = -2))
window_position$tot<-NULL
window_position$window<-rownames(window_position)
write_csv2(x = window_position,path = "window_position.csv")
