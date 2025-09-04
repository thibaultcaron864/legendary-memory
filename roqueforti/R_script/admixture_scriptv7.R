################# Initialization#################"
#Set working directory
setwd("/home/thibault/ownCloud/Thèse/admixture")
#Uploading of libraries
`%nin%` = Negate(`%in%`)
library(gdata)
library(reshape2)
library(ggplot2)
library(ggsci)
library(stringr)
library(dplyr)
library(zoo)
library(randomcoloR)
library(tidyr)


pop<-read.table(header = F,"strains_order")# This list strains contains strains used in the analysis
colnames(pop)<-c("strains")
pop$strains=str_replace(pop$strains,pattern =  "NT12_ESE00548",replacement = "NT12")
pop$strains=str_replace(pop$strains,pattern =  "VTTC4559_ESE00549",replacement = "VTTC4559")


#### Define your range of K that you used, Kmin and Kmax
kmin<-2
kmax<-7
# Read inferred admixture proportions file
temp2<-read.table("results_qopt/NGSadmix_K2_repeat_1.qopt")
temp3<-read.table("results_qopt/NGSadmix_K3_repeat_1.qopt")
temp4<-read.table("results_qopt/NGSadmix_K4_repeat_1.qopt")
temp5<-read.table("results_qopt/NGSadmix_K5_repeat_1.qopt")
temp6<-read.table("results_qopt/NGSadmix_K6_repeat_1.qopt")
temp7<-read.table("results_qopt/NGSadmix_K7_repeat_1.qopt")
# temp8<-read.table("results_qopt/NGSadmix_K8_repeat_1.qopt")
# temp9<-read.table("results_qopt/NGSadmix_K9_repeat_1.qopt")
# temp10<-read.table("results_qopt/NGSadmix_K10_repeat_1.qopt")
# temp11<-read.table("all_wout_conta_wout_poly_K11_repeat_1.qopt")
# temp12<-read.table("all_wout_conta_wout_poly_K12_repeat_1.qopt")
# temp13<-read.table("all_wout_conta_wout_poly_K13_repeat_1.qopt")
# temp14<-read.table("all_wout_conta_wout_poly_K14_repeat_1.qopt")
# temp15<-read.table("all_wout_conta_wout_poly_K15_repeat_1.qopt")
# temp16<-read.table("all_wout_conta_wout_poly_K16_repeat_21.qopt")
# temp17<-read.table("all_wout_conta_wout_poly_K17_repeat_2.qopt")
# temp18<-read.table("all_wout_conta_wout_poly_K18_repeat_61.qopt")
# temp19<-read.table("all_wout_conta_wout_poly_K19_repeat_10.qopt")
# temp20<-read.table("all_wout_conta_wout_poly_K20_repeat_38.qopt")

# First step : Find the order of strains that order the better the plot ###########""
#Create a big dataframe with all data for each K and population ID
df<-cbind.data.frame(mget(paste("temp",kmin:kmax,sep = "")))

# Create the distance matrix
d <- dist(df, method = "euclidean")
# Hierarchical clustering
hc1 <-hclust(d, method = "ward.D2" )
# Plot the obtained dendrogram
#library(dendextend)
#order.hclust(hc1,d)

plot(hc1, cex = 0.6, hang = -1)

#test=dendextend::sort_dist_mat(d,by_cols = F)
# as.matrix(test)
# 
# library(ggdendro)
# # basic option
# ggdendrogram(hc1, size = 4, theme_dendro = FALSE)
# 
# 
# #reordering levels of strains to follow the hclust tree order
# pop$strains <- factor(pop$strains, levels = pop[hc1$order,"strains"])
# #pop$strains<-pop[order(pop$strains),]
#reordering levels of strains to follow the MLtree_order
ml_order=read.table("order_of_strains.csv")$V1
pop$strains <- factor(pop$strains, levels = ml_order)

# Second step : Identify which colors goes where and in which order
# The difficulty is that for different K you have a different number of colors and a different order because columns are independantly ordered
colonne_name<-paste("K",rep(kmin:kmax, times=kmin:kmax),"_",colnames(df),sep = "")
colnames(df)<-colonne_name
df2<-as.data.frame(t(df))
# Dissimilarity matrix
d <- as.matrix(dist(df2, method = "euclidean"))
save<-c()
for (i in kmin:(kmax-1)) {
n<-i
beginr<-n*(n+1)/2-kmin+2
endr<-beginr+n
beginc<-n*(n-1)/2-kmin+2
endc<-beginc+n-1
test<-d[beginc:endc,beginr:endr]
#numeric_vector<-paste(apply(test, 1, which.min))
save_duplicates<-list()
for (j in 1:n) {
  save_duplicates[[j]]<-apply(test, 1,function(x){which(x==sort(x,partial=j)[j])})
}
save_duplicates<-do.call(rbind,save_duplicates)
save_duplicates[1,]<-replace(save_duplicates[1,], duplicated(save_duplicates[1,]), NA)
save_duplicates<-na.locf(save_duplicates,fromLast=T)
for (j in 2:n) {
  save_duplicates[1:j,]<-t(apply(save_duplicates[1:j,], 1, function(x) replace(x, duplicated(x), NA)))
  save_duplicates<-na.locf(save_duplicates,fromLast=T)
}
numeric_vector<-save_duplicates[1,]

while (any( duplicated(numeric_vector))) {
  numeric_vector[duplicated(numeric_vector)]<-sample(which(1:(n+1)%nin%numeric_vector),1)
}
save<-c(save,paste(numeric_vector,collapse = ",",sep = ""))
}

#From Paul Tol: https://personal.sron.nl/~pault/
Tol_muted <- c("blue", '#117733', "gray60", "#00CC99", '#999933',"gray20",  "Maroon 1", '#AA4499', '#DDDDDD','#E41A1C','black','cyan','darkred','darkgoldenrod1','#44AA99')
palette_geo=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")

#Big palette
#Tol_muted<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#000000' , '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff','#fffac8')

#scales::show_col(Tol_muted)
palette<-Tol_muted[1:kmax]
save_palette<-paste(Tol_muted[1:kmax],sep = "",collapse = ",")
for (j in (kmax-1):kmin) {
  palette<- palette[  as.numeric(str_split(save[j-1],pattern = ",",simplify = T)[1,])  ]
  save_palette<-c(save_palette,paste(palette,collapse = ",",sep = ""))
}
save_palette<-rev(save_palette)

#######New color for this K function
newcolor<-c()
for (i in 1:(kmax-kmin)) {
  temp<-str_split(save_palette[i+1],pattern = ',',simplify = T)[1,]
  templag<-str_split(save_palette[i],pattern = ',',simplify = T)[1,]
  newcolor<-c(newcolor,which(temp%nin%templag))
}
############### plot##################
plot_histo <- function(pop,temp,palette) {
  df<-cbind.data.frame(pop,temp)
  df<- melt(df, id.vars = c("strains"))
  k_number<-length(unique(df$variable))
  mycol<-str_split(save_palette[k_number-kmin+1],pattern = ",",simplify = T)
  df<-df%>%arrange(value)%>%group_by(strains)%>%mutate(ordre= row_number())
  #mycol<-as.numeric(str_split(palette[palette$K==k_number,2],pattern = ",",simplify = T))
  #mycol<-pal_npg(palette =  	"nrc")(10)[mycol]
  if(k_number-kmin!=0){
    p<-ggplot(data=df, aes(x=strains, y=value,fill=variable,group=ordre)) +
      geom_bar(stat="identity")+
      theme_void()+ 
      theme(axis.text.y = element_text( size=12))+
      scale_fill_manual(label ="new color",name=paste("K = ",k_number,sep = ""),values = mycol,breaks = levels(df$variable)[newcolor[k_number-kmin]])+
      theme(legend.title.align = 0.5,
            plot.margin = unit(c(1, 5.5, 1, 5.5), "points"),
            text=element_text(size=21))
  }
  if(k_number==kmin){
    p<-ggplot(data=df, aes(x=strains, y=value,fill=variable,group=ordre)) +
      geom_bar(stat="identity")+
      theme_void()+
      theme(axis.text.y = element_text( size=12),
            plot.title = element_text(hjust = 0.5,size=21))+
      ggtitle("Membership of each strains to K clusters")+
      scale_fill_manual(label =rep("new color",length(levels(df$variable))),name=paste("K = ",k_number,sep = ""),values = mycol)+
      theme(legend.title.align = 0.5,
            plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
            text=element_text(size=21))
  }
  return(p)  
}

plot_histo_alone <- function(pop,temp,palette) {
  df<-cbind.data.frame(pop,temp)
  df<- melt(df, id.vars = c("strains"))
  k_number<-length(unique(df$variable))
  mycol<-str_split(save_palette[k_number-kmin+1],pattern = ",",simplify = T)
  df<-df%>%arrange(value)%>%group_by(strains)%>%mutate(ordre= row_number())
  df$strains=forcats::fct_rev(df$strains)
  p<-ggplot(data=df, aes(y=strains, x=value,fill=variable,group=ordre)) +
      geom_bar(stat="identity")+
      theme(axis.text.y = element_text( size=12),
            plot.title = element_text(hjust = 0.5,size=21))+
      ggtitle("Membership of each strains to K clusters")+
      scale_fill_manual(label=c("Cheese_1","Cheese_2","Cheese_3","GeoB","GeoC"),values = palette_geo)+
      theme(legend.title.align = 0.5,
            plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
            text=element_text(size=21))
  return(p)  
}

p2<-plot_histo(pop,temp2,save_palette)
p3<-plot_histo(pop,temp3,save_palette)
p4<-plot_histo(pop,temp4,save_palette)
p5<-plot_histo(pop,temp5,save_palette)
p6<-plot_histo(pop,temp6,save_palette)
p7<-plot_histo(pop,temp7,save_palette)
 #p8<-plot_histo(pop,temp8,save_palette)
 #p9<-plot_histo(pop,temp9,save_palette)
 #p10<-plot_histo(pop,temp10,save_palette)
# p11<-plot_histo(pop,temp11,save_palette)
# p12<-plot_histo(pop,temp12,save_palette)
# p13<-plot_histo(pop,temp13,save_palette)
# p14<-plot_histo(pop,temp14,save_palette)
# p15<-plot_histo(pop,temp15,save_palette)
# p16<-plot_histo(pop,temp16,save_palette)
# p17<-plot_histo(pop,temp17,save_palette)
# p18<-plot_histo(pop,temp18,save_palette)
# p19<-plot_histo(pop,temp19,save_palette)
# p20<-plot_histo(pop,temp20,save_palette)

############Adding informations###########################
# Get ID and pop info for each individual tab delimited 
strains_info<-read.csv(file = "../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t")
pop2=merge(pop,strains_info)
pop3=pop2[order(pop2$strains),]

palette_geo=c("#00CC99", "blue", "Maroon 1", "gray60", "gray20","darkslategrey")
pop3$population=factor(pop3$population,levels =c("Cheese_1","Cheese_2","Cheese_3","GeoB","GeoC") )
leg_strains<-ggplot(data = pop3)+
  geom_tile(aes(fill=population,x=strains,y=1),color="black")+
  scale_fill_manual(name="Population infered on NJ tree with SNP data",values =palette_geo,na.value = "white")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"),
        text=element_text(size=21))+
  scale_y_continuous( expand = c(0, 0)) +
  scale_x_discrete(labels = pop3$strains)

pal_milk=c("goldenrod1","purple","firebrick","darkgreen","tan4","black")
leg5<-ggplot(data = pop3)+
  geom_tile(aes(fill=milk.type,x=strains,y=1),color="black")+theme_void()+
  scale_fill_manual(name="Milk type",values = pal_milk)+
  theme(plot.margin = unit(c(1, 5.5, 1, 5.5), "points"),
        text=element_text(size=21),
        legend.position = "bottom")+
  scale_y_continuous( expand = c(0, 0))+
  guides(fill = guide_legend(nrow = 1,override.aes = list(size = 1)))

leg5<-ggplot(data = pop3)+
  geom_tile(aes(fill=milk.type,x=strains,y=1),color="black")+theme_void()+
  scale_fill_manual(name="Milk type",values = pal_milk)+
  theme(plot.margin = unit(c(1, 5.5, 1, 5.5), "points"),
        text=element_text(size=21),
        legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1,override.aes = list(size = 1)))

plot_list<-mget(c(paste0("p",kmin:kmax),"leg5","leg_strains"))
#egg::ggarrange(plots=plot_list,ncol = 1,  heights = c(rep(3,kmax-kmin+1),1,1,1)  )

ggsave(plot =egg::ggarrange(plots=plot_list,ncol = 1,  heights = c(rep(3,kmax-kmin+1),1,1)  ),
       filename = "NGSadmix_results_ward.D2_euclidean.png",units="cm", width=60, height=50, dpi = 300)

# multiplot ---------------------------------------------------------------

for (f in 2:10) {

######################## Checking likeliness of different run for same k  ###############"
K<-f
temp = list.files(path = "results_qopt",pattern=paste("all_wout_conta_with_poly_K",K,"_repeat_[0-9]*.qopt",sep = ""),full.names = T)
myfiles = lapply(temp, read.table)
test<-do.call(cbind,myfiles)

#colnames(test)<-c(rep("FIRST",each=K),rep(2:100,each=K))
colnames(test)<-rep(1:100,each=K)#paste("rep_",rep(1:100,each=K),"_V_",rep(1:K,times=K),sep="")
d <- dist(t(test), method = "euclidean")
# Hierarchical clustering
hc1 <- hclust(d, method = "ward.D2" )
# Plot the obtained dendrogram
jpeg(filename = paste("tree_K",K,sep = ""),units="cm", width=70, height=10, res = 300)
plot(hc1, cex = 0.6, hang = -1,cex=0.3)
rect.hclust(hc1,k = K,border = "blue")
rect.hclust(hc1,k = K+1,border = "red")
dev.off()


# Create the distance matrix
d <- dist(test, method = "euclidean")
# Hierarchical clustering
hc1 <- hclust(d, method = "ward.D2" )
# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)
#reordering levels of strains to follow the hclust tree order
pop$strains <- factor(pop$strains, levels = pop[hc1$order,"strains"])

df2<-as.data.frame(t(test))
rownames(df2)<-paste("rep_",rep(1:100,each=K),"_V_",rep(1:K,times=100),sep="")
# Dissimilarity matrix
d <- as.matrix(dist(df2, method = "euclidean"))
save<-c()
for (i in 1:99) {
  n<-i
  beginc<-n*K +1
  endc<-n*K+K
  beginr<-1
  endr<-K
  temp<-d[beginc:endc,beginr:endr]
  #numeric_vector<-paste(apply(temp, 1, which.min))
  
  save_duplicates<-list()
  for (j in 1:K) {
    save_duplicates[[j]]<-apply(temp, 1,function(x){which(x==sort(x,partial=j)[j])})
  }
  save_duplicates<-do.call(rbind,save_duplicates)
  save_duplicates<-t(apply(save_duplicates,1,function(x){replace(x,duplicated(x),NA)}))
  save_duplicates<-na.locf(save_duplicates,fromLast=T)

  numeric_vector<-save_duplicates[1,]
  

  while (any( is.na(numeric_vector))) {
    numeric_vector[is.na(numeric_vector)]<-sample(which(1:K%nin%numeric_vector),1)
  }
  save<-c(save,paste(numeric_vector,collapse = ",",sep = ""))
}

#Big palette
Tol_muted<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#000000' , '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff','#fffac8')
#scales::show_col(Tol_muted)
palette<-Tol_muted[1:K]
save_palette<-paste(Tol_muted[1:K],sep = "",collapse = ",")
for (j in 1:99){
  numeric_vector<-as.numeric(str_split(save[j],pattern = ",",simplify = T)[1,])
  palette_temp<-palette[numeric_vector]
save_palette<-c(save_palette,paste(palette_temp,collapse = ","))
}


save<-c()
for (i in 1:99) {
  n<-i
  beginc<-n*K +1
  endc<-n*K+K
  beginr<-1
  endr<-K
  temp<-d[beginc:endc,beginr:endr]
  #numeric_vector<-paste(apply(temp, 2, which.min))
  
  save_duplicates<-list()
  for (j in 1:K) {
    save_duplicates[[j]]<-apply(temp, 2,function(x){which(x==sort(x,partial=j)[j])})
  }
  save_duplicates<-do.call(rbind,save_duplicates)
  save_duplicates<-t(apply(save_duplicates,1,function(x){replace(x,duplicated(x),NA)}))
  save_duplicates<-na.locf(save_duplicates,fromLast=T)
  
  numeric_vector<-save_duplicates[1,]
  
  while (any( is.na(numeric_vector))) {
    numeric_vector[is.na(numeric_vector)]<-sample(which(1:K%nin%numeric_vector),1)
  }
  
  save<-c(save,paste(numeric_vector,collapse = ",",sep = ""))
}
#####Get repeat order that is more convenient #####
caca<-c(1:K,as.numeric(str_split(paste0(save,collapse = ","),pattern = ",",simplify = T)[1,]) )
caca2<-rep(0:99,each=K)*K
caca3<-caca+caca2
df3<-df2[caca3,]
df3$rep<-rep(1:100,each=K)
df3$var<-rep(1:K,times=100)
rownames(df3)<-paste("rep",df3$rep,"var",df3$var)
list_distance_matrix<-list()
for (i in 1:K) {
  matrix_temp<-dist(df3[df3$var==i,colnames(df3)%nin%c("var","rep")], method = "euclidean")
 matrix_temp<- matrix_temp/max(matrix_temp)
list_distance_matrix[[i]]<-matrix_temp
}

bla<-Reduce('+', list_distance_matrix)

hc2 <- hclust(bla, method = "ward.D2" )
# Plot the obtained dendrogram
plot(hc2, cex = 0.6, hang = -1)
order_for_plot<-hc2$order


new<-as.data.frame(t(test))
new$rep<-rep(1:100,each=K)
new$var<-rep(1:K,times=100)
new<-new%>%gather(individual,value,-c(rep,var))
new$individual<-as.numeric(str_sub(start=2,new$individual))
new$individual<-as.factor(new$individual)
levels(new$individual)<-pop$strains
new$individual <- factor(new$individual, levels = pop[hc1$order,"strains"])
new<-new%>%arrange(value)%>%group_by(individual)%>%mutate(ordre= row_number())
new$var<-as.factor(new$var)

myplot_list<-list()
for (i in 1:100) {
my_col<-str_split(save_palette[i],pattern = ",",simplify = T)

myplot_list[[i]]<-ggplot(data=new[new$rep==i,], aes(x=individual, y=value,fill=var,group=ordre)) +#,group=ordre
  geom_bar(stat="identity")+
  scale_fill_manual(values=my_col)+
  theme_void()+
  theme(legend.title.align = 0.5,
        axis.title.y = element_text(),
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
        text=element_text(size=21))+guides(fill=F)+
  ylab(paste("rep",i))

# assign(paste0("test",i),ggplot(data=new[new$rep==i,], aes(x=individual, y=value,fill=var,group=ordre)) +#,group=ordre
#   geom_bar(stat="identity")+
# scale_fill_manual(values=my_col)+
#   theme_void()+
#   theme(legend.title.align = 0.5,
#         plot.margin = unit(c(5.5, 5.5, 1, 5.5), "points"),
#         text=element_text(size=21))+guides(fill=F))
}

pop$strains <- factor(pop$strains, levels = pop[hc1$order,"strains"])

strains_info<-read.table(file = "../pop_info/Gcandidum summary - Feuille 1.tsv",header = T,sep = "\t")
pop2=merge(pop,strains_info)
pop3=pop2[order(pop2$strains),]

leg_strains<-ggplot(data = pop3)+
  geom_tile(aes(fill=newcluster,x=strains,y=1),color="black")+
  scale_fill_manual(name="Population infered on NJ tree with SNP data",values =pal_d3("category10")(10))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position="bottom",
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points"),
        text=element_text(size=21))+
  scale_y_continuous( expand = c(0, 0)) +
  scale_x_discrete(labels = pop3$strains)
leg_strains


myplot_list_reordered<-myplot_list[order_for_plot]
myplot_list_reordered[[101]]<-leg_strains
ggsave(filename = paste("Plot_rep_","K",K,sep = ""),device = "png",
       plot = egg::ggarrange(plots =myplot_list_reordered ,ncol = 1),
       width=20, height=50, dpi = 300,limitsize = F)
}

