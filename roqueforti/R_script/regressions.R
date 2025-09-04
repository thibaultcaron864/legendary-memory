#thibault CARON, 21/09/23, figure REG

library(ggplot)
library(ggpubr)

#data
setwd("/home/thibault/Documents/WORK/CROQ/QTL/phenotypes/lipoProteo")
regLipo <- read.csv("./2_regression/lipoReg_all.csv")
regProteo <- read.csv("./2_regression/proteoReg_all.csv")
lipoProteo <- read.csv("./1_raw/lipolysis_proteolysis_all.csv")
setwd("/home/thibault/Documents/WORK/CROQ/QTL/phenotypes/growth")
regGrowth <- read.csv("./reglin/reglin.csv")

#medians for lipo et proteo
median(regLipo$rsquared, na.rm = T) #take the forst one I0033
median(regProteo$rsquared, na.rm = T) #take the only one I0919
median(regGrowth$rsquared, na.rm = T) #take I1497 -> box 4879

growData <- read.csv("./csv_corrected/19-04-2018/boite_4879_colonies_detaillees.csv", header=T, sep=";")
proteoData <- lipoProteo[which(lipoProteo$strain=="I0919"),c(3:6)]
lipoData <- lipoProteo[which(lipoProteo$strain=="I0919"),c(7:10)]

rm(lipoProteo, regLipo, regProteo, regGrowth)

setwd("/home/thibault/Documents/WORK/CROQ/QTL/phenotypes")
growData$Time <- as.numeric(gsub(",",".",gsub(" ","",gsub("\\+","",gsub("h","",growData$Temps)))))

#first try on growth
# par(mar=c(4,4,1,1))
# plot(growData$Diametre~growData$Time, xlab="",ylab="")
# title(xlab="time (hours)",cex.lab=1.5,line=2.5)
# #title(ylab="colony surface (pixels)")
# title(ylab="colony diameter (pixels)",cex.lab=1.5,line=2.5)
# abline(lm(growData$Diametre~growData$Time), lwd=2, col="red")
# a<-round(summary(lm(growData$Diametre~growData$Time))$coefficients[2,1],digits=2)
# eq<-paste("y = ",a,"x ",b,sep="")
# r2<-round(summary(lm(growData$Diametre~growData$Time))$r.squared, digits=4)
# mtext(text=eq,side=3,line=-1.5,at=49, cex=1.5,col="red")
# mtext(text=expression("r"^2*" = 0.985"),side=3,line=-3.5,at=45.3,cex=1.5,col="red")

#growth
lm_eqn <- function(ord, abs, df){
  m <- lm(ord ~ abs, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}
growPlot <- ggplot(data = growData, aes(x = Time, y = Diametre)) + 
  geom_line(linewidth = 2) +
  labs(x = "Time (hours)", y = "Diameter (pixels)") +
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  theme_classic() +
  theme(plot.title = element_text(size=20), axis.title=element_text(size=15)) +
  geom_text(x=75, y=max(growData$Diametre), label = lm_eqn(growData$Diametre, growData$Time, growData), parse = TRUE, cex= 5)

#lipo
lipoData <- as.data.frame(cbind(t(lipoData),c(7,14,21,28)))
colnames(lipoData) <- c("Lipolysis", "Time")
lipoPlot <- ggplot(data = lipoData, aes(x = Time, y = Lipolysis)) + 
  geom_line(linewidth = 2) +
  labs(x = "Time (days)", y = "Lipolysis front (mm)") +
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  theme_classic() +
  theme(plot.title = element_text(size=20), axis.title=element_text(size=15)) +
  geom_text(x=16, y=max(lipoData$Lipolysis), label = lm_eqn(lipoData$Lipolysis, lipoData$Time, lipoData), parse = TRUE, cex= 5)

#proteo
proteoData <- as.data.frame(cbind(t(proteoData),c(7,14,21,28)))
colnames(proteoData) <- c("Proteolysis", "Time")
proteoPlot <- ggplot(data = proteoData, aes(x = Time, y = Proteolysis)) + 
  geom_line(linewidth = 2) +
  labs(x = "Time (days)", y = "Proteolysis front (mm)") +
  geom_smooth(method = "lm", se=FALSE, color="blue", formula = y ~ x) +
  theme_classic() +
  theme(plot.title = element_text(size=20), axis.title=element_text(size=15)) +
  geom_text(x=16, y=max(proteoData$Proteolysis), label = lm_eqn(proteoData$Proteolysis, proteoData$Time, proteoData), parse = TRUE, cex= 5)

REG <- ggarrange(growPlot, lipoPlot, proteoPlot, labels = c("A", "B", "C"), ncol = 2, nrow = 2)

png("reg.png", width = 1300, height = 1300, units = "px", pointsize = 12, bg = "white", res = 120)
#svg("reg.svg", width = 10, height = 10)
REG
dev.off()
