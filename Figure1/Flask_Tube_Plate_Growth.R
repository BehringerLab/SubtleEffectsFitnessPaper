
##Requires GrowthCurve Package from Git Hub to plot curves (https://github.com/LennonLab/Growth_Curves)
library(dplyr)
library(growthcurver)
library(purrr)
library(reshape2)
library(ggplot2)
library(ggimage)
library(rstatix)
library(bblme)
library(cowplot)

setwd("~/Desktop/Github/Growth_Curves/output")

# Load Dependencies and Growth Curve R function
source("../bin/modified_Gomp.R")
source("../bin/read.synergy.R")

###Flask Data adapted to simulate BioTek Output

growth.modGomp("<PATH TO>/TrialFlaskAdapted.txt", "TrialFlask", skip = 86)


#Smooth Flask Curves
Param<-read.table("<PATH TO>/TrialFlask.txt",sep=",",header=TRUE)
time<-list(a=1:16)
Param[r,2]

for (r in 1:nrow(Param)){
  my.gomp <- function(t){
    b0 <-   Param[r,2]						
    A <-    Param[r,3]			
    umax <- Param[r,4]						
    L <-    Param[r,5]
    b0+A*exp(-exp(umax*exp(1)*(L-t)/A+1))
  }
  outfile<-lapply(time,my.gomp)
  write.table(as.data.frame(outfile$a),file=paste("<PATH TO>/Flask_Param",r,".txt"))
}

Curve_Fits_Flask<-read.table("<PATH TO>Flask_Param_Tab.txt",sep="\t",header = TRUE)
Curve_Fits_Flask

Curve_Fits_Flask_melt<-reshape2::melt(Curve_Fits_Flask,id="Time")

WT<-Curve_Fits_Flask_melt$variable %in% c("WT_1","WT_2","WT_3")
Curve_Fits_Flask_melt$Background[WT] <- "WT"

Rho<-Curve_Fits_Flask_melt$variable %in% c("Rho_1","Rho_2","Rho_3")
Curve_Fits_Flask_melt$Background[Rho] <- "M1"

RY<-Curve_Fits_Flask_melt$variable %in% c("RY_1","RY_2","RY_3")
Curve_Fits_Flask_melt$Background[RY] <- "M1+2"

YdcI<-Curve_Fits_Flask_melt$variable %in% c("YdcI_1","YdcI_2","YdcI_3")
Curve_Fits_Flask_melt$Background[YdcI] <- "M2"

FitFlaskPlot<-ggplot(data=Curve_Fits_Flask_melt, aes(x=Time, y=value,col=Background))+
  stat_summary(fun.y=mean, geom="point", pch=20,size=1)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",width=0.25)+
  stat_summary(fun.y=mean, geom="line", size=1)+
  xlab("Time (h)")+
  ylab("Absorbance (600nm)")+
  scale_y_continuous(limits=c(0,1.5))+
  scale_x_continuous(limits=c(0,15))+
  ggtitle("Culture flask")+
  scale_color_manual(values=c("#AA0000","#AA00AA","#0000AA","#ABABAB"),name="Strain")+
  theme_bw()+theme(legend.position = "none",legend.background = element_rect(fill = "white", color = "black"),legend.text = element_text(size=8),legend.title = element_text(size=8, face="bold"),axis.title = element_text(size=10,face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(size=10,face="bold",hjust = 0.1, vjust = -8))

FitFlaskPlot
##Calculate Flask growth Parameters
FlaskGC<-read.synergy("<PATH TO>/TrialFlaskAdapted.txt",skip = 86)
FlaskGC.sub = subset(FlaskGC, select = -c(Temp))
FlaskGC.sub.dup<-FlaskGC.sub
FlaskGC.summarize<-SummarizeGrowthByPlate(FlaskGC.sub)
FlaskGC.sub.melt<-reshape2::melt(FlaskGC.sub, id="Time")

FlaskGC.sub.dup$WT<-(FlaskGC.sub.dup$A1+FlaskGC.sub.dup$A2+FlaskGC.sub.dup$A3)/3
FlaskGC.sub.dup$M1<-(FlaskGC.sub.dup$A4+FlaskGC.sub.dup$A5+FlaskGC.sub.dup$A6)/3
FlaskGC.sub.dup$M2<-(FlaskGC.sub.dup$A7+FlaskGC.sub.dup$A8+FlaskGC.sub.dup$A9)/3
FlaskGC.sub.dup$M1_2<-(FlaskGC.sub.dup$A10+FlaskGC.sub.dup$A11+FlaskGC.sub.dup$A12)/3
FlaskGC.average = subset(FlaskGC.sub.dup, select = c(Time,WT,M1,M2,M1_2))
FlaskGC.average.melt<-reshape2::melt(FlaskGC.average, id="Time")


WT<-FlaskGC.summarize$sample %in% c("A1","A2","A3")
FlaskGC.summarize$Background[WT] <- "WT"

Rho<-FlaskGC.summarize$sample %in% c("A4","A5","A6")
FlaskGC.summarize$Background[Rho] <- "M1"

RY<-FlaskGC.summarize$sample %in% c("A10","A11","A12")
FlaskGC.summarize$Background[RY] <- "M1+2"

YdcI<-FlaskGC.summarize$sample %in% c("A7","A8","A9")
FlaskGC.summarize$Background[YdcI] <- "M2"

FlaskGC.summarize

comp <- list(c("WT", "M1"), c("WT", "M2"), c("WT","M1+2"))

pairwise_t_test(FlaskGC.summarize,r~Background,comparisons=comp)
pairwise_t_test(FlaskGC.summarize,k~Background,comparisons=comp)
pairwise_t_test(FlaskGC.summarize,auc_e~Background,comparisons=comp)

###Tube Data adapted to simulate BioTek Output

growth.modGomp("<PATH TO>/TrialTubeAdapted.txt", "TrialTube", skip = 86)

##Smooth Tube Curves
TubeParam<-read.table("<PATH TO>/TrialTube.txt",sep=",",header=TRUE)
time<-list(a=1:16)
TubeParam[r,2]


for (r in 1:nrow(TubeParam)){
  my.gomp <- function(t){
    b0 <-   TubeParam[r,2]						
    A <-    TubeParam[r,3]			
    umax <- TubeParam[r,4]						
    L <-    TubeParam[r,5]
    b0+A*exp(-exp(umax*exp(1)*(L-t)/A+1))
  }
  outfile<-lapply(time,my.gomp)
  write.table(as.data.frame(outfile$a),file=paste("<PATH TO>/Tube_Param",r,".txt"))
}

Curve_Fits_Tube<-read.table("<PATH TO>/Tube_Param_Tab.txt",sep="\t",header = TRUE)
Curve_Fits_Tube

Curve_Fits_Tube_melt<-reshape2::melt(Curve_Fits_Tube,id="Time")

WT<-Curve_Fits_Tube_melt$variable %in% c("WT_1","WT_2","WT_3")
Curve_Fits_Tube_melt$Background[WT] <- "WT"

Rho<-Curve_Fits_Tube_melt$variable %in% c("Rho_1","Rho_2","Rho_3")
Curve_Fits_Tube_melt$Background[Rho] <- "M1"

RY<-Curve_Fits_Tube_melt$variable %in% c("RY_1","RY_3")
Curve_Fits_Tube_melt$Background[RY] <- "M1+2"

YdcI<-Curve_Fits_Tube_melt$variable %in% c("YdcI_1","YdcI_2","YdcI_3")
Curve_Fits_Tube_melt$Background[YdcI] <- "M2"
Curve_Fits_Tube_melt
FitTubePlot<-ggplot(data=Curve_Fits_Tube_melt, aes(x=Time, y=value,col=Background))+
  #geom_point()+
  stat_summary(fun.y=mean, geom="point", pch=20,size=1)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",width=0.25)+
  stat_summary(fun.y=mean, geom="line", size=1)+
  xlab("Time (h)")+
  ylab("Absorbance (600nm)")+
  scale_y_continuous(limits=c(0,1.5))+
  scale_x_continuous(limits=c(0,15))+
  ggtitle("Culture tube")+
  scale_color_manual(values=c("#AA0000","#AA00AA","#0000AA","#ABABAB"))+
  theme_bw()+theme(legend.position = "none",axis.title = element_text(size=10,face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(size=10,face="bold",hjust = 0.1, vjust = -8))
FitTubePlot
##Calculate Tube growth Parameters
TubeGC<-read.synergy("<PATH TO>/TrialTubeAdapted.txt",skip = 86)

TubeGC.sub = subset(TubeGC, select = -c(Temp))
TubeGC.summarize<-SummarizeGrowthByPlate(TubeGC.sub)
TubeGC.sub.melt<-reshape2::melt(TubeGC.sub, id="Time")

TubeGC.sub.dup<-TubeGC.sub
TubeGC.sub.dup$WT<-(TubeGC.sub.dup$A1+TubeGC.sub.dup$A2+TubeGC.sub.dup$A3)/3
TubeGC.sub.dup$M1<-(TubeGC.sub.dup$A4+TubeGC.sub.dup$A5+TubeGC.sub.dup$A6)/3
TubeGC.sub.dup$M2<-(TubeGC.sub.dup$A7+TubeGC.sub.dup$A8+TubeGC.sub.dup$A9)/3
TubeGC.sub.dup$M1_2<-(TubeGC.sub.dup$A10+TubeGC.sub.dup$A11+TubeGC.sub.dup$A12)/3
TubeGC.average = subset(TubeGC.sub.dup, select = c(Time,WT,M1,M2,M1_2))
TubeGC.average.melt<-reshape2::melt(TubeGC.average, id="Time")


WT<-TubeGC.summarize$sample %in% c("A1","A2","A3")
TubeGC.summarize$Background[WT] <- "WT"

Rho<-TubeGC.summarize$sample %in% c("A4","A5","A6")
TubeGC.summarize$Background[Rho] <- "M1"

RY<-TubeGC.summarize$sample %in% c("A10","A11","A12")
TubeGC.summarize$Background[RY] <- "M1+2"

YdcI<-TubeGC.summarize$sample %in% c("A7","A8","A9")
TubeGC.summarize$Background[YdcI] <- "M2"


pairwise_t_test(TubeGC.summarize,r~Background,comparisons=comp)
pairwise_t_test(TubeGC.summarize,k~Background,comparisons=comp)
pairwise_t_test(TubeGC.summarize,auc_e~Background,comparisons=comp)


#####GrowthPlate
PlateGC <- read.synergy("<PATH TO>/LB_GrowthCurvePlate.txt",skip=86)
PlateGC.sub = subset(PlateGC, select = -c(Temp))
PlateGC.summarize<-SummarizeGrowthByPlate(PlateGC.sub)
PlateGC.sub.melt<-reshape2::melt(PlateGC.sub, id="Time")

PlateGC.sub.dup<-PlateGC.sub
PlateGC.sub.dup$WT<-(PlateGC.sub.dup$A1+PlateGC.sub.dup$A2+PlateGC.sub.dup$A3)/3
PlateGC.sub.dup$M1<-(PlateGC.sub.dup$A4+PlateGC.sub.dup$A5+PlateGC.sub.dup$A6)/3
PlateGC.sub.dup$M2<-(PlateGC.sub.dup$A7+PlateGC.sub.dup$A8+PlateGC.sub.dup$A9)/3
PlateGC.sub.dup$M1_2<-(PlateGC.sub.dup$A10+PlateGC.sub.dup$A11+PlateGC.sub.dup$A12)/3
PlateGC.average = subset(PlateGC.sub.dup, select = c(Time,WT,M1,M2,M1_2))
PlateGC.average.melt<-reshape2::melt(PlateGC.average, id="Time")

FitPlatePlot<-ggplot(PlateGC.average.melt, aes(x=Time,y=value,col=variable))+
  geom_point(pch=20,size=1)+
  geom_line(size=1)+
  xlab("Time (h)")+
  ylab("Absorbance (600nm)")+
  ggtitle("96-well plate")+
  scale_y_continuous(limits=c(0,1.5))+
  scale_x_continuous(limits=c(0,15))+
  scale_color_manual(values=c("#ABABAB","#AA0000","#0000AA","#AA00AA"))+
  theme_bw()+theme(legend.position = "none",legend.background = element_rect(fill = "white", color = "black"),legend.text = element_text(size=8),legend.title = element_text(size=8, face="bold"),axis.title = element_text(size=10,face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(size=10,face="bold",hjust = 0.1, vjust = -8))


WT<-PlateGC.summarize$sample %in% c("A1","A2","A3")
PlateGC.summarize$Background[WT] <- "WT"

Rho<-PlateGC.summarize$sample %in% c("A4","A5","A6")
PlateGC.summarize$Background[Rho] <- "M1"

RY<-PlateGC.summarize$sample %in% c("A10","A11","A12")
PlateGC.summarize$Background[RY] <- "M1+2"

YdcI<-PlateGC.summarize$sample %in% c("A7","A8","A9")
PlateGC.summarize$Background[YdcI] <- "M2"


pairwise_t_test(PlateGC.summarize,r~Background,comparisons=comp)
pairwise_t_test(PlateGC.summarize,k~Background,comparisons=comp)
pairwise_t_test(PlateGC.summarize,auc_e~Background,comparisons=comp)




GrowthCurvePlots_Review_Top<-plot_grid(FitPlatePlot,FitFlaskPlot,FitTubePlot,nrow=1, labels = c("a","",""))

GrowthCurvePlots_Review_Top


RY_Fit_D1<-read.table("<PATH TO>/RY_Day1Fit.txt",sep="\t",header = TRUE)
RY_Fit_D1$Competition<-factor(RY_Fit_D1$Competition, levels=c("M1","M2","M1+2"))

RY_Fit_D1 %>% group_by(Competition)%>% summarise(mean.W=mean(W_Flask))
RY_Fit_D1 %>% group_by(Competition)%>% summarise(mean.W=mean(W_Tubes))

M1_Fit<-RY_Fit_D1[which(RY_Fit_D1$Competition=="M1"),]
M2_Fit<-RY_Fit_D1[which(RY_Fit_D1$Competition=="M2"),]
M12_Fit<-RY_Fit_D1[which(RY_Fit_D1$Competition=="M1+2"),]

t.test(M1_Fit$W_Flask,M1_Fit$Control)
t.test(M12_Fit$W_Flask,M12_Fit$Control)
t.test(M2_Fit$W_Flask,M2_Fit$Control)

t.test(M1_Fit$W_Tubes,M1_Fit$Control)
t.test(M12_Fit$W_Tubes,M12_Fit$Control)
t.test(M2_Fit$W_Tubes,M2_Fit$Control)

M12_Fit


RelFitTube<-ggplot(data=RY_Fit_D1, aes(x=Competition, y=W_Tubes,col=Competition))+
  geom_point(alpha=0.3)+
  stat_summary(fun.y=mean, geom="point", pch=21,size=2)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",width=0.25,size=1)+
  geom_hline(yintercept = 1,linetype="dashed",color="black")+
  xlab("")+
  ylab("Relative Fitness")+
  ggtitle("Culture tube")+
  scale_y_continuous(limits=c(0.7,1.4))+
  scale_color_manual(values=c("#AA0000","#0000AA","#AA00AA"))+
  theme_bw()+theme(legend.position = "none",axis.title = element_text(size=10,face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(size=10,face="bold",hjust = 0.1, vjust = -8))

RelFitFlask<-ggplot(data=RY_Fit_D1, aes(x=Competition, y=W_Flask,col=Competition))+
  geom_point(alpha=0.3)+
  stat_summary(fun.y=mean, geom="point", pch=21,size=2)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",width=0.25,size=1)+
  geom_hline(yintercept = 1,linetype="dashed",color="black")+
  xlab("")+
  ylab("Relative Fitness")+
  ggtitle("Culture flask")+
  scale_y_continuous(limits=c(0.7,1.4))+
  scale_color_manual(values=c("#AA0000","#0000AA","#AA00AA"))+
  theme_bw()+theme(legend.position = "none",axis.title = element_text(size=10,face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(size=10,face="bold",hjust = 0.1, vjust = -8))

FitTable<-ggdraw() + 
  draw_image("<PATH TO>/Fig1Table.png")

GrowthCurvePlots_Review_Bottom<-plot_grid(FitTable,RelFitFlask,RelFitTube,nrow=1, rel_widths = c(1,0.5,0.5), labels=c("","b",""))
GrowthCurvePlots_Review_Bottom
plot_grid(GrowthCurvePlots_Review_Top,GrowthCurvePlots_Review_Bottom,nrow=2,rel_heights = c(1,0.85))

