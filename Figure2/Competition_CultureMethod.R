library(ggplot2)
library(ggpubr)
library(reshape2)
library(cowplot)
library(dbplyr)
library(scales)
library(magick)


scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}

M12_Comp_R<-read.table("<PATH_TO>/FlaskTubeFit.txt",sep="\t",header=TRUE,check.names = FALSE)
M12_Comp_R
M12_Comp_R$Comp2<-factor(M12_Comp_R$Comp2,levels=c("Rho","ydcI","Rho_ydcI"))
M12_Comp_R$Competition<-factor(M12_Comp_R$Competition,levels=c("WT_Rho","WFL_ydcI","WT_RY","ydcI_RY"),labels=c("WT vs M1","WT vs M2","WT vs M1+2","M2 vs M1+2"))

M12_R$Method<-factor(M12_Comp_R$Method,levels=c("SingleFlask","SingleTube","MultipleTubes"))
M12_Comp_R$Replicate<-as.factor(M12_Comp_R$Replicate)

as.data.frame(M12_Comp_R %>% group_by(Competition, Day) %>% summarise(strain1=var(Strain_1_CFU),strain2=var(Strain_2_CFU)))
TukeyHSD(aov(data=M12_Comp_R, Strain_2_CFU~Method))

M12_CompFit_Flask<-ggplot(data=M12_Comp_R,aes(shape=Replicate))+
  geom_point(aes(x=Day, y=Strain_2_CFU, col=Comp2,group=Comp2),alpha=0.2)+
  geom_point(aes(x=Day, y=Strain_1_CFU,col=Comp1,group=Comp1),alpha=0.2)+
  geom_line(aes(x=Day, y=Strain_2_CFU, col=Comp2),alpha=0.2, linetype="dashed")+
  geom_line(aes(x=Day, y=Strain_1_CFU,col=Comp1),alpha=0.2, linetype="dashed")+
  stat_summary(aes(x=Day, y=Strain_2_CFU, col=Comp2,group=Comp2),fun.y=mean, geom="line", linetype="dashed")+
  stat_summary(aes(x=Day, y=Strain_2_CFU,col=Comp2,group=Comp2),fun.data="mean_se",geom="errorbar", size=0.5, width=0.2)+
  stat_summary(aes(x=Day, y=Strain_1_CFU, col=Comp1,group=Comp1),fun.y=mean, geom="line", linetype="dashed")+
  stat_summary(aes(x=Day, y=Strain_1_CFU,col=Comp1,group=Comp1),fun.data="mean_se",geom="errorbar", size=0.5, width=0.2)+
  scale_color_manual(values =c("#DE2D26","#660066","#888888","#08306B"))+
  scale_y_log10(label=scientific)+
  scale_x_continuous(limits=c(0,14))+
  ylab("CFU/mL")+
  xlab("Time (d)")+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5,legend.position="bottom", strip.text =element_text(size=8,face="bold"),strip.background =element_rect(color="white", fill="white"),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
  facet_grid(Method~Competition)


M12_s<-read.table("<PATH TO>/Flask_Tube_selection.txt",sep="\t",header=TRUE,check.names = FALSE)
M12_s$Comp2<-factor(M12_s$Comp2,levels=c("Rho","ydcI","Rho_ydcI"))
M12_s$Competition<-factor(M12_s$Competition,levels=c("WT_Rho","WFL_ydcI","WT_RY","ydcI_RY"),labels=c("WT vs M1","WT vs M2","WT vs M1+2","M2 vs M1+2"))
M12_s$Method<-factor(M12_s$Method,levels=c("SingleFlask","SingleTube","MultipleTubes"))
M12_s$Replicate<-as.factor(M12_s$Replicate)
M12_s$Day<-factor(M12_s$Day,levels=c("Day 1","Day 4","Day 14"))
M12_s_plot<- M12_s[which(M12_s$Day=="Day 4"|M12_s$Day=="Day 14"),]

as.data.frame(compare_means(data=M12_s_plot,s~Method,method="t.test",p.adjust.method = "fdr", group.by=c("Competition","Day")))

M12_selectionPlot<-ggplot(data=M12_s_plot,aes(shape=Replicate))+
  geom_point(aes(x=Method, y=s, col=Competition,group=Method),alpha=0.2)+
  stat_summary(aes(x=Method, y=s, col=Competition,group=Method),fun.y=mean, geom="point", pch=21)+
  stat_summary(aes(x=Method, y=s,col=Competition,group=Method),fun.data="mean_se",geom="errorbar", size=0.5, width=0.2)+
   scale_color_manual(values =c("#AA0000","#0000AA","#AA00AA","#660066"))+
  geom_hline(yintercept = 0,color="#888888",linetype="dashed")+
  ylab("Selection rate (s)")+
  xlab("Competition")+
  theme_bw()+theme(axis.title = element_text(face="bold",size=8),legend.background = element_rect(color="black"), legend.title = element_text(face="bold",size=8), legend.title.align=0.5,legend.position="bottom", strip.text =element_text(size=8, face="bold"),strip.background =element_rect(color="white", fill="white"), axis.text.x = element_text(angle=35,hjust=1,size=7),panel.grid.major = element_blank(),panel.grid.minor = element_blank() )+
  theme(legend.position = "none")+
  facet_grid(Day~Competition)


SingleFlask<-ggdraw() + 
  draw_image("<PATH TO>/SingleFlask_Fig2.png")

SingleTube<-ggdraw() + 
  draw_image("<PATH TO>/SingleTube_Fig2.png")

MultipleTubes<-ggdraw() + 
  draw_image("<PATH TO>/MultipleTube_Fig2.png")

MethodCulture<-plot_grid(SingleFlask,SingleTube,MultipleTubes, nrow=3)

plot_grid(MethodCulture,M12_CompFit_Flask,M12_selectionPlot,nrow=1,rel_widths = c(0.3,1,0.7),labels = c("a","b","c"))
