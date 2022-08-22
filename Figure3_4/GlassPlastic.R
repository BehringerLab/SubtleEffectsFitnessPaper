library(ggplot2)
library(scales)
library(cowplot)
library(ggpubr)

review_pH<-read.table("<PATH TO>/GlassPlastic_10d_pH.txt",sep="\t",header=TRUE)

pH.labs <- c("Clone 1","Clone 2","Clone 3")
names(pH.labs) <- c("403-1","403-2","403-5")

day.labs <- c("Day 0","Day 1","Day 4", "Day 10")
names(day.labs) <- c("0","1","4","10")

DayPair.labs<- c("Day 1","Day 4", "Day 10")
names(DayPair.labs) <- c("1","4","10")

Fit_plot<-ggplot(data=review_pH, aes(x=CultureDay,col=CultureTube))+
  #geom_point()+
  stat_summary(aes(y=WTCount),fun.y="mean", geom="line", color="#ABABAB")+
  stat_summary(aes(y=WTCount),fun.data = "mean_cl_boot", geom = "errorbar",color="#ABABAB")+
  stat_summary(aes(y=EvoCount),fun.y="mean", geom="line")+
  stat_summary(aes(y=EvoCount),fun.data = "mean_cl_boot", geom = "errorbar")+
  scale_color_brewer(palette = "Dark2")+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
  xlab("Time (d)")+
  ylab("CFU/mL")+
  scale_x_continuous(breaks=c(0,1,4,10))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",,strip.background.y = element_blank(),strip.text.y = element_blank(),axis.title = element_text(size=8,face="bold"),strip.text.x = element_text(face="bold"),strip.background.x =element_rect(color="white", fill="white"))+
  facet_grid(Clone~CultureTube,labeller = labeller(Clone = pH.labs))

pH_plot<-ggplot(data=review_pH, aes(x=CultureTube,fill=CultureTube,y=pH))+
  geom_boxplot(aes(fill=CultureTube))+
  scale_fill_brewer(palette = "Dark2")+
  xlab("Tube Material")+
  scale_x_discrete(labels=c("Glass"="G","Plastic"="P"))+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),,legend.position = "none",,strip.background.y = element_blank(),strip.text.y = element_blank(),axis.title = element_text(size=8,face="bold"),strip.text.x = element_text(face="bold"),strip.background.x =element_rect(color="white", fill="white"))+
  facet_grid(Clone~CultureDay,labeller = labeller(Clone = pH.labs, CultureDay=day.labs))

plot_grid(Fit_plot,pH_plot)

as.data.frame(compare_means(pH~CultureTube, review_pH,group.by = "CultureDay"))

review_selection<-read.table("/Users/megangrace/Box/Behringer_Lab_Box_Drive/Manuscripts/In_Progress/ReviewPaper/GlassPlastic_10d_pH_selection.txt",sep="\t",header=TRUE)

Day10_pH_s<-review_selection[which(review_selection$DayPair==10),]
Day10_pH_s_C1<-Day10_pH_s[which(Day10_pH_s$Clone=="403-1"),]
Day10_pH_s_C2<-Day10_pH_s[which(Day10_pH_s$Clone=="403-2"),]
Day10_pH_s_C3<-Day10_pH_s[which(Day10_pH_s$Clone=="403-5"),]

cor.test(Day10_pH_s_C1$pH,Day10_pH_s_C1$Selection)
cor.test(Day10_pH_s_C2$pH,Day10_pH_s_C2$Selection)
cor.test(Day10_pH_s_C3$pH,Day10_pH_s_C3$Selection)

selectionPlot<-ggplot(data=review_selection, aes(y=pH,x=Selection,col=pH))+
  geom_point(aes(shape=CultureTube))+
  geom_smooth(method="lm", color="black",linetype="dashed", se=FALSE, size=0.5)+
  scale_color_continuous(type = "viridis")+
  #scale_color_gradient(low="blue", high="red")+
  labs(shape = "Tube")+
  xlab("Selection rate (s)")+
  guides(shape = guide_legend(override.aes = list(size = 1)),color = guide_legend(override.aes = list(size = 1)),)+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "right",legend.title = element_text(size = 8), legend.text = element_text(size = 8),legend.key.size = unit(0.1, "lines"),axis.title = element_text(size=8,face="bold"),strip.text = element_text(face="bold"),strip.background =element_rect(color="white", fill="white"))+
  facet_grid(Clone~DayPair,labeller = labeller(Clone = pH.labs, DayPair=DayPair.labs))
selectionPlot

plot_grid(Fit_plot,pH_plot,selectionPlot,nrow=1,rel_widths = c(0.75,0.7,1), labels=c("a","b","c"))
  

still_shake<-read.table("<PATH TO>/StillShake_pH.txt",sep="\t",header=TRUE)

still_shake
still_shake$Material<-factor(still_shake$Material, levels=c("Plastic","Glass","Flask"))
still_shake$Shaking<-factor(still_shake$Shaking, levels=c("Still","Shaking"))


ggplot(data=still_shake,aes(x=Material,y=pH,col=Material))+
  geom_jitter(alpha=0.3)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  scale_color_manual(values=c("#d95f02","#1b9e77","#7570b3"))+
  scale_x_discrete(labels=c("Plastic"="Tube (P)","Glass"="Tube (G)","Flask" = "Flask"))+
  xlab("Culture Vessel")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",axis.title = element_text(size=8,face="bold"),axis.text.x = element_text(angle=22,hjust=1),strip.text = element_text(face="bold"),strip.background =element_rect(color="white", fill="white"))+
  facet_wrap(~Shaking)

TukeyHSD(aov(data=still_shake, pH~Material*Shaking))$`Material:Shaking`



