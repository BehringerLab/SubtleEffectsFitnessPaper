library(ggplot2)
library(ggthemes)
Late403<-read.table("<PATH TO>/LateBloom403_1.txt",sep="\t",header=TRUE)

Late403$Clone<-factor(Late403$Clone, levels=c("403-P","403-1","403-12","403-15","403-125","403-2","403-5","403-25"))
Late403$Clone1<-factor(Late403$Clone1, levels=c("yes","no"))

late403.labs <- c("Population","Clone 1","Clone 1+2","Clone 1+3","Clone 1+2+3","Clone 2","Clone 3","Clone 2+3")
names(late403.labs) <- c("403-P","403-1","403-12","403-15","403-125","403-2","403-5","403-25")

LateBloomPlot<-ggplot(data=Late403,aes(x=CultureDay,y=proportion,col=Clone1))+
  geom_point(alpha=0.3)+
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar")+
  scale_x_continuous(breaks=c(0,1,4,10))+
  ylab("Proportion(Late Colonies)")+
  xlab("Time (d)")+
  labs(color = "Contains Clone 1")+
  scale_color_wsj(palette = "colors6")+
  theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "top",strip.background.y = element_blank(),strip.text.y = element_blank(),axis.title = element_text(size=8,face="bold"),legend.title = element_text(size=8,face="bold"),strip.text =element_text(size=8, face="bold"),strip.background =element_rect(color="white", fill="white"))+
  facet_grid(~Clone,labeller = labeller(Clone = late403.labs))
Fig5_Top<-plot_grid(LateBloomPlot, labels = c("a"))

LateBloom_24h<-ggdraw() + 
  draw_image("<PATH TO>/LateBloomPhoto1.png")

LateBloom_48h<-ggdraw() + 
  draw_image("<PATH TO>/LateBloomPhoto2.png")

Fig5_Bottom<-plot_grid(LateBloom_24h,LateBloom_48h, labels = c("b","c"))


plot_grid(Fig5_Top,Fig5_Bottom,nrow=2,rel_heights = c(1,0.7))
