####figure2D 
###Author Jinpei Lin
###Times 2024.09
###R4

library(ggplot2)
m <- readRDS("2D.RDS")
mycolor <- c('stage1'='#a7d3d4','stage2'='#009b9e','stage3'='#e4c1d9','stage4'='#c75dab')

####point plot
ggplot(data = m, aes(gsva_p_all,M_M_SMC))+
geom_point(size=4.5,aes(color=stage_ud))+ 
scale_color_manual(values = mycolor)+
  theme_bw()+theme(panel.grid = element_blank())+
  labs(title= "plaue_score",x="plaue_score",y="Mac_Msmc_SMC")+theme_classic()+theme(axis.text = element_text (size = 20))+
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+
    theme(legend.title = element_text(size=20),
          legend.text = element_text(size=15))

######Ternary plot
library(ggtern)
library(RColorBrewer)
library(grid)
library(scales)


ggtern(data=m,aes(x=SMC,y=Macrophage,z=Modulated_SMC))+geom_mask()+
  geom_point(aes(size=gsva_p_all,color=stage_ud),alpha=0.8)+
  scale_colour_manual(values = mycolor)+
  guides(color = guide_legend(override.aes = list(size = 4)))+
  labs(title = "Ternary stage cell")+theme_rgbw()+theme_custom(col.T = "#df65b0",col.L ="#1f78b4" ,col.R = "#A1CFFA",
  tern.panel.background = "white")