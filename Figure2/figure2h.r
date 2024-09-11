####2H
###Author Jinpei Lin
###Times 2024.09
###R4
library(ggplot2)
m <- readRDS("2h.RDS")
cols <- c("Macrophage"='#df65b0',"Modulated_SMC"="#A1CFFA",
           "SMC"='#1f78b4')
ggplot(m,aes(x=gsva,y = percent,fill=celltype,color=celltype))+theme_classic ()+
labs(subtitle = "")+geom_smooth(method = "loess",span=1,se = T,level=0.95,alpha=0.15)+scale_fill_manual(values = cols)+scale_color_manual(values = cols)+
labs(x="Plaque_score",y="Proportion")+theme(axis.text = element_text (size = 20))+#调整坐标轴字体大小
    theme(axis.title.x=element_text(vjust=0, size=25,face = "plain"))+#调整xlab字体大小
    theme(axis.title.y=element_text(vjust=0, size=25,face = "plain"))+#调整ylab字体大小
    theme(plot.title = element_text(size = 25, face = "bold"))+
    theme(axis.line = element_line(color = "black",linewidth = 0.8))+#调整坐标轴样式
    theme(legend.title = element_text(size=20),  # 调整图例标题的大小
          legend.text = element_text(size=15))