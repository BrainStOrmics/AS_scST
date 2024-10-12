###3h ECM deg
###Author Jinpei Lin
###Times 2024.09
###R4.3
library(scRNAtoolVis)
m <- readRDS("3h.RDS")
jjVolcano(diffData = m,
          tile.col = c("#DCFFEB","#88C7ED","#ED7A90"),
          size  = 3.5,topGeneN = 5,
          fontface = 'italic',pSize=2.5,legend.position=c(0.9,0.97))