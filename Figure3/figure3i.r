####3G msmc vs fibro deg base on stage3 and stage4
###Author Jinpei Lin
###Times 2024.09
###R4.3
library(Seurat)
library(EnhancedVolcano) 
m <- readRDS("3g.RDS")
deg = FindMarkers(m,ident.1 = 'MSMC_100',
            ident.2 = 'Fibro_25')
head(deg[order(deg$p_val),])
EnhancedVolcano(deg,gridlines.major = FALSE,gridlines.minor = FALSE,
                lab = rownames(deg),
                x = 'avg_log2FC',
                y = 'p_val_adj')