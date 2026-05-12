###Author Jinpei Lin
###Times 2024.07

library(raster)
library(Seurat)
dgz <- "/DATA/User/liqian7/pienapple/00AS_project/DATA_all/00human_data/00gem/"
dout <- "/DATA/User/liqian7/pienapple/00AS_project/DATA_all/00human_data/02RDS_raw/"
for(s in Sys.glob(file.path(dgz,"*.gem.gz"))){
    name <- gsub(".gem.gz","",basename(s))
    input <- paste0(dgz,name,".gem.gz")
    binsize <- 50
    output <- dout
    source(file.path("/DATA/User/liqian7/pienapple/00AS_project/09gem_to_bin50/bin/bin1_to_bin50.R"))
    seurat_spatialObj <- LoadBGI_Spatial(input, outdir = output, bin_data = T, bin_size = binsize, cell_mask = F, area_mask = F, pro_name = name)
}