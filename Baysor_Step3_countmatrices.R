samples<-1:60
smaples<-paste("/rds/projects/c/croftap-celldive01/CosmX_analysis/basor_segmentation_from_raw/S2/baysor_outs/FOV", samples, ".csv", sep="")

smaples_area<-paste("/rds/projects/c/croftap-celldive01/CosmX_analysis/basor_segmentation_from_raw/S2/baysor_outs/FOV", samples, "_cell_stats.csv", sep="")

library(splitstackshape)
1:length(smaples)
data = list()
for (i in 1:length(smaples)) {
    data[[i]] <- read.csv(smaples[[i]]);
    data[[i]]<-data[[i]][data[[i]]$is_noise== "false",];
    data[[i]]<-cSplit(data[[i]], splitCols = "cell", sep="-")
}

data_plot<-rbindlist(data)

#calculate area of cells (can be used to re run baysor or attach to seurat obj for further metrics)
data_area = list()
for (i in 1:length(smaples)) {
    data_area[[i]] <- read.csv(smaples_area[[i]]);
    data_area[[i]]<-cSplit(data_area[[i]], splitCols = "cell", sep="-");
    data_area[[i]]<-data_area[[i]][data_area[[i]]$cell_2 %in%  data[[i]]$cell_2];
    data_area[[i]]$radius<-sqrt(data_area[[i]]$area/3.14);
    data_area[[i]]<-data_area[[i]][data_area[[i]]$radius > 20 & data_area[[i]]$radius < 60, ]
}


for (i in 1:length(smaples)) {
   data[[i]]<-data[[i]][data[[i]]$cell_2 %in%  data_area[[i]]$cell_2]
}

#turn tx counts to count matrices
tx_to_counts <- function(genes, cells, remove_bg = TRUE) {
    if (remove_bg) {
        idx <- which(cells != 0)
        cells <- cells[idx]
        genes <- genes[idx]
    }
    genes <- factor(genes)
    cells <- factor(cells)
    counts <- Matrix::sparseMatrix(
        i = as.integer(genes), 
        j = as.integer(cells), 
        x = rep(1, length(genes)),
        dims = c(length(levels(genes)), length(levels(cells)))
    )
    rownames(counts) <- levels(genes)
    colnames(counts) <- levels(cells)
    return(counts)
}

#this might take some time to run
count_matrices=list()
for (i in 1:length(smaples)) {
  count_matrices[[i]]<-tx_to_counts(genes =data[[i]]$gene, data[[i]]$cell_2)
  }

smaples_fov<-paste("FOV", samples, sep="")
seurat_obj<-list()
for (i in 1:length(smaples)) {
  seurat_obj[[i]]<-CreateSeuratObject(counts =count_matrices[[i]], project = smaples_fov[[i]])
}

for (i in 1:length(smaples)) {
  seurat_obj[[i]]$radius<-data_area[[i]]$radius
}

slide2_merged<-merge(x=seurat_obj[[1]], y=seurat_obj[c(2:60)])


hist(slide2_merged$radius)

max(slide2_merged$radius)

#follow normal seurat processing to annotate

median(slide2_merged$nCount_RNA)  #19
median(slide2_merged$nFeature_RNA)  #17

hist(slide2_merged$nCount_RNA)  
hist(slide2_merged$nFeature_RNA)

slide2_merged_f <- subset(slide2_merged, subset = nFeature_RNA > 10 & nCount_RNA > 10)
median(slide2_merged$nCount_RNA)

slide2_merged_f <- slide2_merged_f %>% 
  NormalizeData(scale.factor = 19) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() %>% 
  RunUMAP(dims=1:40)



#plot FOV txs
library(splitstackshape)
1:length(smaples)
data = list()
for (i in 1:length(smaples)) {
    data[[i]] <- read.csv(smaples[[i]]);
    data[[i]]<-data[[i]][data[[i]]$is_noise== "false",];
    data[[i]]<-cSplit(data[[i]], splitCols = "cell", sep="-")
}

data_plot<-rbindlist(data)

data_plot[data_plot$fov %in% c("40"),] %>% ggplot(aes(x=x,y=y)) + 
  geom_point(alpha=0.3, shape=".") +
  geom_point(data=data_plot[data_plot$fov %in% c("40") & data_plot$gene %in% c("VWF", "IGFBP7", "PECAM1", "CDH5"),],
             aes(x=x,y=y, color=gene), 
                   size=1)+ scale_color_brewer(palette="Dark2")

