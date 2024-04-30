#select genes from scRNAseq base don xenium

#read in scRNAseq marker genes
markers_named <- read_csv("/rds/projects/c/croftap-mapjagdata/MAPJAGv2/Chris/markers_named.csv")

#filter
fibs_MAPJAG_marker_genes_f <- markers_named %>% filter(cluster %in% c("CD34/MFAP5"   ,  "POSTN"     ,     "SFRP/CXCL12" ,        "LL" ,          
 "CXCL14"  ,       "SOX5_CDH11" ,    "MMPs"   ) & avg_log2FC >0.5 & p_val_adj < 0.05)


##this step is optional and was not that useful for fibs
#subset xenium to just fibs (remove contaminating clusters)
Idents(fibs_X) <- 'named_sub'
fibs_clean_X <- subset(fibs_X, idents=c("IFD1_COL5A2" , "PDGFRA",  "TNC1_FBN1" ,  "LL" ,  "MFAP5_CD34",  "THY1_MFAP5"))


#keep overlapping marker genes with scRNAseq and genes in xenium
genes_use <- fibs_MAPJAG_marker_genes_f$gene
genes_use <- genes_use[genes_use %in% rownames(fibs_clean_X)]
genes_use <- na.omit(genes_use)


cosmxdata_celltype<-list()
counts<-Seurat::GetAssayData(object = fibs_clean_X, slot = "counts")
cosmxdata_celltype$counts <- counts[genes_use, ]
cosmxdata_celltype$metadata<-fibs_clean_X[[]]
cosmxdata_celltype$metadata$ID<-'xenium'


mapjagdata<-list()
counts <- Seurat::GetAssayData(object = fibs, slot = "counts")
mapjagdata$counts <- counts[rownames(cosmxdata_celltype$counts), ]
mapjagdata$metadata<-fibs@meta.data
mapjagdata$metadata$ID<-'mapjag'


fibs_X_xenium_adj<-CreateSeuratObject(counts=cosmxdata_celltype$counts, meta.data = cosmxdata_celltype$metadata)
fibs_MAPJAG_xenium_adj<-CreateSeuratObject(counts=mapjagdata$counts, meta.data = mapjagdata$metadata)

fibs_X_xenium_adj <- fibs_X_xenium_adj %>% 
  NormalizeData(scale.factor = median(fibs_X_xenium_adj$nCount_RNA)) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()

fibs_MAPJAG_xenium_adj <- fibs_MAPJAG_xenium_adj %>% 
  NormalizeData(scale.factor = median(fibs_MAPJAG_xenium_adj$nCount_RNA)) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() 
fibs_X_xenium_adj$ID <- "xenium"
fibs_MAPJAG_xenium_adj$ID <- 'mapjag'

fibs_X_xenium_adj$clusters <- fibs_X_xenium_adj$named_sub

features <- SelectIntegrationFeatures(object.list = c(fibs_MAPJAG_xenium_adj,fibs_X_xenium_adj ))
int.anch <- FindIntegrationAnchors(object.list = c(fibs_MAPJAG_xenium_adj,fibs_X_xenium_adj ), anchor.features = features)
int <- IntegrateData(anchorset = int.anch)

int <- int %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()  %>% 
  RunUMAP(dims=1:20)

DimPlot(int, group.by = "ID")

DimPlot(int, group.by = "clusters")

int$ID_clusters <- paste(int$ID, int$clusters, sep="_")


int_H<-RunHarmony(int, group.by.vars = c("ID_clusters") )
int_H <- RunUMAP(int_H, reduction="harmony", dims=1:10)
DimPlot(int_H, group.by = "ID_clusters")


DimPlot(int_H, group.by = "clusters") #now check for overlaping/non-overlaping clusters (this was not good for fibs)


####
Idents(fibs_clean_X) <- 'named_sub'
fibs_clean_X_markers <- FindAllMarkers(fibs_clean_X, only.pos = T)
fibs_clean_X_markers <- fibs_clean_X_markers %>% filter(p_val_adj < 0.05)
fibs_MAPJAG_marker_genes_f_xenium_genes <- fibs_MAPJAG_marker_genes_f[fibs_MAPJAG_marker_genes_f$gene %in% rownames(fibs_clean_X),]


fibs_clean_X_markers$cluster %>% unique()


xenium_markers <- list()
xenium_clusters_overlap <- list()
df_list <- list()

for (i in 1:length(unique(fibs_clean_X_markers$cluster))){
xenium_markers[[i]] <- fibs_clean_X_markers %>% filter(cluster == unique(fibs_clean_X_markers$cluster)[[i]])
xenium_clusters_overlap[[i]] <- fibs_MAPJAG_marker_genes_f_xenium_genes[fibs_MAPJAG_marker_genes_f_xenium_genes$gene %in% xenium_markers[[i]]$gene,]

df <- table(xenium_clusters_overlap[[i]]$cluster) %>% as.data.frame()
df2 <- table(fibs_MAPJAG_marker_genes_f_xenium_genes$cluster) %>% as.data.frame()
df <- merge(x = df, y = df2, by = "Var1", all = TRUE)
df$pct <- df$Freq.x/df$Freq.y*100
df_list[[i]] <- df %>% replace(is.na(.), 0)
df_list[[i]]$cluster <- unique(fibs_clean_X_markers$cluster)[[i]]
print(ggplot(df_list[[i]], aes(x=Freq.x, y=pct)) + 
  geom_point()+theme_ArchR()+xlim(0,17)+ geom_text_repel(aes(label = df$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_list[[i]]$Freq.y)/2, linetype='dotted', col = 'red', size=0.5)+geom_vline(xintercept = max(df_list[[i]]$Freq.x)/2, linetype="dotted", 
                color = "red", size=0.5)+ggtitle(unique(fibs_clean_X_markers$cluster)[[i]]))



}

names(xenium_clusters_overlap) <- unique(fibs_clean_X_markers$cluster)
#use this output to name and merge xenium clusters




