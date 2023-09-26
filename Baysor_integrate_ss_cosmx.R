#need to add selection of comon gene

#this method worked best for me
fibs$named_new <- fibs$clusters
mapjagdata_seurat<-CreateSeuratObject(Seurat::GetAssayData(object = fibs, slot = "counts"), meta.data = fibs@meta.data)
cosmx_seurat<-CreateSeuratObject(Seurat::GetAssayData(object = fibs_s, slot = "counts"), meta.data = fibs_s@meta.data)

head(colnames(mapjagdata$counts))
head(rownames(mapjagdata$metadata))

median(cosmx_seurat$nCount_RNA)
median(mapjagdata_seurat$nCount_RNA)

cosmx_seurat <- cosmx_seurat %>% 
  NormalizeData(scale.factor = 47) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()

mapjagdata_seurat <- mapjagdata_seurat %>% 
  NormalizeData(scale.factor = 9170) %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA() 
cosmx_seurat$ID <- "cosmx"
mapjagdata_seurat$ID <- 'mapjag'

features <- SelectIntegrationFeatures(object.list = c(cosmx_seurat,mapjagdata_seurat ))
int.anch <- FindIntegrationAnchors(object.list = c(cosmx_seurat,mapjagdata_seurat ), anchor.features = features)
int <- IntegrateData(anchorset = int.anch)

int <- int %>% 
  ScaleData() %>% 
  FindVariableFeatures() %>% 
  RunPCA()  %>% 
  RunUMAP(dims=1:40)

DimPlot(int, group.by = "ID")

int_H<-RunHarmony.Seurat_CM(int, group.by.vars = c("named_new"), max.iter.harmony = 40, max.iter.cluster = 70 )
int_H <- RunUMAP(int_H, reduction="harmony", dims=1:40)
DimPlot(int_H, group.by = "named_new")
