##Seurat 3 processing of snRNA data for both 12wks and 8 wks
#Set up
library(dplyr)
library(Seurat)
library(gplots)
library(ggplot2)
library(gridExtra)
setwd("/dfs4/som/seurat3/")
r62str.counts <- Read10X(data.dir = "aggr_all12/outs/filtered_feature_bc_matrix/")
#meta data for samples
meta=read.csv(file="metadata.csv",header=F)
names(meta)=c("ID","name","geno","age","batch")
r62str <- CreateSeuratObject(counts = r62str.counts)
m=data.frame(orig.ident=as.numeric(substr(colnames(r62str),18,19)))
rownames(m)=colnames(r62str)
r62str$orig.ident <- m
x=data.frame(ID=r62str$orig.ident)
y=merge(x,meta,by="ID",all.x=T)
r62str$batch <- y$batch
r62str$geno <- y$geno
r62str$age <- y$age
list=SplitObject(r62str,split.by="age")
#subset to 12wk samples
r62str=list[[1]]
########Individual sample processing using standard Seurat
r62str <- subset(r62str, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 2)
r62str <- NormalizeData(object = r62str)
r62str <- FindVariableFeatures(r62str, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(r62str)
r62str <- ScaleData(r62str, features = all.genes)
r62str <- RunPCA(r62str, features = VariableFeatures(object = r62str))
VizDimLoadings(r62str, dims = 1:2, reduction = "pca")
print(r62str[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(r62str, reduction = "pca")
DimHeatmap(r62str, dims = 1, cells = 500, balanced = TRUE)
r62str <- FindNeighbors(r62str, dims = 1:20)
r62str <- FindClusters(r62str, resolution = 0.05)
r62str <- RunTSNE(object = r62str,dims=1:20)
DimPlot(object = r62str, reduction = "tsne")
r62str <- RunUMAP(r62str, dims = 1:20)
DimPlot(r62str, reduction = "umap")
###############finder markers for each cluster 
for (i in 0:12) {
  x=FindMarkers(r62str,ident.1=i,min.pct=0.25)
  fn=paste0("12wks_cluster",i,"vsall.csv")
  write.csv(as.data.frame(x),file=fn)
}
for (i in 0:12) {
  x=FindMarkers(r62str,subset.ident=i,group.by = "geno",ident.1 = "HD", min.pct=0.25)
  fn=paste0("12wks_cluster",i,"HDvsCTR.csv")
  write.csv(as.data.frame(x),file=fn)
}
###############
#recluster #2, 7 in 12wks for oligo and OPC
r62str1=subset(r62str,idents=c(2,7))
r62str1 <- NormalizeData(object = r62str1)
r62str1 <- FindVariableFeatures(r62str1, selection.method = "vst", nfeatures = 1500)
all.genes <- rownames(r62str1)
r62str1 <- ScaleData(r62str1, features = all.genes)
r62str1 <- RunPCA(r62str1, features = VariableFeatures(object = r62str1))
DimPlot(r62str1, reduction = "pca")
DimHeatmap(r62str1, dims = 1, cells = 500, balanced = TRUE)
r62str1 <- FindNeighbors(r62str1, dims = 1:11)
r62str1$saved.idents=Idents(r62str1)
summary(Idents(r62str1))
r62str1 <- FindClusters(r62str1, resolution = 0.15)
r62str1 <- RunTSNE(object = r62str1,dims=1:11)
DimPlot(object = r62str1, reduction = "tsne")
r62str1 <- RunUMAP(r62str1, dims = 1:11)
DimPlot(r62str1, reduction = "umap")
saved=as.factor(r62str1$saved.idents)
###label clusters
a="Excitatory neurons 1
D2+ MSN
Oligodendrocytes
D1/D2+ de-identified
D1+ MSN
Astrocyte
Inhibitory neurons 1 (foxp2&pbx3)
OPC
Excitatory neurons 2
Microglia
Excitatory neurons 3
Endo, peri, and vSMC
Cholinergic interneurons"
x=strsplit(a,"\n")
new.cluster.ids <- x[[1]]
names(new.cluster.ids) <- levels(r62str)
r62str <- RenameIdents(r62str, new.cluster.ids)
