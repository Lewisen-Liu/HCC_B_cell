#GSE174739 data for Breg, Frontiers in Immunology####

#Fig S14a#####
library(Seurat)
library(dplyr)

counts <- read.delim(".../GSE174739_RAW/GSE174739_filtered_10X_B_counts.txt",
                     header = TRUE, row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE)

meta   <- read.delim(".../GSE174739_RAW/GSE174739_filtered_10X_B_metadata.txt",
                     header = TRUE, row.names = 1,
                     stringsAsFactors = FALSE)

counts <- counts[, rownames(meta)]

srt <- CreateSeuratObject(
  counts   = counts,
  meta.data= meta,
  project  = "GSE174739_B"
)

head(srt)


# seurat function######
set.seed(1)
seurat<-function(x){
  x<-NormalizeData(x,normalization.method = "LogNormalize",scale.factor = 10000)
  x<-FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  x<-ScaleData(x,features = VariableFeatures(x))
  x<-RunPCA(x, features = VariableFeatures(object = x),npcs = 50, verbose = TRUE)
  ElbowPlot(x,ndims = 50)
  x<-FindNeighbors(x, reduction = "pca",dims = 1:30)
  x<-FindClusters(x, resolution = 0.2)
  x<-RunUMAP(x,dims = 1:30,reduction = "pca")
  return(x)
}

srt=seurat(srt)


FeaturePlot(srt, features=c('Cd19','Ms4a1','Il10','Il12a','Ebi3','Havcr1'),pt.size = 1,reduction = "umap",order = T,cols = c("#DFDFDF","#FF0000"),raster=FALSE)
ggsave(".../mouse scRNAseq/Frontiers_mouse_Cd19_and_Tim1.pdf",
       width = 5,height = 6,units = "in")


