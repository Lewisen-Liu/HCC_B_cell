library(SingleCellExperiment)
library(Seurat)
library(dplyr)
library(glmGamPoi)
library(ggsci)
library(ggplot2)
library(ggpubr)

#GSE206325 dataset from GEO database
load(".../GSE206325/GSE206325_data_SingleCellExperiment_object.Rda")
cluster_annotation=read.csv(".../GSE206325/GSE206325_full_HCC_cluster_annotation.csv")
sample_annotation=read.csv(".../GSE206325/GSE206325_sample_annots_Liver_Treated_patients.csv")

# convert sce to seurat file (https://github.com/satijalab/seurat/issues/4556)
df <- as.data.frame(sce@colData)
rownames(df) <- df$barcodes
df <- df[colnames(sce@assays$data@listData[["counts"]]),]
srt <- CreateSeuratObject(sce@assays$data@listData[["counts"]] , meta.data = df)

colnames(sample_annotation)[1]='cell_to_sample_ID'
sample_annotation$cell_to_sample_ID=as.character(sample_annotation$cell_to_sample_ID)
srt@meta.data=left_join(srt@meta.data, sample_annotation, by = "cell_to_sample_ID")
rownames(srt@meta.data)=srt@meta.data$barcodes

useful=unique(srt@meta.data$cell_to_sample_ID) %>% intersect(sample_annotation$cell_to_sample_ID)

srt=srt %>% subset(subset=cell_to_sample_ID %in% useful)
srt
rm(list = c('useful'))

colnames(cluster_annotation)[5]='cell_to_cluster'
cluster_annotation$cell_to_cluster=as.character(cluster_annotation$cell_to_cluster)

srt@meta.data=left_join(srt@meta.data, cluster_annotation, by = "cell_to_cluster")
rownames(srt@meta.data)=srt@meta.data$barcodes # barcode name is necessary for subset function
rm(list = c('cluster_annotation','sample_annotation','df','sce'))

#start analyzing
library(ggplot2)
library(sctransform)

# store mitochondrial percentage in object meta data
srt <- PercentageFeatureSet(srt, pattern = "^MT-", col.name = "percent.mt")
head(srt)

# seurat function
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


mycolor=c('#46A5DA','#FAEC95','#5080CC','#EB4C20','#DEA077','#BF76CC','#8DADC4','#40BFCC','#CC6F62','#CC3E74',
          '#7373EA','#FEDDAA','#6BCC50','#BBB8CC','#D7DE5A','#8A52CC','#CC6823','#E288D1','#48CCAD','#ACCCA3',
          '#C936C2','#E55F99','#B9CEDB','#8F8FCC','#C19FCC','#F22727','#3975AC','#F5BC1D','#F1A9BA')


srt_tumor=srt %>% subset(subset=tissue %in% c('Tumor')) #339590
srt_tumor <- subset(srt_tumor, cells = rownames(srt_tumor@meta.data[is.na(srt_tumor@meta.data$exclude), ])) #336998 cells

rm(list = c('srt'))

set.seed(1)
srt_tumor=seurat(srt_tumor)


#Fig S16b for immune cell prop######

meta <- srt_tumor@meta.data %>% as.data.frame()

totals <- meta %>%
  group_by(patient_ID, treatment_Resp) %>%
  summarise(total_cells = n(), .groups = "drop")

counts <- meta %>%
  group_by(patient_ID, treatment_Resp, group) %>%
  summarise(cells = n(), .groups = "drop")

df_prop <- counts %>%
  left_join(totals, by = c("patient_ID", "treatment_Resp")) %>%
  mutate(prop = cells / total_cells)

desired_order <- c("T","B",  "MonoMac", "NK", 'Stromal', "DC", "ILC",'Mast')  
df_prop$group <- factor(df_prop$group, levels = desired_order)


pvals <- df_prop %>%
  group_by(group) %>%
  summarise(
    p     = t.test(
      prop[treatment_Resp == "anti-PD1_NR"],
      prop[treatment_Resp == "anti-PD1_R"]
    )$p.value,
    y_max = max(prop)
  ) %>%
  ungroup() %>%
  mutate(
    label = paste0("p=", signif(p, 2)),
    ypos  = y_max + 0.05
  )



ggplot(df_prop, aes(x = group, y = prop, fill = treatment_Resp)) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.5
  ) +
  geom_jitter(
    aes(color = treatment_Resp),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
    size = 1.5, alpha = 0.8
  ) +
  stat_summary(
    fun = mean, geom = "point",
    shape = 18, size = 3, color = "black",
    position = position_dodge(width = 0.8)
  ) +
  geom_text(
    data = pvals,
    aes(x = group, y = ypos, label = label),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_fill_manual(
    values = c("anti-PD1_NR" = "#2E9FDF", "anti-PD1_R" = "#FF6B6B")
  ) +
  scale_color_manual(
    values = c("anti-PD1_NR" = "#2E9FDF", "anti-PD1_R" = "#FF6B6B"),
    guide = "none"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(
    title = "anti-PD1_NR vs. anti-PD1_R",
    x     = "Cell Type",
    y     = "Proportion per Patient",
    fill  = "Treatment Response"
  ) +
  theme_bw() +
  theme(
    plot.title  = element_text(size = 12, hjust = 0.5),
    axis.title  = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5,color = 'black'),
    axis.text.y = element_text(size = 10,color = 'black'),
    legend.position = "right"
  )

ggsave(".../human scRNAseq/proportion_compare.pdf",width = 7,height = 4,units = "in")







#B cell######
srt_tumor_Bcell=srt_tumor %>% subset(subset=group %in% c('B')) #19248 cells
srt_tumor_Bcell_noplasma=srt_tumor_Bcell %>% subset(subset=subgroup %in% c('Memory','Naive')) #14137 cells

set.seed(1)
srt_tumor_Bcell_noplasma=seurat(srt_tumor_Bcell_noplasma)


#Fig S17a#####
DimPlot(srt_tumor_Bcell_noplasma, reduction = "umap",label = F, group.by = "subgroup")+scale_color_manual(values = mycolor)
ggsave(".../human_scRNA/B_cell_cluster.pdf",width = 5.5,height = 5,units = "in")

#Fig S17b#####
DimPlot(srt_tumor_Bcell_noplasma, reduction = "umap",label = F, group.by = "treatment_Resp")+scale_color_manual(values = mycolor[c(3,4)])
ggsave(".../human_scRNA/B_cell_cluster_response.pdf",width = 6,height = 5,units = "in")


markers.srt_tumor_Bcell_noplasma<-FindMarkers(srt_tumor_Bcell_noplasma,
                                              ident.1 = "anti-PD1_NR", 
                                              ident.2 = "anti-PD1_R",
                                              group.by='treatment_Resp')


#Fig S14b######
FeaturePlot(srt_tumor_Bcell_noplasma, 
            features=c('CD19','MS4A1','IL10','IL12A','EBI3','HAVCR1'),reduction = "umap",order = T,cols = c("#DFDFDF","#FF0000"),raster=FALSE,pt.size = 1)

ggsave(".../Nature_Medicine_human_Cd19_and_Tim1.pdf",
       width = 5,height = 6,units = "in")


#Fig S17h TIM1+ vs TIM1-######

havcr1_expr <- GetAssayData(srt_tumor_Bcell_noplasma, 
                            assay = "RNA", 
                            slot = "data")["HAVCR1", ]

# divede into 2 groups
srt_tumor_Bcell_noplasma$TIM1_group <- ifelse(havcr1_expr > 0, 
                                              "TIM-1_positive", 
                                              "TIM-1_negative")


VlnPlot(srt_tumor_Bcell_noplasma, 
        features = c("CD40", 'CD74','HLA-DMA','HLA-DMB'), 
        group.by = "TIM1_group", 
        pt.size = 0, ncol = 4)
ggsave(".../antigen_markers.pdf",width = 10,height = 4,units = "in")






#enrichment analysis#####
#Fig S17c,d,e####
library(ggpubr)
library(ggplot2)
library(ggthemes)

markers.srt_tumor_Bcell_noplasma_filter=markers.srt_tumor_Bcell_noplasma %>% filter(p_val_adj<0.05)


markers_up=markers.srt_tumor_Bcell_noplasma_filter %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
markers_down=markers.srt_tumor_Bcell_noplasma_filter %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))

markers_up=markers.srt_tumor_Bcell_noplasma %>% filter(avg_log2FC<0) %>% arrange(avg_log2FC)
markers_down=markers.srt_tumor_Bcell_noplasma %>% filter(avg_log2FC>0) %>% arrange(desc(avg_log2FC))



library(DESeq2)
library(clusterProfiler)
library(msigdbr)
library(stringr)

geneset_GOBP=msigdbr(species = 'Homo sapiens',
                     category = 'C5',
                     subcategory = 'GO:BP') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_GOBP$gs_name=str_to_title(geneset_GOBP$gs_name)
head(geneset_GOBP)

geneset_reactome=msigdbr(species = 'Homo sapiens',
                         category = 'C2',
                         subcategory = 'CP:REACTOME') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_reactome$gs_name=str_to_title(geneset_reactome$gs_name)
head(geneset_reactome)

geneset_KEGG=msigdbr(species = 'Homo sapiens',
                     category = 'C2',
                     subcategory = 'CP:KEGG') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_KEGG$gs_name=str_to_title(geneset_KEGG$gs_name)
head(geneset_KEGG)

geneset_hallmark=msigdbr(species = 'Homo sapiens',
                         category = 'H') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_hallmark$gs_name=str_to_title(geneset_hallmark$gs_name)
head(geneset_hallmark)



#for markers_up
ranks= -markers_up$avg_log2FC
names(ranks)=rownames(markers_up)

#for markers_down
ranks= markers_down$avg_log2FC
names(ranks)=rownames(markers_down)

head(ranks)
ranks=sort(ranks,decreasing = T)

#####GSEA-GOBP#####
set.seed(1)
Res=GSEA(ranks, TERM2GENE = geneset_GOBP,seed = T)
dotplot(Res, x = "NES",showCategory=10)
ggsave(".../up_Gobp_GSEA.pdf",width = 6,height = 7,units = "in")

#####GSEA-REACTOME#####
set.seed(1)
Res=GSEA(ranks, TERM2GENE = geneset_reactome,seed = T)
dotplot(Res, x = "NES",showCategory=10)
ggsave(".../up_Reactome_GSEA.pdf",width = 6,height = 7,units = "in")


#####enricher-KEGG#####
set.seed(1)
enrich_KEGG <- enricher(names(ranks),TERM2GENE=geneset_KEGG)
barplot(enrich_KEGG,showCategory=10)
ggsave(".../up_Kegg_enricher.pdf",width = 6,height = 7,units = "in")




#Fig S17f####
antigen_pathway=geneset_KEGG %>% filter(gs_name%in%c('Kegg_antigen_processing_and_presentation'))

srt_tumor_Bcell_noplasma=AddModuleScore(
  srt_tumor_Bcell_noplasma,
  features = list(unique(antigen_pathway$gene_symbol)),
  name='antigen_processing_and_presentation'
)
head(srt_tumor_Bcell_noplasma)

VlnPlot(srt_tumor_Bcell_noplasma, features = c('antigen_processing_and_presentation1'),group.by = 'treatment_Resp',pt.size = 0) +
  scale_fill_manual(values = mycolor) + stat_compare_means(method = "t.test")
ggsave(".../compare_antigen.pdf",width = 5,height = 6,units = "in")

#Fig S17g####
b_cell_activation=geneset_GOBP %>% filter(gs_name%in%c('Gobp_b_cell_activation'))

srt_tumor_Bcell_noplasma=AddModuleScore(
  srt_tumor_Bcell_noplasma,
  features = list(unique(b_cell_activation$gene_symbol)),
  name='b_cell_activation'
)
head(srt_tumor_Bcell_noplasma)

VlnPlot(srt_tumor_Bcell_noplasma, features = c('b_cell_activation1'),group.by = 'treatment_Resp',pt.size = 0) +
  scale_fill_manual(values = mycolor) + stat_compare_means(method = "t.test")
ggsave(".../compare_b_cell_activation.pdf",width = 5,height = 6,units = "in")






#CellChat analysis#####
head(srt_tumor)
Idents(srt_tumor) <- "group"
srt_tumor$samples=srt_tumor$patient_ID
srt_tumor$samples=as.factor(srt_tumor$samples)
levels(Idents(srt_tumor))

srt_tumor_response=srt_tumor %>% subset(treatment_Resp %in% c('anti-PD1_R'))
srt_tumor_nonresponse=srt_tumor %>% subset(treatment_Resp %in% c('anti-PD1_NR'))

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)




#Fig S16c######
#1.response####
data.input <- srt_tumor_response[["RNA"]]$data # normalized data matrix
labels <- Idents(srt_tumor_response)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellChat_R <- createCellChat(object = srt_tumor_response, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact","Secreted Signaling"), key = "annotation") # use Secreted Signaling

cellChat_R@DB <- CellChatDB.use

set.seed(1)
cellChat_R <- subsetData(cellChat_R) 
future::plan("multisession", workers = 4) # do parallel
cellChat_R <- identifyOverExpressedGenes(cellChat_R)

library(future)
options(future.globals.maxSize = 10 * 1024^3) # 10 GB
ptm = Sys.time()
cellChat_R <- identifyOverExpressedInteractions(cellChat_R)


execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
cellChat_R <- computeCommunProb(cellChat_R, type = "triMean")
cellChat_R <- filterCommunication(cellChat_R, min.cells = 10)
cellChat_R <- computeCommunProbPathway(cellChat_R)



#Calculate the aggregated cell-cell communication network

cellChat_R <- aggregateNet(cellChat_R)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellChat_R@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat_R@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat_R@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#save as 8*8


mat <- cellChat_R@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#save as 8*8




#2.nonresponse####
data.input <- srt_tumor_nonresponse[["RNA"]]$data # normalized data matrix
labels <- Idents(srt_tumor_nonresponse)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellChat_NR <- createCellChat(object = srt_tumor_nonresponse, group.by = "ident", assay = "RNA")


CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Cell-Cell Contact","Secreted Signaling"), key = "annotation") # use Secreted Signaling

cellChat_NR@DB <- CellChatDB.use

set.seed(1)
cellChat_NR <- subsetData(cellChat_NR) 
future::plan("multisession", workers = 4) # do parallel
cellChat_NR <- identifyOverExpressedGenes(cellChat_NR)

library(future)
options(future.globals.maxSize = 10 * 1024^3) # 10 GB
ptm = Sys.time()
cellChat_NR <- identifyOverExpressedInteractions(cellChat_NR)

execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
cellChat_NR <- computeCommunProb(cellChat_NR, type = "triMean")
cellChat_NR <- filterCommunication(cellChat_NR, min.cells = 10)
cellChat_NR <- computeCommunProbPathway(cellChat_NR)

cellChat_NR <- aggregateNet(cellChat_NR)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))

ptm = Sys.time()
groupSize <- as.numeric(table(cellChat_NR@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellChat_NR@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellChat_NR@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
#save as 8*8

mat <- cellChat_NR@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#save as 8*8



#Fig S16d#####
netVisual_bubble(cellChat_R,
                 sources.use = c('B'), targets.use = c('T','MonoMac','NK','Stromal','DC','ILC','Mast'), 
                 #pairLR.use = df,
                 remove.isolate = FALSE)

netVisual_bubble(cellChat_NR,
                 sources.use = c('B'), targets.use = c('T','MonoMac','NK','Stromal','DC','ILC','Mast'), 
                 #pairLR.use = df,
                 remove.isolate = FALSE)

vec <- c("PTPRC_MRC1","PPIA_BSG","HLA-F_LILRB1","HLA-E_KLRC1",   #"MIF_(CD74+CD44)","IL16_CD4",
         "HLA-E_CD94:NKG2A",       #"HLA-E_CD8B","HLA-B_CD8A",
         "CLEC2C_KLRB1","CD99_PILRA","CD55_ADGRE5")

df <- data.frame(interaction_name = vec, stringsAsFactors = FALSE)
df

netVisual_bubble(cellChat_R,
                 sources.use = c('B'), targets.use = c('T','MonoMac','NK','Stromal','DC','ILC','Mast'), 
                 pairLR.use = df,
                 remove.isolate = FALSE)
#save as 4*4

netVisual_bubble(cellChat_NR, 
                 sources.use = c('B'), targets.use = c('T','MonoMac','NK','Stromal','DC','ILC','Mast'), 
                 pairLR.use = df,
                 remove.isolate = FALSE)
#save as 4*4








###TIM-1 correlation analysis#####
#Fig R9####

srt_tumor_Bcell=srt_tumor %>% subset(subset=group %in% c('B')) #19248 cells

breg_genes <- list(c(
  "IL10",
  "TGFB1", "TGFB2",
  "IL12A", "EBI3",
  "CD24","CD38","CD27","FCGR2B","CD22","FCRL3","CD274",
  "HAVCR1","CD40LG")
)

srt_tumor_Bcell <- AddModuleScore(srt_tumor_Bcell, features = breg_genes, name = "BregScore")

meta_value <- srt_tumor_Bcell@meta.data$BregScore1

gene_expr <- FetchData(srt_tumor_Bcell, vars = "HAVCR1")[,1]

cor.test(meta_value, gene_expr, method = "pearson")

library(ggplot2)

df <- data.frame(
  meta_value = meta_value,
  gene_expr = gene_expr
)

ggplot(df, aes(x = meta_value, y = gene_expr)) +
  geom_point(size = 1.5, alpha = 0.6, color = "#1f78b4") +   
  geom_smooth(method = "lm", color = "#e31a1c", se = TRUE) + 
  theme_classic(base_size = 14) +                            
  labs(
    x = "Breg score",
    y = "HAVCR1 expression",
    title = "Correlation: Breg score ~ HAVCR1"
  ) +
  stat_cor(method = "pearson",                             
           label.x.npc = "left", 
           label.y.npc = "top",
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")))

ggsave(".../TIM1_correlation.pdf",width = 4,height = 4,units = "in")





