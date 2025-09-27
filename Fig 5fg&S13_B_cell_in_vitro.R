#read in counts data#####
df_Bcell=readRDS(".../Bcell_counts.rds")
TPM_Bcell=readRDS(".../Bcell_TPM.rds")

#Fig 5f####
pca_df=as.data.frame(t(TPM_Bcell))
pca_df[1:12,1:4]
pca_df=pca_df[,-1] 
pca_df$Sample=rownames(pca_df)
pca_df$Group[pca_df$Sample%in%c("IgG_1","IgG_2","IgG_3")]='IgG'
pca_df$Group[pca_df$Sample%in%c("aTIM1_1","aTIM1_2","aTIM1_3")]='aTIM1'
pca_df$Group[pca_df$Sample%in%c("STING_1.0_IgG_1","STING_1.0_IgG_2","STING_1.0_IgG_3")]='STING+IgG'
pca_df$Group[pca_df$Sample%in%c("STING_1.0_aTIM1_1","STING_1.0_aTIM1_2","STING_1.0_aTIM1_3")]='STING+aTIM1'

pca_df[,colSums(pca_df[,1:(dim(pca_df)[2]-2)])==0]=NULL 
com1=prcomp(pca_df[,1:(dim(pca_df)[2]-2)], center=T, scale.=T)

df1<-com1$x
head(df1)

df1<-data.frame(df1,pca_df$Sample, pca_df$Group)
colnames(df1)[c(dim(df1)[1]+1,dim(df1)[1]+2)]=c('Sample','Group')
head(df1)

summ<-summary(com1)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")

set.seed(1)
library(vegan)
bray=adonis2(pca_df[,1:(dim(pca_df)[2]-2)]~Group, data=pca_df,permutations = 999,method = 'bray')
adonis=paste0('adonis R^2: ',round(bray$R2,3),'; P-value: ',bray$`Pr(>F)`)

library(ggord)
library(ggsci)
library(ggplot2)
library(pairwiseAdonis)
pca_df$Group <- factor(pca_df$Group, levels = c('IgG','aTIM1',
                                                'STING+IgG','STING+aTIM1'))

#DESeq2 PCA
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(stringr)

colnames(df_Bcell)
df_Bcell[1:3,1:4]
df_Bcell=df_Bcell[rowMeans(df_Bcell[,1:12])>0,] #20626
df_Bcell=df_Bcell[-1,]

group_list <- c('IgG','IgG','IgG',
                'aTIM1','aTIM1','aTIM1',
                'STING+IgG','STING+IgG','STING+IgG',
                'STING+aTIM1','STING+aTIM1','STING+aTIM1')
group_list <- factor(group_list,levels = c('IgG','aTIM1',
                                           'STING+IgG','STING+aTIM1'))
                     
condition = group_list
coldata <- data.frame(row.names = colnames(df_Bcell), condition) 
dds <- DESeqDataSetFromMatrix(countData = df_Bcell,
                              colData = coldata,
                              design = ~condition)  
dds$condition<- relevel(dds$condition, ref = "IgG") 
dds <- DESeq(dds)


#DESeq2 PCA analysis
vsd <- vst(dds, nsub=1000)
plotPCA(vsd)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)

pcaData$Antibodies[pcaData$group%in%c("IgG","STING+IgG")]='IgG' 
pcaData$Antibodies[pcaData$group%in%c("aTIM1","STING+aTIM1")]='aTIM1' 
pcaData$Antibodies=factor(pcaData$Antibodies, levels = c('IgG','aTIM1'))

pcaData$Treatment[pcaData$group%in%c("IgG","aTIM1")]='PBS'
pcaData$Treatment[pcaData$group%in%c("STING+aTIM1","STING+IgG")]='STING'
pcaData$Treatment=factor(pcaData$Treatment, levels = c('PBS','STING'))

percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 1)

ggplot(pcaData, aes(PC1, PC2, color=Treatment,shape=Antibodies)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c("#1A6291","#C66500","#D41919"))+
  theme_bw()+
  theme(
    aspect.ratio = 1
  )
ggsave(".../PCA_4_groups.pdf",width = 4,height = 3,units = "in")


#prepare the database#####
options(timeout = 10000000)
database=msigdbr(species = "Mus musculus")
collections=msigdbr_collections()
unique(database$gs_cat)
unique(database$gs_subcat)

geneset_GOBP=msigdbr(species = 'Mus musculus',
                     category = 'C5',
                     subcategory = 'GO:BP') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_GOBP$gs_name=str_to_title(geneset_GOBP$gs_name)
head(geneset_GOBP)

geneset_reactome=msigdbr(species = 'Mus musculus',
                         category = 'C2',
                         subcategory = 'CP:REACTOME') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_reactome$gs_name=str_to_title(geneset_reactome$gs_name)
head(geneset_reactome)

geneset_KEGG=msigdbr(species = 'Mus musculus',
                     category = 'C2',
                     subcategory = 'CP:KEGG') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_KEGG$gs_name=str_to_title(geneset_KEGG$gs_name)
head(geneset_KEGG)

geneset_hallmark=msigdbr(species = 'Mus musculus',
                         category = 'H') %>% dplyr::select(gs_name,gene_symbol) #entrez_gene,gene_symbol
geneset_hallmark$gs_name=str_to_title(geneset_hallmark$gs_name)
head(geneset_hallmark)


#Fig 5g####
res=results(dds, contrast = c('condition','STING+IgG','IgG'))

nrDEG_DESeq2 <- as.data.frame(res)
rld <- rlog(dds)
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),]

plotdf=nrDEG_DESeq2
plotdf=plotdf %>% dplyr::arrange(dplyr::desc(log2FoldChange))

ranks=plotdf$log2FoldChange
names(ranks)=rownames(plotdf)
head(ranks)
ranks=sort(ranks,decreasing = T)


#####GSEA GOBP
set.seed(1)
Res_GSEA_GOBP=GSEA(ranks, TERM2GENE = geneset_GOBP,seed = T)

#####GSEA REACTOME
set.seed(1)
Res_GSEA_REACTOME=GSEA(ranks, TERM2GENE = geneset_reactome,seed = T)

#####GSEA KEGG
set.seed(1)
Res_GSEA_KEGG=GSEA(ranks, TERM2GENE = geneset_KEGG,seed = T)


#####GSEA Hallmark
set.seed(1)
Res_GSEA_HALLMARK=GSEA(ranks, TERM2GENE = geneset_hallmark,seed = T)



library(GSVA)
#database
library(msigdbr)
library(stringr)
library(dplyr)
geneset_all=geneset_GOBP %>% bind_rows(geneset_reactome) %>% bind_rows(geneset_KEGG) %>% bind_rows(geneset_hallmark)

# Extract gene set names that contain 'b_cell'
b_cell_gene_sets <- grep("b_cell", geneset_all$gs_name, value = TRUE)
head(b_cell_gene_sets)

# To retrieve the corresponding gene sets
b_cell_gene_sets_content <- geneset_all %>% filter(gs_name %in% unique(b_cell_gene_sets))

gs <- split(b_cell_gene_sets_content$gene_symbol,b_cell_gene_sets_content$gs_name)

gsvaPar <- ssgseaParam(as.matrix(df_Bcell), gs)

gsva_es <- gsva(gsvaPar, verbose=FALSE)

colnames(gsva_es)
colnames(gsva_es)=c("IgG_1","IgG_2","IgG_3",
                    "aTIM1_1","aTIM1_2","aTIM1_3",
                    "STING+IgG_1","STING+IgG_2","STING+IgG_3",
                    "STING+aTIM1_1","STING+aTIM1_2","STING+aTIM1_3")

library(pheatmap)
pheatmap(gsva_es,
         color = colorRampPalette(c("#0F7B9F", "white", "#D83215"))(50),
         cellwidth = 15,cellheight = 9,
         cluster_cols = F, cluster_rows = T,
         border_color = 'white',
         scale = 'row',
         fontsize = 9,
         #treeheight_row = 10,
         #clustering_method = "complete",
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         angle_col = '45',
         #cutree_rows = 6,
         #main='control vs com'
         gaps_col = c(3,6,9),
         gaps_row = c(35,36))
#width 10, height 7






#Fig S13####
res=results(dds, contrast = c('condition','STING+aTIM1','STING+IgG'))

nrDEG_DESeq2 <- as.data.frame(res)
rld <- rlog(dds)
normal_gset <- assay(rld) 
nrDEG_DESeq2 = nrDEG_DESeq2[order(nrDEG_DESeq2$log2FoldChange),]

plotdf=nrDEG_DESeq2
plotdf=plotdf %>% dplyr::arrange(dplyr::desc(log2FoldChange))

ranks=plotdf$log2FoldChange
names(ranks)=rownames(plotdf)
head(ranks)
ranks=sort(ranks,decreasing = T)



#####GSEA GOBP
set.seed(1)
Res_GSEA_GOBP=GSEA(ranks, TERM2GENE = geneset_GOBP,seed = T,pvalueCutoff = 1)

#####GSEA REACTOME
set.seed(1)
Res_GSEA_REACTOME=GSEA(ranks, TERM2GENE = geneset_reactome,seed = T,pvalueCutoff = 1)

#####GSEA KEGG
set.seed(1)
Res_GSEA_KEGG=GSEA(ranks, TERM2GENE = geneset_KEGG,seed = T,pvalueCutoff = 1)

#####GSEA Hallmark
set.seed(1)
Res_GSEA_HALLMARK=GSEA(ranks, TERM2GENE = geneset_hallmark,seed = T,pvalueCutoff = 1)

library(GSVA)

#database
library(msigdbr)
library(stringr)
library(dplyr)
geneset_all=geneset_GOBP %>% bind_rows(geneset_reactome) %>% bind_rows(geneset_KEGG) %>% bind_rows(geneset_hallmark)

b_cell_gene_sets <- grep("antigen_processing_and_presentation", geneset_all$gs_name, value = TRUE)

# To retrieve the corresponding gene sets
b_cell_gene_sets_content <- geneset_all %>% filter(gs_name %in% unique(b_cell_gene_sets))

gs <- split(b_cell_gene_sets_content$gene_symbol,b_cell_gene_sets_content$gs_name)

gsvaPar <- ssgseaParam(as.matrix(df_Bcell), gs)

gsva_es <- gsva(gsvaPar, verbose=FALSE)

colnames(gsva_es)
colnames(gsva_es)=c("IgG_1","IgG_2","IgG_3",
                    "aTIM1_1","aTIM1_2","aTIM1_3",
                    "STING+IgG_1","STING+IgG_2","STING+IgG_3",
                    "STING+aTIM1_1","STING+aTIM1_2","STING+aTIM1_3")

library(pheatmap)
pheatmap(gsva_es, 
         color = colorRampPalette(c("#0F7B9F", "white", "#D83215"))(50),
         cellwidth = 15,cellheight = 9,
         cluster_cols = F, cluster_rows = T,
         border_color = 'white',
         scale = 'row', 
         fontsize = 9,
         #treeheight_row = 10,
         #clustering_method = "complete",
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         angle_col = '45',
         #cutree_rows = 6,
         #main='control vs com'
         gaps_col = c(3,6,9),
         gaps_row = c(35,36))
#width 11, height 7






