set.seed(1)
final_df=readRDS(".../merged_counts.rds")
final_TPM=readRDS(".../merged_TPM.rds")

#DESeq2
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(stringr)

final_df=final_df[rowMeans(final_df[,1:63])>0,]

group_list <- c('Copa_Ctrl','Copa_Ctrl','Copa_Ctrl',
                'Copa_PD1','Copa_PD1','Copa_PD1','Copa_PD1',
                'Nano_Ctrl','Nano_Ctrl','Nano_Ctrl',
                'Nano_PD1','Nano_PD1','Nano_PD1',
                'Sting_RIL175_Ctrl','Sting_RIL175_Ctrl','Sting_RIL175_Ctrl',
                'Sting_RIL175_STING','Sting_RIL175_STING','Sting_RIL175_STING',
                'Sting_RIL175_DP','Sting_RIL175_DP','Sting_RIL175_DP',
                'Sting_RIL175_SDP','Sting_RIL175_SDP','Sting_RIL175_SDP',
                'Sting_HCA1_Ctrl','Sting_HCA1_Ctrl','Sting_HCA1_Ctrl',
                'Sting_HCA1_STING','Sting_HCA1_STING','Sting_HCA1_STING',
                'JITC_Ctrl','JITC_Ctrl','JITC_Ctrl',
                'JITC_PD1','JITC_PD1','JITC_PD1',
                'Hepatology_Ctrl','Hepatology_Ctrl','Hepatology_Ctrl','Hepatology_Ctrl',
                'Hepatology_DP','Hepatology_DP','Hepatology_DP','Hepatology_DP','Hepatology_DP',
                'Hepatology_PD1','Hepatology_PD1','Hepatology_PD1','Hepatology_PD1',
                'CXCR4_Ctrl','CXCR4_Ctrl','CXCR4_Ctrl','CXCR4_Ctrl',
                'CXCR4_PD1','CXCR4_PD1','CXCR4_PD1',
                'GC_ICB_Ctrl','GC_ICB_Ctrl','GC_ICB_Ctrl',
                'GC_ICB_PD1_CTLA4','GC_ICB_PD1_CTLA4','GC_ICB_PD1_CTLA4')


condition = group_list
coldata <- data.frame(row.names = colnames(final_df), condition) 
dds <- DESeqDataSetFromMatrix(countData = final_df,
                              colData = coldata,
                              design = ~condition)  

#DESeq2 PCA analysis (Fig 1a,1b)#####
library(ggplot2)
library(ggsci)
vsd <- vst(dds, nsub=1000)
plotPCA(vsd)

pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 1)

pcaData$Projects[pcaData$group%in%c("Copa_Ctrl","Copa_PD1")]='HCA-1 cohort'
pcaData$Projects[pcaData$group%in%c("Nano_Ctrl","Nano_PD1")]='Xiao et al., 2022'
pcaData$Projects[pcaData$group%in%c("Sting_RIL175_Ctrl","Sting_RIL175_STING",'Sting_RIL175_DP','Sting_RIL175_SDP')]='RIL-175 cohort'
pcaData$Projects[pcaData$group%in%c("Sting_HCA1_Ctrl","Sting_HCA1_STING")]='HCA-1 cohort'
pcaData$Projects[pcaData$group%in%c("JITC_Ctrl","JITC_PD1")]='Shigeta et al., 2020b'
pcaData$Projects[pcaData$group%in%c("Hepatology_Ctrl","Hepatology_DP",'Hepatology_PD1')]='Shigeta et al., 2020a'
pcaData$Projects[pcaData$group%in%c("CXCR4_Ctrl","CXCR4_PD1")]='Morita et al., 2024'
pcaData$Projects[pcaData$group%in%c("GC_ICB_Ctrl","GC_ICB_PD1_CTLA4")]='Chen et al., 2024'


pcaData$Treatment[pcaData$group%in%c("Copa_Ctrl","Nano_Ctrl",'Sting_RIL175_Ctrl','Sting_HCA1_Ctrl','JITC_Ctrl','Hepatology_Ctrl',
                                     'CXCR4_Ctrl','GC_ICB_Ctrl')]='Control'
pcaData$Treatment[pcaData$group%in%c("Copa_PD1","Nano_PD1",'Sting_RIL175_STING','Sting_RIL175_DP','Sting_RIL175_SDP',
                                     'Sting_HCA1_STING','JITC_PD1','Hepatology_DP','Hepatology_PD1','CXCR4_PD1','GC_ICB_PD1_CTLA4')]='Immunoactivator'


#Correct the batch
library(sva)

#use RNA-Seq count data
batch <- pcaData$Projects
cov1 <- as.numeric(factor(pcaData$Treatment))
adjusted_counts <- ComBat_seq(as.matrix(final_df), 
                              batch=batch,
                              group=cov1)
adjusted_counts=as.data.frame(adjusted_counts)

#After removing batch effect, PCA
library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(stringr)

colnames(adjusted_counts)
adjusted_counts[1:3,1:4]
adjusted_counts=adjusted_counts[rowMeans(adjusted_counts[,1:63])>0,] 

group_list <- c('Copa_Ctrl','Copa_Ctrl','Copa_Ctrl',
                'Copa_PD1','Copa_PD1','Copa_PD1','Copa_PD1',
                'Nano_Ctrl','Nano_Ctrl','Nano_Ctrl',
                'Nano_PD1','Nano_PD1','Nano_PD1',
                'Sting_RIL175_Ctrl','Sting_RIL175_Ctrl','Sting_RIL175_Ctrl',
                'Sting_RIL175_STING','Sting_RIL175_STING','Sting_RIL175_STING',
                'Sting_RIL175_DP','Sting_RIL175_DP','Sting_RIL175_DP',
                'Sting_RIL175_SDP','Sting_RIL175_SDP','Sting_RIL175_SDP',
                'Sting_HCA1_Ctrl','Sting_HCA1_Ctrl','Sting_HCA1_Ctrl',
                'Sting_HCA1_STING','Sting_HCA1_STING','Sting_HCA1_STING',
                'JITC_Ctrl','JITC_Ctrl','JITC_Ctrl',
                'JITC_PD1','JITC_PD1','JITC_PD1',
                'Hepatology_Ctrl','Hepatology_Ctrl','Hepatology_Ctrl','Hepatology_Ctrl',
                'Hepatology_DP','Hepatology_DP','Hepatology_DP','Hepatology_DP','Hepatology_DP',
                'Hepatology_PD1','Hepatology_PD1','Hepatology_PD1','Hepatology_PD1',
                'CXCR4_Ctrl','CXCR4_Ctrl','CXCR4_Ctrl','CXCR4_Ctrl',
                'CXCR4_PD1','CXCR4_PD1','CXCR4_PD1',
                'GC_ICB_Ctrl','GC_ICB_Ctrl','GC_ICB_Ctrl',
                'GC_ICB_PD1_CTLA4','GC_ICB_PD1_CTLA4','GC_ICB_PD1_CTLA4')

condition = group_list
coldata <- data.frame(row.names = colnames(adjusted_counts), condition) 
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                              colData = coldata,
                              design = ~condition)  

library(ggplot2)
library(ggsci)
vsd <- vst(dds, nsub=1000)
plotPCA(vsd)
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"), digits = 1)

pcaData$Projects[pcaData$group%in%c("Copa_Ctrl","Copa_PD1")]='HCA-1 cohort'
pcaData$Projects[pcaData$group%in%c("Nano_Ctrl","Nano_PD1")]='Xiao et al., 2022'
pcaData$Projects[pcaData$group%in%c("Sting_RIL175_Ctrl","Sting_RIL175_STING",'Sting_RIL175_DP','Sting_RIL175_SDP')]='RIL-175 cohort'
pcaData$Projects[pcaData$group%in%c("Sting_HCA1_Ctrl","Sting_HCA1_STING")]='HCA-1 cohort'
pcaData$Projects[pcaData$group%in%c("JITC_Ctrl","JITC_PD1")]='Shigeta et al., 2020b'
pcaData$Projects[pcaData$group%in%c("Hepatology_Ctrl","Hepatology_DP",'Hepatology_PD1')]='Shigeta et al., 2020a'
pcaData$Projects[pcaData$group%in%c("CXCR4_Ctrl","CXCR4_PD1")]='Morita et al., 2024'
pcaData$Projects[pcaData$group%in%c("GC_ICB_Ctrl","GC_ICB_PD1_CTLA4")]='Chen et al., 2024'

pcaData$Projects <- factor(pcaData$Projects,
                        levels = c("HCA-1 cohort", "RIL-175 cohort",'Shigeta et al., 2020a','Shigeta et al., 2020b',
                                   'Xiao et al., 2022','Morita et al., 2024','Chen et al., 2024'))


ggplot(pcaData, aes(PC1, PC2, color=Projects)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  scale_color_manual(values = c('#46A5DA',
                                '#FAEC95',
                                '#5080CC',
                                '#EB4C20',
                                '#DEA077',
                                '#BF76CC',
                                '#9FCCA2'))+
  theme_bw()+
  theme(
    aspect.ratio = 1
  )
ggsave(".../PCA_based_on_DESeq2_by_project_after_combat_new.pdf",width = 7,height = 5.5,units = "in")



#draw by treatment
pcaData$Treatment[pcaData$group%in%c("Copa_Ctrl","Nano_Ctrl",'Sting_RIL175_Ctrl','Sting_HCA1_Ctrl','JITC_Ctrl','Hepatology_Ctrl',
                                     'CXCR4_Ctrl','GC_ICB_Ctrl')]='Control'
pcaData$Treatment[pcaData$group%in%c("Copa_PD1","Nano_PD1",'Sting_RIL175_STING','Sting_RIL175_DP','Sting_RIL175_SDP',
                                     'Sting_HCA1_STING','JITC_PD1','Hepatology_DP','Hepatology_PD1','CXCR4_PD1','GC_ICB_PD1_CTLA4')]='Immunoactivator'

t.test(PC1 ~ Treatment, data = pcaData)
t.test(PC2 ~ Treatment, data = pcaData)

ggplot(pcaData, aes(PC1, PC2, color=Treatment)) + 
  geom_point(size=3) +
  stat_ellipse(type = "t", level = 0.95, linetype = "dashed") +  
  xlab(paste0("PC1: ",percentVar[1],"% variance (P=0.003)")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance (P=0.075)")) +
  scale_color_igv()+
  theme_bw()+
  theme(
    aspect.ratio = 1
  )
ggsave(".../PCA_based_on_DESeq2_by_treatment_after_combat_new_with_P_value.pdf",width = 7,height = 5.5,units = "in")






#after batch removal, convert to TPM####
#read in the gene_name and gene_id conversion table
ensembl_file=read.table(".../mart_export_Ensembl_Genes_107_Mouse_genes_GRCm39.txt",sep='\t',header=T,check.names = T)
ensembl_file = ensembl_file[!duplicated(ensembl_file$Gene.stable.ID),] 
colnames(ensembl_file)[1]='Geneid'
ensembl_file = ensembl_file %>% filter(Gene.name%in%rownames(adjusted_counts))
ensembl_file <- ensembl_file[!duplicated(ensembl_file$Gene.name),] #38620 genes

library(dplyr)
#merge data and conversion table
adjusted_counts$Gene.name=rownames(adjusted_counts)
merged_adj_count=left_join(adjusted_counts, ensembl_file, by='Gene.name')

colnames(merged_adj_count)
gene_info=merged_adj_count[,c(64,66,67)]
gene_info$length = gene_info[,3] - gene_info[,2]

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

adj_TPM <- apply(merged_adj_count[,1:63], 2, function(x) tpm(x, gene_info$length))
dim(adj_TPM) 
colSums(adj_TPM)
adj_TPM=as.data.frame(adj_TPM)
adj_TPM$gene_name=merged_adj_count$Gene.name

length(unique(adj_TPM$gene_name)) #38620

rownames(adj_TPM)=adj_TPM$gene_name
class(adj_TPM)
adj_TPM=adj_TPM[,-64]
adj_TPM=adj_TPM[rowMeans(adj_TPM[,1:63])>0,]  #38620
colSums(adj_TPM[,1:63])

#Fig S1b, markers for B cells#######
colnames(adj_TPM)
compare_TPM=adj_TPM[,c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60,
                       4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)]
colnames(compare_TPM)
compare_TPM=compare_TPM[,c(1:3,10:12,7:9,16:19,13:15,4:6,20:26,
                           27:30,43:45,34:42,54:57,49:53,46:48,31:33,58:63)]
compare_Bcell_heatmap=compare_TPM[,c(27:30,31:33,43:46,52:54,55:57,58:60,34:36)]
colnames(compare_Bcell_heatmap)
library(pheatmap)
pheatmap(as.matrix(compare_Bcell_heatmap[c('Cd19','Ms4a1','Cd79a','Cd79b',
                                           'Cd27','Cd38','Cd69','Cd80','Cd86','Icosl',
                                           'Havcr1','Tigit','Il10','Il12a','Ebi3','Cd274','Pdcd1lg2','Tnfsf18',
                                           'Aicda','Bcl6','Fas',
                                           'Fcrl5', 'Zbtb20','Ighg1','Igha','Ighm'
                                           ),]), 
         cellwidth = 12,cellheight = 12,
         cluster_cols = F, cluster_rows = F,
         border_color = 'white',
         scale = 'row', 
         fontsize = 8,
         angle_col = '45',
         gaps_col = c(4,7,20),
         gaps_row = c(4,10,18,21))

#save #width 7, height 9
#".../B_cell_marker_compare.pdf"





#merged_human is TPM data with human genes#####
#convert mouse gene to human
library(biomaRt)
human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org")
mouse.gene=rownames(adj_TPM)
m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                values = mouse.gene,mart = mouse,
                attributesL = c("hgnc_symbol"),
                martL = human,uniqueRows = T)

adj_TPM$MGI.symbol=rownames(adj_TPM)
library(dplyr)
adj_TPM=left_join(adj_TPM,m2h.g,by='MGI.symbol')
TPM_human=na.omit(adj_TPM) #20476 genes

merged_human <- TPM_human %>%
  group_by(HGNC.symbol) %>%
  summarise(across(1:63, ~ sum(.x, na.rm = TRUE)))
colSums(merged_human[,2:64])

merged_human=as.data.frame(merged_human)
rownames(merged_human)=merged_human$HGNC.symbol
merged_human=merged_human[,-1]


#xCell####
library(xCell) 
xCell.data=xCellAnalysis(merged_human)
xCell_res=as.data.frame(xCell.data)

colnames(xCell_res)
colSums(xCell_res)

rownames(xCell_res)
xCell_res_Bcell=xCell_res[c('B-cells','Class-switched memory B-cells',
                            'Memory B-cells','naive B-cells','Plasma cells',
                            'pro B-cells'),]

res_xCell_Bcell=xCell_res_Bcell%>% as.data.frame()
res_xCell_Bcell$method='xCell'
res_xCell_Bcell$subtype=rownames(res_xCell_Bcell)

colnames(res_xCell_Bcell)[c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60)]
colnames(res_xCell_Bcell)[c(4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)]

res_xCell_Bcell=res_xCell_Bcell[,c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60,
                 4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)]

colnames(res_xCell_Bcell)
res_xCell_Bcell=res_xCell_Bcell[,c(1:3,10:12,7:9,16:19,13:15,4:6,20:26,
                                   27:30,43:45,34:42,54:57,49:53,46:48,31:33,58:63)]
res_xCell_Bcell$method='xCell'
res_xCell_Bcell$subtype=rownames(res_xCell_Bcell)

#Fig S1a, xCell: compare PD1 vs STING######
compare_Bcell=res_xCell_Bcell[,c(27:30,31:33,43:46,52:54,55:57,58:60,34:36)]
colnames(compare_Bcell)

library(tidyverse)

df_compare_Bcell <- compare_Bcell %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = -Sample,
    names_to = "CellType",
    values_to = "Abundance"
  )

df_compare_Bcell <- df_compare_Bcell %>%
  mutate(Group = case_when(
    grepl("^Copa_PD1", Sample) ~ "HCA1_PD1",
    grepl("^Sting_HCA1", Sample) ~ "HCA1_STING",
    grepl("^Hepatology", Sample) ~ "RIL175_PD1",
    grepl("^JITC", Sample) ~ "RIL175_PD1",
    grepl("^Nano", Sample) ~ "RIL175_PD1",
    grepl("^CXCR4", Sample) ~ "RIL175_PD1",
    grepl("^Sting_RIL175", Sample) ~ "RIL175_STING",
    TRUE ~ "Other"
  ))

table(df_compare_Bcell$Group)

library(ggplot2)
library(ggpubr)   
library(ggsci)    

my_comparisons <- list(
  c("HCA1_PD1", "HCA1_STING"),
  c("RIL175_PD1", "RIL175_STING")
)

ggplot(df_compare_Bcell, aes(x = Group, y = Abundance)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA, alpha = 0.5, color = "black") +
  geom_jitter(aes(color = Group), width = 0.2, size = 2) +
  scale_color_igv() +
  facet_wrap(~ CellType, scales = "free_y",strip.position = "top") +
  stat_compare_means(
    method = "wilcox.test",     
    label = "p.format",         
    label.y.npc = "top",        
    comparisons = my_comparisons,         
    size = 3
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 10, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 10, color = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),  
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10)
  )+
  labs(
    title = "Enrichment score of each B cell subtype across treatment groups",
    x = "", y = "Enrichment score (xCell)"
  )

ggsave(".../bulk compare/Enrichment_score_compare.pdf",width = 7,height = 7,units = "in")







#CIBERSORT######

library(immunedeconv)
#download files for CIBERSORT
set_cibersort_binary(".../CIBERSORT/CIBERSORT (1).R")
set_cibersort_mat(".../CIBERSORT/LM22.txt")

res_cibersort=deconvolute_cibersort(as.matrix(merged_human),FALSE,absolute = F,abs_method = "sig.score")

rownames(res_cibersort)
res_cibersort_Bcell=res_cibersort[c('B cells naive','B cells memory',
                                    'Plasma cells'),] %>% as.data.frame()
res_cibersort_Bcell$method='CIBERSORT'
res_cibersort_Bcell$subtype=rownames(res_cibersort_Bcell)


#methdos in immunedeconv package
human_method=deconvolution_methods
mouse_method=deconvolution_methods_mouse

#CIBERSORT Abs####
res_cibersort_abs=deconvolute_cibersort(as.matrix(merged_human),FALSE,absolute = TRUE,abs_method = "sig.score")
res_cibersort_abs_Bcell=res_cibersort_abs[c('B cells naive','B cells memory',
                                            'Plasma cells'),] %>% as.data.frame()
res_cibersort_abs_Bcell$method='CIBERSORT_Abs'
res_cibersort_abs_Bcell$subtype=rownames(res_cibersort_abs_Bcell)

#MCPcounter#####
res_MCPcounter = deconvolute_mcp_counter(as.matrix(merged_human))
res_MCPcounter_Bcell=res_MCPcounter['B lineage',]
res_MCPcounter_Bcell=as.data.frame(t(res_MCPcounter_Bcell)) 
res_MCPcounter_Bcell$method='MCPcounter'
res_MCPcounter_Bcell$subtype='B lineage'

#EPIC#####
res_EPIC = deconvolute_epic(as.matrix(merged_human),tumor=TRUE,scale_mrna=TRUE)
res_EPIC_Bcell=res_EPIC['Bcells',]
res_EPIC_Bcell=as.data.frame(t(res_EPIC_Bcell)) 
res_EPIC_Bcell$method='EPIC'
res_EPIC_Bcell$subtype='Bcells'

#quanTIseq####
res_quanTIseq = deconvolute_quantiseq(as.matrix(merged_human),tumor=TRUE,arrays = F,scale_mrna=TRUE)
res_quanTIseq_Bcell=res_quanTIseq['B.cells',]
res_quanTIseq_Bcell=as.data.frame(t(res_quanTIseq_Bcell)) 
res_quanTIseq_Bcell$method='quanTIseq'
res_quanTIseq_Bcell$subtype='B.cells'

#TIMER####
res_TIMER = deconvolute_timer(as.matrix(merged_human),indications = c(rep('lihc',57),rep('chol',6)))
res_TIMER_Bcell=res_TIMER['B_cell',]
res_TIMER_Bcell=as.data.frame(t(res_TIMER_Bcell)) 
res_TIMER_Bcell$method='TIMER'
res_TIMER_Bcell$subtype='B_cell'

#ConsensusTME####
res_ConsensusTME = deconvolute_consensus_tme(as.matrix(merged_human),
                                             indications = c(rep('lihc',57),rep('chol',6)),
                                             method = 'gsva')
res_ConsensusTME_Bcell=res_ConsensusTME[c('B_cells','Plasma_cells'),]%>% as.data.frame()
res_ConsensusTME_Bcell$method='ConsensusTME'
res_ConsensusTME_Bcell$subtype=rownames(res_ConsensusTME_Bcell)

#ABIS####
res_ABIS = deconvolute_abis(as.matrix(merged_human))
res_ABIS_Bcell=res_ABIS[c('B Naive','B Memory','Plasmablasts'),]%>% as.data.frame()
res_ABIS_Bcell$method='ABIS'
res_ABIS_Bcell$subtype=rownames(res_ABIS_Bcell)


#below is for mouse gene####
#mMCPcounter####
res_mMCPcounter = deconvolute_mmcp_counter(as.matrix(final_TPM),log2=F)
res_mMCPcounter_Bcell=res_mMCPcounter[c('B derived','Memory B cells'),]%>% as.data.frame()
res_mMCPcounter_Bcell$method='mMCPcounter'
res_mMCPcounter_Bcell$subtype=rownames(res_mMCPcounter_Bcell)

#seqImmuCC####
#count data, use raw counts
res_seqImmuCC_SVR = deconvolute_seqimmucc(as.matrix(final_df),algorithm = c("SVR"))
res_seqImmuCC_SVR_Bcell=res_seqImmuCC_SVR['B Cells',]
res_seqImmuCC_SVR_Bcell=as.data.frame(t(res_seqImmuCC_SVR_Bcell)) 
res_seqImmuCC_SVR_Bcell$method='seqImmuCC_SVR'
res_seqImmuCC_SVR_Bcell$subtype='B Cells'

res_seqImmuCC_LLSR = deconvolute_seqimmucc(as.matrix(final_df),algorithm = c("LLSR"))
res_seqImmuCC_LLSR_Bcell=res_seqImmuCC_LLSR['B Cells',]
res_seqImmuCC_LLSR_Bcell=as.data.frame(t(res_seqImmuCC_LLSR_Bcell)) 
res_seqImmuCC_LLSR_Bcell$method='seqImmuCC_LLSR'
res_seqImmuCC_LLSR_Bcell$subtype='B Cells'

#DCQ####
res_DCQ = deconvolute_dcq(as.matrix(final_TPM))
res_DCQ_Bcell=res_DCQ[c('B_cells','Dev_B_cells'),]%>% as.data.frame()
res_DCQ_Bcell$method='DCQ'
res_DCQ_Bcell$subtype=rownames(res_DCQ_Bcell)

#BASE####
res_BASE = deconvolute_base_algorithm(as.matrix(final_TPM))
res_BASE_Bcell=res_BASE[c('B_cells','Dev_B_cells'),]%>% as.data.frame()
res_BASE_Bcell$method='BASE'
res_BASE_Bcell$subtype=rownames(res_BASE_Bcell)



#merge results from different methods#####
plotdf=res_MCPcounter_Bcell %>% rbind(res_EPIC_Bcell) %>% rbind(res_quanTIseq_Bcell) %>% 
  rbind(res_xCell_Bcell) %>% rbind(res_cibersort_Bcell) %>% rbind(res_cibersort_abs_Bcell) %>% 
  rbind(res_TIMER_Bcell) %>% rbind(res_ConsensusTME_Bcell) %>% rbind(res_ABIS_Bcell) %>%
  rbind(res_mMCPcounter_Bcell) %>% rbind(res_seqImmuCC_SVR_Bcell) %>% rbind(res_seqImmuCC_LLSR_Bcell) %>% 
  rbind(res_DCQ_Bcell) %>% rbind(res_BASE_Bcell)

plotdf$name=paste0(plotdf$method,'_',plotdf$subtype)
rownames(plotdf)=plotdf$name


colnames(plotdf)[c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60)]
colnames(plotdf)[c(4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)]


plotdf=plotdf[,c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60,
                 4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)]
plotdf=plotdf[,c(1:3,10:12,7:9,16:19,13:15,4:6,20:26,
                 27:30,43:45,34:42,54:57,49:53,46:48,31:33,58:63)]



rownames(plotdf)
rownames(plotdf)[c(7,11,12,14,15,18,20,21,25,27,28,29)]
my_colors <- colorRampPalette(c("blue", "white", "red"))(100)


plotdf_ht=plotdf[-c(7,11,12,14,15,18,20,21,25,27,28,29),]
plotdf_ht=plotdf_ht[order(rownames(plotdf_ht)), ]

pheatmap(as.matrix(plotdf_ht), 
         cellwidth = 8,cellheight = 12,
         cluster_cols = F, cluster_rows = F,
         border_color = 'white',
         scale = 'row', 
         fontsize = 5,
         color = my_colors,
         #treeheight_row = 10,
         #clustering_method = "complete",
         #clustering_distance_rows = "euclidean",
         #clustering_distance_cols = "euclidean",
         angle_col = '45',
         #cutree_rows = 6,
         #main='control vs com'
         gaps_col = c(26),
         gaps_row = c())


library(ComplexHeatmap)

colnames(plotdf_ht)

#Hepatology: Shigeta et al., 2020a
#JITC: Shigeta et al., 2020b

# notes
Group_vector <- factor(c(rep("Control", 26), rep("Immunoactivator", 37)), levels = c("Control", "Immunoactivator"))

Cohort_vector <- factor(c(rep('HCA-1 cohort',6),
                          rep('RIL-175 cohort',3),
                          rep('Shigeta et al., 2020a',4),
                          rep('Shigeta et al., 2020b',3),
                          rep('Xiao et al., 2022',3),
                          rep('Morita et al., 2024',4),
                          rep('Chen et al., 2024',3),
                          rep('HCA-1 cohort',7),
                          rep('RIL-175 cohort',9),
                          rep('Shigeta et al., 2020a',9),
                          rep('Shigeta et al., 2020b',3),
                          rep('Xiao et al., 2022',3),
                          rep('Morita et al., 2024',3),
                          rep('Chen et al., 2024',3)),
                        levels = c("HCA-1 cohort", "RIL-175 cohort",'Shigeta et al., 2020a','Shigeta et al., 2020b',
                                   'Xiao et al., 2022','Morita et al., 2024','Chen et al., 2024'))

Treament_vector <- factor(c(rep('IgG/PBS',26),
                            rep('aPD1',4),
                            rep('STINGa',6),
                            rep('aPD1+aVEGFR2',3),
                            rep('STINGa+aPD1+aVEGFR2',3),
                            rep('aPD1',4),
                            rep('aPD1+aVEGFR2',5),
                            rep('aPD1',9),
                            rep('aPD1+aCTLA4',3)),
                          levels = c('IgG/PBS','aPD1','STINGa','aPD1+aVEGFR2','STINGa+aPD1+aVEGFR2','aPD1+aCTLA4'))

col_annotation <- HeatmapAnnotation(
  Group = Group_vector,  
  Cohort = Cohort_vector,
  Treament = Treament_vector,
  col = list(Group = c("Control" = "lightblue", "Immunoactivator" = "pink"), 
             Cohort= c('HCA-1 cohort'='#46A5DA',
                       'RIL-175 cohort'='#FAEC95',
                       'Shigeta et al., 2020a'='#5080CC',
                       'Shigeta et al., 2020b'='#EB4C20',
                       'Xiao et al., 2022'='#DEA077',
                       'Morita et al., 2024'='#BF76CC',
                       'Chen et al., 2024'='#9FCCA2'),
             Treament=c('IgG/PBS'='#40BFCC',
                        'aPD1'='#CC6F62',
                        'STINGa'='#F5BC1D',
                        'aPD1+aVEGFR2'='#C19FCC',
                        'STINGa+aPD1+aVEGFR2'='#F22727',
                        'aPD1+aCTLA4'='#3975AC')) 
)


library(ggsci)
library(scales)

colors <- pal_d3("category20")(20)
show_col(colors)

row_annotation <- rowAnnotation(
  Method = c('ABIS','CIBERSORT','CIBERSORT_Abs','ConsensusTME','DCQ','EPIC','MCPcounter',
             'mMCPcounter','mMCPcounter','quanTIseq','seqImmuCC_SVR','TIMER',rep('xCell',5)),  
  col = list(Method = c("ABIS" = '#1F77B4', 
                        "CIBERSORT" = '#FF7F0E',
                        'CIBERSORT_Abs'='#2CA02C',
                        'ConsensusTME'='#D62728',
                        'DCQ'='#9467BD',
                        'EPIC'='#8C564B',
                        'MCPcounter'='#E377C2',
                        'mMCPcounter'='#7F7F7F',
                        'quanTIseq'='#BCBD22',
                        'seqImmuCC_SVR'='#17BECF',
                        'TIMER'='#AEC7E8',
                        'xCell'='#FFBB78'
                        ))  
)


mat_scaled <- t(scale(t(plotdf_ht)))

ht_final=Heatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  name = "Z-score", 
  top_annotation = col_annotation, 
  left_annotation = row_annotation, 
  row_names_side = "left", 
  column_names_side = "bottom",  
  show_column_names=F,
  #border = 'white',
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width = width * 0.9, height = height * 0.9, 
              gp = gpar(fill = fill, col = NA))
    grid.rect(x, y, width = width * 0.97, height = height * 0.97, 
              gp = gpar(col = "white", fill = NA, lwd = 1))
  },
  column_split=Group_vector

)

# export PDF file
#Fig 1e#####
pdf(".../immune_decov_Bcells_new.pdf", width = 13, height = 5)  
draw(ht_final)
dev.off()  



#Cohen's d value#####

#Cohen's d is better to describe the effect size
library(effsize)

res_MCPcounter %>% rbind(res_EPIC) %>% rbind(res_quanTIseq) %>% 
  rbind(xCell_res) %>% rbind(res_cibersort) %>% rbind(res_cibersort_abs) %>% 
  rbind(res_TIMER) %>% rbind(res_ConsensusTME) %>% rbind(res_ABIS) %>%
  rbind(res_mMCPcounter) %>% rbind(res_seqImmuCC_SVR) %>%
  rbind(res_DCQ)

colnames(res_MCPcounter)
control=c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60)
treated=c(4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)

#MCPcounter
res_MCPcounter_FC=as.data.frame(res_MCPcounter)
cohen_d_values <- apply(res_MCPcounter_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_MCPcounter_FC$Cohen_d <- cohen_d_values

#EPIC
res_EPIC_FC=as.data.frame(res_EPIC)
cohen_d_values <- apply(res_EPIC_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_EPIC_FC$Cohen_d <- cohen_d_values

#quanTIseq
res_quanTIseq_FC=as.data.frame(res_quanTIseq)
cohen_d_values <- apply(res_quanTIseq_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_quanTIseq_FC$Cohen_d <- cohen_d_values

#xCell
res_xCell_FC=as.data.frame(xCell_res)
cohen_d_values <- apply(res_xCell_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_xCell_FC$Cohen_d <- cohen_d_values



#cibersort
res_cibersort_FC=as.data.frame(res_cibersort)
cohen_d_values <- apply(res_cibersort_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_cibersort_FC$Cohen_d <- cohen_d_values

#cibersort_abs
res_cibersort_abs_FC=as.data.frame(res_cibersort_abs)
cohen_d_values <- apply(res_cibersort_abs_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_cibersort_abs_FC$Cohen_d <- cohen_d_values

#TIMER 
res_TIMER_FC=as.data.frame(res_TIMER)
cohen_d_values <- apply(res_TIMER_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_TIMER_FC$Cohen_d <- cohen_d_values

#ConsensusTME 
res_ConsensusTME_FC=as.data.frame(res_ConsensusTME)
cohen_d_values <- apply(res_ConsensusTME_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_ConsensusTME_FC$Cohen_d <- cohen_d_values


#ABIS
res_ABIS_FC=as.data.frame(res_ABIS)
cohen_d_values <- apply(res_ABIS_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_ABIS_FC$Cohen_d <- cohen_d_values


#mMCPcounter
res_mMCPcounter_FC=as.data.frame(res_mMCPcounter)
cohen_d_values <- apply(res_mMCPcounter_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_mMCPcounter_FC$Cohen_d <- cohen_d_values


#seqImmuCC_SVR  
res_seqImmuCC_SVR_FC=as.data.frame(res_seqImmuCC_SVR)
cohen_d_values <- apply(res_seqImmuCC_SVR_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_seqImmuCC_SVR_FC$Cohen_d <- cohen_d_values


#DCQ  
res_DCQ_FC=as.data.frame(res_DCQ)
cohen_d_values <- apply(res_DCQ_FC, 1, function(row) {
  group1 <- row[treated]
  group2 <- row[control]
  cohen.d(group1, group2)$estimate
})
res_DCQ_FC$Cohen_d <- cohen_d_values



#merge
empty_df <- as.data.frame(matrix(ncol = 12, nrow = 10))
colnames(empty_df)=c('MCPcounter','EPIC','quanTIseq','xCell',
                     'CIBERSORT','CIBERSORT_Abs','TIMER','ConsensusTME',
                     'ABIS','mMCPcounter','seqImmuCC_SVR','DCQ')
rownames(empty_df)=c('CD4_T','CD8_T','B_cell','NK','Macrophage',
                     'Monocyte','Neutrophil','DC','Fibroblast','Endothelial')


empty_df$ABIS=c(res_ABIS_FC['T CD4 Memory',"Cohen_d"],
                res_ABIS_FC['T CD8 Naive',"Cohen_d"],
                res_ABIS_FC['B Naive','Cohen_d'],
                res_ABIS_FC['NK','Cohen_d'],
                NA,
                res_ABIS_FC['Monocytes C','Cohen_d'],
                res_ABIS_FC['Neutrophils LD','Cohen_d'],
                res_ABIS_FC['mDCs','Cohen_d'],
                NA,
                NA)



empty_df$CIBERSORT=c(res_cibersort_FC['T cells CD4 memory activated',"Cohen_d"],
                     res_cibersort_FC['T cells CD8',"Cohen_d"],
                     res_cibersort_FC['B cells naive','Cohen_d'],
                     res_cibersort_FC['NK cells activated','Cohen_d'],
                     res_cibersort_FC['Macrophages M2','Cohen_d'],
                     res_cibersort_FC['Monocytes','Cohen_d'],
                     res_cibersort_FC['Neutrophils','Cohen_d'],
                     res_cibersort_FC['Dendritic cells activated','Cohen_d'],
                     NA,
                     NA)


empty_df$CIBERSORT_Abs=c(res_cibersort_abs_FC['T cells CD4 memory activated',"Cohen_d"],
                         res_cibersort_abs_FC['T cells CD8',"Cohen_d"],
                         res_cibersort_abs_FC['B cells naive','Cohen_d'],
                         res_cibersort_abs_FC['NK cells activated','Cohen_d'],
                         res_cibersort_abs_FC['Macrophages M2','Cohen_d'],
                         res_cibersort_abs_FC['Monocytes','Cohen_d'],
                         res_cibersort_abs_FC['Neutrophils','Cohen_d'],
                         res_cibersort_abs_FC['Dendritic cells activated','Cohen_d'],
                         NA,
                         NA)


empty_df$ConsensusTME=c(res_ConsensusTME_FC['T_cells_CD4',"Cohen_d"],
                        res_ConsensusTME_FC['T_cells_CD8',"Cohen_d"],
                        res_ConsensusTME_FC['B_cells','Cohen_d'],
                        res_ConsensusTME_FC['NK_cells','Cohen_d'],
                        res_ConsensusTME_FC['Macrophages_M2','Cohen_d'],
                        res_ConsensusTME_FC['Monocytes','Cohen_d'],
                        res_ConsensusTME_FC['Neutrophils','Cohen_d'],
                        res_ConsensusTME_FC['Dendritic_cells','Cohen_d'],
                        res_ConsensusTME_FC['Fibroblasts','Cohen_d'],
                        res_ConsensusTME_FC['Endothelial','Cohen_d'])

empty_df$DCQ=c(res_DCQ_FC['CD4_Naive',"Cohen_d"],
               res_DCQ_FC['CD8_Naive',"Cohen_d"],
               res_DCQ_FC['B_cells','Cohen_d'],
               res_DCQ_FC['NK_cells','Cohen_d'],
               res_DCQ_FC['Macrophages','Cohen_d'],
               res_DCQ_FC['Monocytes','Cohen_d'],
               res_DCQ_FC['Granulocytes','Cohen_d'],
               res_DCQ_FC['Dendritic_cells','Cohen_d'],
               NA,
               NA)


empty_df$EPIC=c(res_EPIC_FC['CD4_Tcells',"Cohen_d"],
                res_EPIC_FC['CD8_Tcells',"Cohen_d"],
                res_EPIC_FC['Bcells','Cohen_d'],
                res_EPIC_FC['NKcells','Cohen_d'],
                res_EPIC_FC['Macrophages','Cohen_d'],
                NA,
                NA,
                NA,
                res_EPIC_FC['CAFs','Cohen_d'],
                res_EPIC_FC['Endothelial','Cohen_d'])

empty_df$MCPcounter=c(NA,
                      res_MCPcounter_FC['CD8 T cells',"Cohen_d"],
                      res_MCPcounter_FC['B lineage','Cohen_d'],
                      res_MCPcounter_FC['NK cells','Cohen_d'],
                      NA,
                      res_MCPcounter_FC['Monocytic lineage','Cohen_d'],
                      res_MCPcounter_FC['Neutrophils','Cohen_d'],
                      res_MCPcounter_FC['Myeloid dendritic cells','Cohen_d'],
                      res_MCPcounter_FC['Fibroblasts','Cohen_d'],
                      res_MCPcounter_FC['Endothelial cells','Cohen_d'])

empty_df$mMCPcounter=c(NA,
                       res_mMCPcounter_FC['CD8 T cells',"Cohen_d"],
                       res_mMCPcounter_FC['Memory B cells','Cohen_d'],
                       res_mMCPcounter_FC['NK cells','Cohen_d'],
                       res_mMCPcounter_FC['Monocytes / macrophages','Cohen_d'],
                       res_mMCPcounter_FC['Monocytes','Cohen_d'],
                       res_mMCPcounter_FC['Neutrophils','Cohen_d'],
                       NA,
                       res_mMCPcounter_FC['Fibroblasts','Cohen_d'],
                       res_mMCPcounter_FC['Endothelial cells','Cohen_d'])

empty_df$quanTIseq=c(0,
                     res_quanTIseq_FC['T.cells.CD8',"Cohen_d"],
                     res_quanTIseq_FC['B.cells','Cohen_d'],
                     res_quanTIseq_FC['NK.cells','Cohen_d'],
                     res_quanTIseq_FC['Macrophages.M2','Cohen_d'],
                     res_quanTIseq_FC['Monocytes','Cohen_d'],
                     res_quanTIseq_FC['Neutrophils','Cohen_d'],
                     res_quanTIseq_FC['Dendritic.cells','Cohen_d'],
                     NA,
                     NA)

empty_df$seqImmuCC_SVR=c(res_seqImmuCC_SVR_FC['CD4 T Cells',"Cohen_d"],
                         res_seqImmuCC_SVR_FC['CD8 T Cells',"Cohen_d"],
                         res_seqImmuCC_SVR_FC['B Cells','Cohen_d'],
                         res_seqImmuCC_SVR_FC['NK Cells','Cohen_d'],
                         res_seqImmuCC_SVR_FC['Macrophages','Cohen_d'],
                         res_seqImmuCC_SVR_FC['Monocytes','Cohen_d'],
                         res_seqImmuCC_SVR_FC['Neutrophils','Cohen_d'],
                         res_seqImmuCC_SVR_FC['Dendritic Cells','Cohen_d'],
                         NA,
                         NA)

empty_df$TIMER=c(res_TIMER_FC['T_cell.CD4',"Cohen_d"],
                 res_TIMER_FC['T_cell.CD8',"Cohen_d"],
                 res_TIMER_FC['B_cell','Cohen_d'],
                 NA,
                 res_TIMER_FC['Macrophage','Cohen_d'],
                 NA,
                 res_TIMER_FC['Neutrophil','Cohen_d'],
                 res_TIMER_FC['DC','Cohen_d'],
                 NA,
                 NA)


empty_df$xCell=c(res_xCell_FC['CD4+ T-cells',"Cohen_d"],
                 res_xCell_FC['CD8+ T-cells',"Cohen_d"],
                 res_xCell_FC['B-cells','Cohen_d'],
                 res_xCell_FC['NK cells','Cohen_d'],
                 res_xCell_FC['Macrophages','Cohen_d'],
                 res_xCell_FC['Monocytes','Cohen_d'],
                 res_xCell_FC['Neutrophils','Cohen_d'],
                 res_xCell_FC['DC','Cohen_d'],
                 res_xCell_FC['Fibroblasts','Cohen_d'],
                 res_xCell_FC['Endothelial cells','Cohen_d'])





#Fig 1c####

library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)



df_long <- empty_df %>%
  rownames_to_column(var = "Cell_Type") %>%
  pivot_longer(cols = -Cell_Type, names_to = "Method", values_to = "Value") 

df_median <- df_long %>%
  group_by(Cell_Type) %>%
  summarize(median_Value = median(Value, na.rm = TRUE)) %>%
  arrange(median_Value)

df_long$Cell_Type <- factor(df_long$Cell_Type, levels = df_median$Cell_Type)

first_plot=ggplot(df_long, aes(x = Cell_Type, y = Value)) +
  #geom_boxplot(outlier.fill = "white",outlier.color = "white") +
  geom_boxplot(outlier.shape = NA) +  
  #geom_violin(outlier.fill = "white",outlier.color = "white") +
  geom_jitter(aes(color = Method), width = 0.1, size = 2) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(colour = 'black'),
        plot.title = element_text(hjust = 0.5, size = 12)) +
  labs(title = "",
       x = "Cell Type",
       y = "Cohen's d (immunoactivator vs control)") +
  scale_color_d3("category20") +
  coord_flip()




library(FSA)
kruskal_test <- kruskal.test(Value ~ Cell_Type, data = df_long)

dunn_results <- dunnTest(Value ~ Cell_Type, data = df_long, method = "none")

print(dunn_results)

stat_res=dunn_results$res %>% filter(P.adj<0.05)

dunn_p_values <- dunn_results$res$P.adj
dunn_comparisons <- dunn_results$res$Comparison

significant_comparisons <- list()
for (i in seq_along(dunn_p_values)) {
  if (dunn_p_values[i] < 0.05) {
    group_pair <- strsplit(dunn_comparisons[i], " - ")[[1]]
    significant_comparisons[[length(significant_comparisons) + 1]] <- group_pair
  }
}

library(ggsignif)

significant_p_values <- stat_res$P.adj

significant_comparisons=significant_comparisons[c(12,3,10,1,8,7,5,4,11,2,9,6)]
significant_p_values=significant_p_values[c(12,3,10,1,8,7,5,4,11,2,9,6)]

main_plot=first_plot + geom_signif(
  comparisons = significant_comparisons,
  annotations = ifelse(significant_p_values < 0.001, "***", 
                       ifelse(significant_p_values < 0.01, "**", 
                              ifelse(significant_p_values < 0.05, "*", "ns"))),
  step_increase = 0.03,
  tip_length = 0.01,
  vjust = 0.7,
  textsize = 5
) + geom_hline(yintercept=c(0.5),color = "red", linetype = "dashed", size = 0.5)


library(dplyr)
library(ggplot2)
library(cowplot)  

count_above_0.5 <- df_long %>%
  group_by(Cell_Type) %>%
  summarize(count_above_0.5 = sum(Value > 0.5, na.rm = TRUE)) %>%
  arrange(desc(count_above_0.5))


bar_plot <- ggplot(count_above_0.5, aes(x = count_above_0.5, y = Cell_Type)) +
  geom_bar(stat = "identity", fill = "skyblue",width = 0.5) +
  labs(
    title = "", #Count of Cohen's d > 0.5
    x = "Count (Cohen's d > 0.5)",
    y = NULL
  ) +
  theme_bw() +
  #theme_minimal(base_size = 14) +
  theme(
    #plot.title = element_text(hjust = 0.5, size = 16),
    #axis.title = element_text(size = 14),
    #axis.text = element_text(size = 12),
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank()  
  )

library(cowplot)
library(patchwork)
combined_plot <- plot_grid(main_plot, bar_plot, ncol = 2, rel_widths = c(5, 1))

print(combined_plot)

ggsave(".../Cohens_d_plot (immune cell rank).pdf",width = 10,height = 3.5,units = "in")



#Fig 1d#####
colnames(res_ABIS_FC)
control=c(1:3,8:10,14:16,26:28,32:34,38:41,51:54,58:60)
treated=c(4:7,11:13,17:25,29:31,35:37,42:50,55:57,61:63)

library(ggplot2)
library(reshape2)
library(ggpubr)

#average mean
#CD8_T####
CD8_T_df <- as.data.frame(matrix(ncol = 3, nrow = 12))
colnames(CD8_T_df)=c('Control','Immunoactivator','Method')

CD8_T_df[1,]=list(mean(scale(t(res_ABIS_FC['T CD8 Naive',1:63]))[control,]),
                  mean(scale(t(res_ABIS_FC['T CD8 Naive',1:63]))[treated,]),
                  'ABIS')
CD8_T_df[2,]=list(mean(scale(t(res_cibersort_FC['T cells CD8',1:63]))[control,]),
                  mean(scale(t(res_cibersort_FC['T cells CD8',1:63]))[treated,]),
                  'CIBERSORT')
CD8_T_df[3,]=list(mean(scale(t(res_cibersort_abs_FC['T cells CD8',1:63]))[control,]),
                  mean(scale(t(res_cibersort_abs_FC['T cells CD8',1:63]))[treated,]),
                  'CIBERSORT_Abs')
CD8_T_df[4,]=list(mean(scale(t(res_ConsensusTME_FC['T_cells_CD8',1:63]))[control,]),
                  mean(scale(t(res_ConsensusTME_FC['T_cells_CD8',1:63]))[treated,]),
                  'ConsensusTME')
CD8_T_df[5,]=list(mean(scale(t(res_DCQ_FC['CD8_Naive',1:63]))[control,]),
                  mean(scale(t(res_DCQ_FC['CD8_Naive',1:63]))[treated,]),
                  'DCQ')
CD8_T_df[6,]=list(mean(scale(t(res_EPIC_FC['CD8_Tcells',1:63]))[control,]),
                  mean(scale(t(res_EPIC_FC['CD8_Tcells',1:63]))[treated,]),
                  'EPIC')
CD8_T_df[7,]=list(mean(scale(t(res_MCPcounter_FC['CD8 T cells',1:63]))[control,]),
                  mean(scale(t(res_MCPcounter_FC['CD8 T cells',1:63]))[treated,]),
                  'MCPcounter')
CD8_T_df[8,]=list(mean(scale(t(res_mMCPcounter_FC['CD8 T cells',1:63]))[control,]),
                  mean(scale(t(res_mMCPcounter_FC['CD8 T cells',1:63]))[treated,]),
                  'mMCPcounter')
CD8_T_df[9,]=list(mean(scale(t(res_quanTIseq_FC['T.cells.CD8',1:63]))[control,]),
                  mean(scale(t(res_quanTIseq_FC['T.cells.CD8',1:63]))[treated,]),
                  'quanTIseq')
CD8_T_df[10,]=list(mean(scale(t(res_seqImmuCC_SVR_FC['CD8 T Cells',1:63]))[control,]),
                   mean(scale(t(res_seqImmuCC_SVR_FC['CD8 T Cells',1:63]))[treated,]),
                   'seqImmuCC_SVR')
CD8_T_df[11,]=list(mean(scale(t(res_TIMER_FC['T_cell.CD8',1:63]))[control,]),
                   mean(scale(t(res_TIMER_FC['T_cell.CD8',1:63]))[treated,]),
                   'TIMER')
CD8_T_df[12,]=list(mean(scale(t(res_xCell_FC['CD8+ T-cells',1:63]))[control,]),
                   mean(scale(t(res_xCell_FC['CD8+ T-cells',1:63]))[treated,]),
                   'xCell')


CD8_T_df_long <- melt(CD8_T_df, id.vars = "Method", variable.name = "Group", value.name = "Value")

p1=ggplot(CD8_T_df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.5) +  
  geom_jitter(aes(color = Method), position = position_jitter(width = 0), size = 2) +  
  geom_line(aes(group = Method), color = "gray", linetype = "dashed") +  
  stat_compare_means(comparisons = list(c("Control", "Immunoactivator")),
                     method = "wilcox.test", paired = TRUE,
                     tip.length = 0.02,hjust = 0.5,
                     #label.x = 1.5,
                     label = "p.format",
                     aes(label = paste("p =", ..p.format..))  
  ) + 
  scale_fill_manual(values = c('#46A3FF', '#FF5151')) + 
  scale_color_manual(values = c('#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', 
                                '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
                                '#AEC7E8', '#FFBB78')) +  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, colour = 'black',size=9),
        plot.title = element_text(hjust = 0.5, size = 12, colour = "black")) +
  labs(title = "CD8_T", y = "Mean of Z-score", x = "")
p1
ggsave(".../Paired comparison (CD8_T).pdf",width = 4,height = 5,units = "in")



#B_cell######
B_cell_df <- as.data.frame(matrix(ncol = 3, nrow = 12))
colnames(B_cell_df)=c('Control','Immunoactivator','Method')

B_cell_df[1,]=list(mean(scale(t(res_ABIS_FC['B Naive',1:63]))[control,]),
                   mean(scale(t(res_ABIS_FC['B Naive',1:63]))[treated,]),
                   'ABIS')
B_cell_df[2,]=list(mean(scale(t(res_cibersort_FC['B cells naive',1:63]))[control,]),
                   mean(scale(t(res_cibersort_FC['B cells naive',1:63]))[treated,]),
                  'CIBERSORT')
B_cell_df[3,]=list(mean(scale(t(res_cibersort_abs_FC['B cells naive',1:63]))[control,]),
                   mean(scale(t(res_cibersort_abs_FC['B cells naive',1:63]))[treated,]),
                  'CIBERSORT_Abs')
B_cell_df[4,]=list(mean(scale(t(res_ConsensusTME_FC['B_cells',1:63]))[control,]),
                   mean(scale(t(res_ConsensusTME_FC['B_cells',1:63]))[treated,]),
                  'ConsensusTME')
B_cell_df[5,]=list(mean(scale(t(res_DCQ_FC['B_cells',1:63]))[control,]),
                  mean(scale(t(res_DCQ_FC['B_cells',1:63]))[treated,]),
                  'DCQ')
B_cell_df[6,]=list(mean(scale(t(res_EPIC_FC['Bcells',1:63]))[control,]),
                  mean(scale(t(res_EPIC_FC['Bcells',1:63]))[treated,]),
                  'EPIC')
B_cell_df[7,]=list(mean(scale(t(res_MCPcounter_FC['B lineage',1:63]))[control,]),
                  mean(scale(t(res_MCPcounter_FC['B lineage',1:63]))[treated,]),
                  'MCPcounter')
B_cell_df[8,]=list(mean(scale(t(res_mMCPcounter_FC['Memory B cells',1:63]))[control,]),
                  mean(scale(t(res_mMCPcounter_FC['Memory B cells',1:63]))[treated,]),
                  'mMCPcounter')
B_cell_df[9,]=list(mean(scale(t(res_quanTIseq_FC['B.cells',1:63]))[control,]),
                  mean(scale(t(res_quanTIseq_FC['B.cells',1:63]))[treated,]),
                  'quanTIseq')
B_cell_df[10,]=list(mean(scale(t(res_seqImmuCC_SVR_FC['B Cells',1:63]))[control,]),
                    mean(scale(t(res_seqImmuCC_SVR_FC['B Cells',1:63]))[treated,]),
                   'seqImmuCC_SVR')
B_cell_df[11,]=list(mean(scale(t(res_TIMER_FC['B_cell',1:63]))[control,]),
                    mean(scale(t(res_TIMER_FC['B_cell',1:63]))[treated,]),
                   'TIMER')
B_cell_df[12,]=list(mean(scale(t(res_xCell_FC['B-cells',1:63]))[control,]),
                    mean(scale(t(res_xCell_FC['B-cells',1:63]))[treated,]),
                   'xCell')

B_cell_df_long <- melt(B_cell_df, id.vars = "Method", variable.name = "Group", value.name = "Value")

p2=ggplot(B_cell_df_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.5) + 
  geom_jitter(aes(color = Method), position = position_jitter(width = 0), size = 2) + 
  geom_line(aes(group = Method), color = "gray", linetype = "dashed") +  
  stat_compare_means(comparisons = list(c("Control", "Immunoactivator")),
                     method = "wilcox.test", paired = TRUE,
                     tip.length = 0.02,hjust = 0.5,
                     #label.x = 1.5,
                     label = "p.format",
                     aes(label = paste("p =", ..p.format..))  
  ) +  
  scale_fill_manual(values = c('#46A3FF', '#FF5151')) + 
  scale_color_manual(values = c('#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', 
                                '#8C564B', '#E377C2', '#7F7F7F', '#BCBD22', '#17BECF',
                                '#AEC7E8', '#FFBB78')) +  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, colour = 'black',size=9),
        plot.title = element_text(hjust = 0.5, size = 12, colour = "black")) +
  labs(title = "B_cell", y = "Mean of Z-score", x = "")
p2
ggsave(".../Paired comparison (B_cell).pdf",width = 4,height = 5,units = "in")




