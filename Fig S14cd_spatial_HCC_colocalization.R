library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggsci)
library(ggpubr)

#JOH,2023 PMID: 36708811
#https://data.mendeley.com/datasets/skrx2fz79n/1
#P7 P9 P10 responder


P1T <- readRDS(gzfile(".../STData of HCC/P1T_Spatial.rds.gz"))
P3T <- readRDS(gzfile(".../STData of HCC/P3T_Spatial.rds.gz"))
P5T <- readRDS(gzfile(".../STData of HCC/P5T_Spatial.rds.gz"))
P7T <- readRDS(gzfile(".../STData of HCC/P7T_Spatial.rds.gz"))
P8T <- readRDS(gzfile(".../STData of HCC/P8T_Spatial.rds.gz"))
P9T <- readRDS(gzfile(".../STData of HCC/P9T_Spatial.rds.gz"))
P10T <- readRDS(gzfile(".../STData of HCC/P10T_Spatial.rds.gz"))
P11T <- readRDS(gzfile(".../STData of HCC/P11T_Spatial.rds.gz"))

#P1T
keep_img <- "image_P11_T"
P1T@images <- P1T@images[keep_img]
names(P1T@images) <- "image"
P1T


#Figure S14c#####
SpatialFeaturePlot(P11T, features = c("CD79A",'HAVCR1',"IL10"),alpha = (1)) 
ggsave(".../spatial analysis/Represent_No_Response.pdf",width = 8,height = 4,units = "in")



SpatialFeaturePlot(P7T, features = c("CD79A", 'HAVCR1',"IL10")) 
ggsave(".../spatial analysis/Represent_Response.pdf",width = 8,height = 4,units = "in")



#####colocalization, CD79A/B and HAVCR1#####
coexpr_spots <- WhichCells(
  object = P1T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P1T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),    
  cols.highlight  = c("red","lightgrey"),        
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )



coexpr_spots <- WhichCells(
  object = P3T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P3T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),    
  cols.highlight  = c("lightgrey","red"),         
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )


coexpr_spots <- WhichCells(
  object = P5T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P5T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),    
  cols.highlight  = c("lightgrey","red"),          
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )


coexpr_spots <- WhichCells(
  object = P8T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P8T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),    
  cols.highlight  = c("lightgrey","red"),         
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )


coexpr_spots <- WhichCells(
  object = P11T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P11T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),    
  cols.highlight  = c("lightgrey","red"),          
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(".../spatial analysis/Represent_colocalization_P11T.pdf",
       width = 4,height = 4,units = "in")



coexpr_spots <- WhichCells(
  object = P7T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P7T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),   
  cols.highlight  = c("red","lightgrey"),          
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(".../spatial analysis/Represent_colocalization_P7T.pdf",
       width = 4,height = 4,units = "in")



coexpr_spots <- WhichCells(
  object = P9T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P9T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),     
  cols.highlight  = c("red","lightgrey"),         
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )


coexpr_spots <- WhichCells(
  object = P10T,
  expression = (CD79A > 0 | CD79B > 0) & HAVCR1 > 0
)
SpatialPlot(
  object          = P10T,
  cells.highlight = coexpr_spots,
  cols            = c("lightgrey", "lightgrey"),    
  cols.highlight  = c("red","lightgrey"),          
  pt.size.factor  = 1
) + 
  ggtitle("Colocalization of CD79A/B & HAVCR1 Spots")+
  theme(
    plot.title = element_text(hjust = 0.5)
  )



#Figure S14d#####
#statistical comparison#####


#B cell dots
spots_both <- WhichCells(P1T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots:", n_both, "\n")
392

spots_both <- WhichCells(P3T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
133

spots_both <- WhichCells(P5T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
115

spots_both <- WhichCells(P8T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
1576

spots_both <- WhichCells(P11T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
477




spots_both <- WhichCells(P7T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
845

spots_both <- WhichCells(P9T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
923

spots_both <- WhichCells(P10T, expression = CD79A > 0 | CD79B > 0)
n_both <- length(spots_both)
cat("CD79A | CD79B spots：", n_both, "\n")
406


df2=data.frame(
  patient=c('PT1','PT3','PT5','PT8','PT11','PT7','PT9','PT10'),
  colocalization_ratio=c(1/392*100,0,0,0,0,7/845*100,4/923*100,2/406*100),
  group=c(rep('Non-response',5),rep('Response',3))
)




ggplot(df2, aes(x = group, y = colocalization_ratio, color = group)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "gray", linetype = "dashed") +
  stat_compare_means(method = "t.test", label = "p", size = 4, label.x.npc  = "middle") +
  scale_color_manual(values = c("Non-response" = "#2E9FDF", "Response" = "#FF6B6B")) +
  labs(title = "Colocalization of TIM-1 and B cells",
       x = "Group",
       y = "Colocalization Rate among CD79A/B+ Spots (%)",
       color = "Group") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12, color = 'black'),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )+
  NoLegend()

ggsave(".../spatial analysis/Colocalization_2.pdf",width = 3,height = 4.5,units = "in")





