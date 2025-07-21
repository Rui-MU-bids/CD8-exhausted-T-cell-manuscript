library(Seurat)
library(ggplot2)
library(writexl)
library(dplyr)
library(clustree)

exhaustedT <- subset(Tcells, subset = Cd3e > 0 & Cd8a > 0 & Pdcd1 > 0 & Mki67 == 0 & Cd4 == 0, slot = 'counts')
exhaustedT <- FindNeighbors(exhaustedT, dims = 1:30, reduction = 'harmony')
exhaustedT <- FindClusters(exhaustedT, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
clustree(exhaustedT, prefix = "SCT_snn_res.")
exhaustedT <- FindClusters(exhaustedT, resolution = 0.8)

Idents(exhaustedT) <- exhaustedT$SCT_snn_res.0.8
exhaustedT <- RunUMAP(exhaustedT, dims = 1:30, reduction = 'harmony')
DimPlot(exhaustedT, label = TRUE)
DimPlot(exhaustedT, group.by = 'orig.idents')


exhaustedT <- subset(exhaustedT, idents = c("0","1","2","3","4","5","6","7","8","10")) #remove macrophage doublets cluster

exhaustedT <- FindNeighbors(exhaustedT, dims = 1:30, reduction = 'harmony')
exhaustedT <- FindClusters(exhaustedT, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
clustree(exhaustedT, prefix = "SCT_snn_res.")

Idents(exhaustedT) <- exhaustedT$SCT_snn_res.0.8
exhaustedT <- RunUMAP(exhaustedT, dims = 1:30, reduction = 'harmony')
DimPlot(exhaustedT, label = TRUE)

Idents(exhaustedT) <- exhaustedT$SCT_snn_res.0.8
exhaustedT <- subset(exhaustedT, idents = c("0","1","2","3","4","5","6","7","8")) #NK doublets cluster removed

mice_inte_exhaustedT <- FindNeighbors(mice_inte_exhaustedT, dims = 1:30, reduction = 'harmony')
mice_inte_exhaustedT <- FindClusters(mice_inte_exhaustedT, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) #FigS4D
clustree(mice_inte_exhaustedT, prefix = "SCT_snn_res.")

Idents(exhaustedT) <- exhaustedT$SCT_snn_res.0.8

rename<- c("1","4")
exhaustedT@meta.data$humanC<- ifelse(exhaustedT@meta.data$SCT_snn_res.0.8 %in% rename, "C3", NA)
rename<- c("2")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.8 %in% rename, "C2", exhaustedT@meta.data$humanC)
rename<- c("5","6","0")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.8 %in% rename, "C1", exhaustedT@meta.data$humanC)
rename<- c("8","3")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.8 %in% rename, "C6", exhaustedT@meta.data$humanC)
rename<- c("7")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.8 %in% rename, "C5", exhaustedT@meta.data$humanC)

Idents(exhaustedT) <- exhaustedT$humanC

C3 <- subset(exhaustedT, idents = "C3")
C3 <- FindClusters(C3, resolution = 1.2)
C3 <- RunUMAP(C3, dims = 1:30, reduction = 'harmony')
DimPlot(C3, label = TRUE)
VlnPlot(C3,c("Ccl3","Ccl4","Ccl5", "Ccr2"))
VlnPlot(C3, c("Tcf7","Sell"))
Idents(C3) <- "seurat_clusters"
cells_to_move <- WhichCells(C3, idents = c("5"))
exhaustedT$humanC <- exhaustedT$humanC
exhaustedT$humanC[cells_to_move] <- "C1"
DimPlot(exhaustedT, label = TRUE)

Idents(exhaustedT) <- exhaustedT$humanC
C5 <- subset(exhaustedT, idents = "C5")
C5 <- FindClusters(C5, resolution = 1.2)
C5 <- RunUMAP(C5, dims = 1:30, reduction = 'harmony')
DimPlot(C5, label = TRUE)
VlnPlot(C5,c("Tcf7","Sell"),)
FeaturePlot(C5, c("Tcf7","Sell")) 


exhaustedT <- PrepSCTFindMarkers(exhaustedT)
exhaustedT.markers <- FindAllMarkers(exhaustedT, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE)
exhaustedT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write_xlsx(exhaustedT.markers,"mice_exhaustedT_markers.xlsx")



#Fig4B
library(ggplot2)
library(ggthem)
data <- readxl::read_excel("Exhausted T cell percentage.xlsx")
data$Cancer <- factor(data$Cancer,levels = c("Brain tumor","Skin cancer","Pancreatic cancer","Colorectal cancer"))

line_size <- 1
base_size <- 12
axis_text_rel_size <- -1
title_text_rel_size <- 2

plot <- ggplot(data, aes(x = Cancer, y = `CD8 exhausted T cell percentage`)) +
  stat_summary(fun = mean, geom = "bar", fill = "white",color = "black", width = 0.7) +
  geom_jitter(width = 0.1) +
  scale_y_continuous(limits = c(-0.00001, 0.5), expand = expansion(mult = c(0, 0))) + 
  labs(title = "CD8 exhausted T cell percentage", x = "", y = "Percentage") + 
  theme_foundation(base_size = base_size, base_family = "sans") + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    text = element_text(colour = "black"),
    plot.title = element_text(face = "bold", 
                              size = rel((title_text_rel_size + base_size) / base_size), hjust = 0.5),
    axis.line = element_line(colour="black", size = line_size),
    axis.ticks = element_line(colour="black", size = line_size),
    axis.title = element_text(face = "bold", size = rel(1)),
    axis.title.y = element_text(angle = 90, vjust = 2),
    axis.title.x = element_text(vjust = -0.2),
    axis.text = element_text(face = "bold", size = rel((axis_text_rel_size + base_size) / base_size)),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.background = element_blank()
  )
ggsave(filename = "percentage_mice.pdf", plot = plot, device = "pdf", width = 8, height = 6)  



#Fig4C
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("C1","C2","C3","C5","C6"))
plot <- DimPlot(
  exhaustedT,
  combine = TRUE, 
  label = F,
  pt.size = 0.3,
  cols = c("#E69F00", "#56B4E9", "#009E73", "#0072B2","#C26841")
)

ggsave("miceLabeled_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


#Fig4D
Idents(exhaustedT) <- exhaustedT$cancertype
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("FMC","Skin","Pancreatic","Colorectal"))
plot <- DimPlot(exhaustedT,
                reduction = "umap",
                combine = TRUE, 
                label = F,
                pt.size = 0.3,
                cols = c("#D53E4F","#75AADB","#FC8D59","#FEE08B"))
ggsave("miceCancers_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


#Fig4E-F
data <- readxl::read_excel("Cell proportion.xlsx")
data$Cluster <- factor(data$Cluster,levels = c("C1","C2","C3","C5","C6"))
data$Cancer <- factor(data$Cancer,levels = c("Brain tumor","Skin cancer","Pancreatic cancer","Colorectal cancer"))
p1 <- ggplot(data, aes(fill= Cancer, y=Count, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("#D53E4F","#75AADB","#FC8D59","#FEE08B"))+
  theme(panel.background = element_rect(fill = "white"))
ggsave(filename = "Cbar.pdf", plot = p1, device = "pdf", width = 8, height = 6) 

p2 <- ggplot(data, aes(fill = Cluster, y = Count, x = Cancer)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#0072B2","#C26841")) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = "cancer_Cbar.pdf", plot = p2, device = "pdf", width = 8, height = 6)  


#Fig4G
library(ggplot2)
library(dplyr)
library(forcats)

dat <- read_xlsx("pathwayresult.xlsx")
dat <- as.data.frame(dat)
dat$Term <- factor(dat$Term)

custom_palette <- c("#08306B", "#2171B5", "#6BAED6", "#C6DBEF")

p1 <- ggplot(data = dat, aes(x = `Fold Enrichment`, y = fct_reorder(Term, `Fold Enrichment`), fill = FDR)) +
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Fold Enrichment", y = NULL) +
  theme_bw() +
  scale_fill_gradientn(colors = custom_palette) +
  facet_grid(rows = vars(Cluster), scales = "free_y", switch = "y") +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.placement = "outside"
  )

print(p1)
ggsave(filename = "pathway.pdf", plot = p1, device = "pdf", width = 10, height = 18)




#FigS4A-C
p1 <- FeaturePlot(exhaustedT, 'Cd3e', slot = 'data', cols = c("lavender", "#bd7fef"), min.cutoff = 1)
p2 <- FeaturePlot(exhaustedT, 'Cd8a', slot = 'data', cols = c("lavender", "#bd7fef"), min.cutoff = 1)
p3 <- FeaturePlot(exhaustedT, 'Pdcd1', slot = 'data',cols = c("lavender", "#bd7fef"), min.cutoff = 1)
ggsave(filename = "Pdcd1.pdf", plot = p3, device = "pdf", width = 8, height = 6)


#FigS5A
Exhaustion <- c("Klf4","Elf4","Prdm1","Id2","Jun","Tbx21","Tcf7","Ccr5","Cst7","Nkg7","Cd28","Cd27","Gzmk","Vsir","Tnfrsf18",'Tnfrsf9','Tnfrsf4','Cd226','Eomes','Itgae','Alcam','Btla','Cd244','Entpd1','Havcr2','Ctla4','Lag3','Pdcd1','Tox','Tigit')
Idents(exhaustedT) <- exhaustedT$humanC
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("C1","C2","C3","C5","C6"))
plot <- DotPlot(exhaustedT, features = Exhaustion, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()

#FigS5C
subtypes <- c('Ccl3',"Ccl4","Ccl5","Klrd1","Ifng","Cd200r4",
              "S100a4","S100a6","S100a13","Lgals1","Lgals3",
              "Il7r","Tcf7","Ccr7","Ccr8","Sell","Cd83","Il2ra","Slamf6",
              "Ifit1","Ifit2","Ifit3","Ifi208","Ifit3b","Oas3","Oasl2","Isg15","Stat2","Rsad2",
              "Cdca8","Birc5","Cdc45","Mcm3","Cdk1","Tacc3","Gins2","Cks1b","Rad51")
DefaultAssay(exhaustedT) <- 'SCT'
avg_exp <- AverageExpression(exhaustedT, assays = "SCT", features = subtypes, group.by = 'humanC', return.seurat = T)
my_colors <- rev(brewer.pal(n = 7, name = "RdBu"))
DoHeatmap(
  object   = avg_exp,   
  features = subtypes,      
  group.by = "humanC",        
  draw.lines = FALSE           
) +
  
  scale_fill_gradientn(
    colors = my_colors,
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  )

#FigS5D-E
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("C1","C2","C3","C5","C6"))
effect <- c("Runx3","Klrg1","Fas","Spn","Cd44","Osr2","Cxcr1")
others <- c("Cxcr5","Slamf6","Prdm1","Cd69","Cd101")
plot <- DotPlot(exhaustedT, features = effect, dot.scale = 10, scale = F, cols = c("lightgrey", "#bd7fef"),scale.min = 0, scale.max = 100)+coord_flip()+scale_color_gradientn(colors = c("lightgrey", "#bd7fef"), limits = c(0, 5))
plot <- DotPlot(exhaustedT, features = others, dot.scale = 10, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave(filename = "mice_eff.pdf", plot = plot, device = "pdf", width = 7, height = 5) 