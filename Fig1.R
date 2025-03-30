library(Seurat)
library(ggplot2)
install.packages('writexl')
library(writexl)
library(dplyr)
install.packages("clustree")
library(clustree)


exhaustedT <- subset(Tcells_humancancer, subset = CD3E > 0 & CD8A > 0 & PDCD1 > 0 & MKI67 == 0 & CD4 == 0, slot = 'counts')
exhaustedT <- FindNeighbors(exhaustedT, dims = 1:30, reduction = 'harmony')
exhaustedT <- FindClusters(exhaustedT, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2))
clustree(exhaustedT, prefix = "SCT_snn_res.")
exhaustedT <- RunUMAP(exhaustedT, dims = 1:30, reduction = 'harmony')

###Fig1D
Idents(exhaustedT) <- exhaustedT$orig.idents
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("PA","Thyroid","Lung","Pancreatic","Ovarian","CRC","MPNST","Melanoma"))
plot <- DimPlot(exhaustedT,
        reduction = "umap",
        combine = TRUE, 
        label = F,
        pt.size = 0.3,
        cols = c("#D53E4F","#B58CC2","#3288BD","#FC8D59","#99D594","#FEE08B","#F27D8B","#75AADB"))
ggsave("Cancers_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


#FigS1A-C
plot <- FeaturePlot(exhaustedT, 'CD3E', cols = c('lightblue','steelblue'), min.cutoff = 1)
plot <- FeaturePlot(exhaustedT, 'CD8A', cols = c('lightblue','steelblue'), min.cutoff = 1)
plot <- FeaturePlot(exhaustedT, 'PDCD1', cols = c('lightblue','steelblue'), min.cutoff = 1)
ggsave("PDCD1.pdf", plot = plot, width = 6, height = 5, dpi = 300)



exhaustedT <- PrepSCTFindMarkers(exhaustedT)
exhaustedT.markers <- FindAllMarkers(exhaustedT, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE)
exhaustedT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write_xlsx(exhaustedT.markers,"exhaustedT_markers.xlsx")


Idents(exhaustedT) <- exhaustedT$SCT_snn_res.0.5

rename<- c("1","10")
exhaustedT@meta.data$humanC<- ifelse(exhaustedT@meta.data$SCT_snn_res.0.5 %in% rename, "C3", NA)
rename<- c("4")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.5 %in% rename, "C2", exhaustedT@meta.data$humanC)
rename<- c("0","2","3","6","7","8","11")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.5 %in% rename, "C1", exhaustedT@meta.data$humanC)
rename<- c("5")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.5 %in% rename, "C4", exhaustedT@meta.data$humanC)
rename<- c("9")
exhaustedT@meta.data$humanC <- ifelse(exhaustedT@meta.data$SCT_snn_res.0.5 %in% rename, "C5", exhaustedT@meta.data$humanC)


Idents(exhaustedT) <- exhaustedT$humanC
C4 <- subset(exhaustedT, idents = "C4")
C4 <- FindClusters(C4, resolution = 0.8)
C4 <- RunUMAP(C4, dims = 1:30, reduction = 'harmony')
DimPlot(C4, label = TRUE)
FeaturePlot(C4, c('SELL','CCR7','IL7R'))
FeaturePlot(C4, c('HSPA6','HSPA1A','BAG3'))
C4.markers <- FindAllMarkers(C4, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE)
exhaustedT$humanC_original <- exhaustedT$humanC
Idents(C4) <- "seurat_clusters"
cells_to_move <- WhichCells(C4, idents = "3") 
exhaustedT$humanC[cells_to_move] <- "C3"
DimPlot(exhaustedT, label = TRUE)

#Fig1C
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("C1","C2","C3","C4","C5"))
plot <- DimPlot(
  exhaustedT, 
  combine = TRUE, 
  label = F,pt.size = 0.3,
  cols = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
)
ggsave("Labeled_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


table(exhaustedT@active.ident, exhaustedT$humanC)

#Fig1E-F
library(ggplot2)
library(RColorBrewer)
data <- readxl::read_excel("Cell proportion.xlsx")
data$Cluster <- factor(data$Cluster,levels = c("C1","C2","C3","C4","C5"))
data$Cancer <- factor(data$Cancer,levels = c("Brain tumor","Thyroid cancer","Lung cancer","Pancreatic cancer","Ovarian cancer","Colorectal cancer","Nerve sheath cancer","Skin cancer"))
p1 <- ggplot(data, aes(fill= Cancer, y=Count, x=Cluster)) + 
  geom_bar(position="fill", stat="identity") + 
  scale_fill_manual(values = c("#D53E4F","#B58CC2","#3288BD","#FC8D59","#99D594","#FEE08B","#F27D8B","#75AADB"))+
  theme(panel.background = element_rect(fill = "white"))

ggsave(filename = "human_Cbar.pdf", plot = p1, device = "pdf", width = 8, height = 6)  

p2 <- ggplot(data, aes(fill = Cluster, y = Count, x = Cancer)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = "cancer_Cbar.pdf", plot = p2, device = "pdf", width = 8, height = 6)  


library(ggplot2)
library(ggthem)
data <- readxl::read_excel("Exhausted T cell percenrage.xlsx")
data$Cancer <- factor(data$Cancer,levels = c("PA","Thyroid cancer","Lung cancer","Pancreatic cancer","Ovarian cancer","Colorectal cancer","MPNST","Melanoma"))


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


