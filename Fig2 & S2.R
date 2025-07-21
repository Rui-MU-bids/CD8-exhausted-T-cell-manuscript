library(Seurat)
library(ggplot2)
library(writexl)
library(dplyr)
library(clustree)
library(RColorBrewer)

#FigS2D
exhaustedT <- subset(Tcells_humancancer, subset = CD3E > 0 & CD8A > 0 & PDCD1 > 0 & MKI67 == 0 & CD4 == 0, slot = 'counts')
exhaustedT <- FindNeighbors(exhaustedT, dims = 1:30, reduction = 'harmony')
exhaustedT <- FindClusters(exhaustedT, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2))
clustree(exhaustedT, prefix = "SCT_snn_res.")
exhaustedT <- RunUMAP(exhaustedT, dims = 1:30, reduction = 'harmony')

###Fig2B
Idents(exhaustedT) <- exhaustedT$orig.idents
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("PA","Thyroid","Lung","Pancreatic","Ovarian","CRC","MPNST","Melanoma"))
plot <- DimPlot(exhaustedT,
        reduction = "umap",
        combine = TRUE, 
        label = F,
        pt.size = 0.3,
        cols = c("#D53E4F","#B58CC2","#3288BD","#FC8D59","#99D594","#FEE08B","#F27D8B","#75AADB"))
ggsave("Cancers_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


#FigS2A-C
plot <- FeaturePlot(human_inte_exhaustedT, 'CD3E', cols = c('lavender','#bd7fef'), min.cutoff = 1)
plot <- FeaturePlot(exhaustedT, 'CD8A', cols = c('lavender','#bd7fef'), min.cutoff = 1)
plot <- FeaturePlot(exhaustedT, 'PDCD1', cols = c('lavender','#bd7fef'), min.cutoff = 1)
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

#Fig2A
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("C1","C2","C3","C4","C5"))
plot <- DimPlot(
  exhaustedT, 
  combine = TRUE, 
  label = F,pt.size = 0.3,
  cols = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
)
ggsave("Labeled_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


table(exhaustedT@active.ident, exhaustedT$humanC)

#Fig2C-D
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






