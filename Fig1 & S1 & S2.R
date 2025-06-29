library(Seurat)
library(ggplot2)
library(writexl)
library(dplyr)
library(clustree)


#Fig1B
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
ggsave(filename = "percentage_human.pdf", plot = plot, device = "pdf", width = 8, height = 6)  

exhaustedT <- subset(Tcells_humancancer, subset = CD3E > 0 & CD8A > 0 & PDCD1 > 0 & MKI67 == 0 & CD4 == 0, slot = 'counts')
exhaustedT <- FindNeighbors(exhaustedT, dims = 1:30, reduction = 'harmony')
exhaustedT <- FindClusters(exhaustedT, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2))
clustree(exhaustedT, prefix = "SCT_snn_res.")  #FigS1D
exhaustedT <- RunUMAP(exhaustedT, dims = 1:30, reduction = 'harmony')



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


#Fig1D
Idents(exhaustedT) <- exhaustedT$orig.idents
Idents(exhaustedT) <- factor(Idents(exhaustedT),levels = c("PA","Thyroid","Lung","Pancreatic","Ovarian","CRC","MPNST","Melanoma"))
plot <- DimPlot(exhaustedT,
        reduction = "umap",
        combine = TRUE, 
        label = F,
        pt.size = 0.3,
        cols = c("#D53E4F","#B58CC2","#3288BD","#FC8D59","#99D594","#FEE08B","#F27D8B","#75AADB"))
ggsave("Cancers_exh.pdf", plot = plot, width = 6, height = 5, dpi = 300)


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



exhaustedT <- PrepSCTFindMarkers(exhaustedT)
exhaustedT.markers <- FindAllMarkers(exhaustedT, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi = FALSE)
exhaustedT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write_xlsx(exhaustedT.markers,"exhaustedT_markers.xlsx")


#FigS1A-C
plot <- FeaturePlot(exhaustedT, 'CD3E', cols = c("lavender", "#bd7fef"), min.cutoff = 1)
plot <- FeaturePlot(exhaustedT, 'CD8A', cols = c("lavender", "#bd7fef"), min.cutoff = 1)
plot <- FeaturePlot(exhaustedT, 'PDCD1', cols = c("lavender", "#bd7fef"), min.cutoff = 1)
ggsave("PDCD1.pdf", plot = plot, width = 6, height = 5, dpi = 300)

#FigS2A-B
eff <- c("RUNX3","KLRG1","FAS","SPN","CD44","OSR2","CXCR1")
plot <- DotPlot(exhaustedT, features = eff, dot.scale = 10, scale = F, cols = c("lightgrey", "#bd7fef"),scale.min = 0, scale.max = 100)+coord_flip()+scale_color_gradientn(colors = c("lightgrey", "#bd7fef"), limits = c(0, 5))
ggsave("eff.pdf", plot = plot, width = 6, height = 5, dpi = 300)

Development <- c("CXCR5","SLAMF6","PRDM1",'CD69','CD101')
plot <- DotPlot(exhaustedT, features = Development, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("develop.pdf", plot = plot, width = 6, height = 5, dpi = 300)

#FigS2C

CD8T_human<- subset(Tcells_humancancer, subset = CD3E > 0 & CD8A > 0 & CD4 == 0, slot = 'counts')

humanC_values <- exhaustedT@meta.data$humanC
subset_cells <- colnames(exhaustedT)
CD8T_human$subtype[subset_cells] <- humanC_values
Mki67 <- subset(CD8T_human, subset = Mki67 > 0)
subset_cells <- colnames(Mki67)
CD8T_mice$subtype[subset_cells] <- 'Mki67'
CD8T_mice$subtype[is.na(CD8T_mice$subtype)] <- 'others'

table(CD8T_human$subtype)
table(exhaustedT$humanC) #same

plot = DimPlot(CD8T_human,
        reduction = "umap",
        combine = TRUE, 
        label = F,
        pt.size = 0.3,
        cols = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2","hotpink","thistle"))
ggsave("Human_CD8_umap.pdf", plot = plot, width = 6, height = 5, dpi = 300)


#FigS2D
Exhaustion <- c("KLF4","ELF4","PRDM1","ID2","JUN","TBX21","TCF7","CCR5","CST7","NKG7","CD28","CD27","GZMK","VSIR","TNFRSF18",'TNFRSF9','TNFRSF4','CD226','EOMES','ITGAE','ALCAM','BTLA','CD244','ENTPD1','HAVCR2','CTLA4','LAG3','TOX','TIGIT','PDCD1')
plot <- DotPlot(CD8T_human, features = Exhaustion, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("more_exh_markers_human.pdf", plot = plot, width = 8, height = 10, dpi = 300)

#FigS2E
subtypes <- c('CCL3',"CCL4","CXCL13","KLRD1","KIR2DL4","EGR2","IFNG","IGKV1-5","GNLY","LAT2",
              "COL1A1","COL6A1","CCN1","CCN2","CXCL14","CAV1","FN1","MMP9","CTSK","PCOLCE",
              "IL7R","TCF7","CCR7","SELL","S1PR1","CD28",
              "HSPA1A","HSPA1B","HSPA6","HSPA8","HSPH1","HSPD1","DNAJA1","DNAJB1","SERPINH1","BAG3",
              "IFIT1","IFIT2","IFIT3","IFI27","OAS1","OAS3","MX1","MX2","HERC5","ISG15")
DefaultAssay(CD8T_human) <- 'SCT'
avg_exp <- AverageExpression(CD8T_human, assays = "SCT", features = subtypes, group.by = 'subtype', return.seurat = T)
DoHeatmap(avg_exp, features = subtypes,group.by = 'subtype', draw.lines = F)

my_colors <- rev(brewer.pal(n = 7, name = "RdBu"))
DoHeatmap(
  object   = avg_exp,   
  features = subtypes,      
  group.by = "subtype",        
  draw.lines = FALSE           
) +
  
  scale_fill_gradientn(
    colors = my_colors,
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  )

#FigS2F
CD8_labels <- as.character(Idents(CD8T_human))
CD8_labels[CD8_labels %in% c("C1", "C2", "C3", "C4", "C5")] <- "exh"
CD8T_human$CD8_labels <- CD8_labels
Idents(CD8T_human) <- CD8T_human$CD8_labels
Idents(CD8T_human) <- factor(Idents(CD8T_human),levels = c("exh","MKI67",'others'))

CD8T_human <- AddModuleScore(
  object = CD8T_human,
  features = Exhaustion,
  name = 'exhausted_score',
  seed = 123
)

vln <- VlnPlot(CD8T_human, 'exhausted_score1', pt.size = 0) + coord_cartesian(ylim = c(-0.8, 1))
vln + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "red")


exhaustion_data <- data.frame(
  score = CD8T_human$exhausted_score_human1,
  group = CD8T_human$CD8_labels
)

exhaustion_data$group <- factor(exhaustion_data$group, levels = c("exh", "MKI67", "others"))

groups <- levels(exhaustion_data$group)
pairwise_results <- data.frame()
# Test normality for each group
groups <- levels(exhaustion_data$group)
normality_results <- data.frame()

for(group in groups) {
  group_data <- exhaustion_data$score[exhaustion_data$group == group]
  n <- length(group_data)
  
  # Choose appropriate normality test based on sample size
  if(n >= 3 & n <= 5000) {
    test_result <- shapiro.test(group_data)
    test_name <- "Shapiro-Wilk"
  } else if(n > 5000) {
    test_result <- ks.test(group_data, "pnorm", mean(group_data), sd(group_data))
    test_name <- "Kolmogorov-Smirnov"
  } else {
    test_result <- list(p.value = NA)
    test_name <- "Sample too small"
  }
  
  normality_results <- rbind(normality_results, data.frame(
    group = group,
    n = n,
    test = test_name,
    p_value = test_result$p.value,
    normal = ifelse(is.na(test_result$p.value), "Unknown", 
                    ifelse(test_result$p.value > 0.05, "Yes", "No"))
  ))
}

print(normality_results)


wilcox_results <- data.frame()

comparisons <- list(
  c("exh", "MKI67"),
  c("exh", "others"), 
  c("MKI67", "others")
)

for(comp in comparisons) {
  group1_data <- exhaustion_data$score[exhaustion_data$group == comp[1]]
  group2_data <- exhaustion_data$score[exhaustion_data$group == comp[2]]
  
  wilcox_result <- wilcox.test(group1_data, group2_data)
  
  wilcox_results <- rbind(wilcox_results, data.frame(
    comparison = paste(comp[1], "vs", comp[2]),
    median_group1 = median(group1_data),
    median_group2 = median(group2_data),
    wilcox_p = wilcox_result$p.value
  ))
}

wilcox_results$wilcox_p_bonf <- p.adjust(wilcox_results$wilcox_p, method = "bonferroni")
print(wilcox_results)


