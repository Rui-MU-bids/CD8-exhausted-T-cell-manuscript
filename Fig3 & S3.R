#Fig3A
subtypes <- c('CCL3',"CCL4","CXCL13","KLRD1","KIR2DL4","EGR2","IFNG","IGKV1-5","GNLY","LAT2",
              "COL1A1","COL6A1","CCN1","CCN2","CXCL14","CAV1","FN1","MMP9","CTSK","PCOLCE",
              "IL7R","TCF7","CCR7","SELL","S1PR1","CD28",
              "HSPA1A","HSPA1B","HSPA6","HSPA8","HSPH1","HSPD1","DNAJA1","DNAJB1","SERPINH1","BAG3",
              "IFIT1","IFIT2","IFIT3","IFI27","OAS1","OAS3","MX1","MX2","HERC5","ISG15")
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

#Fig3B
Exhaustion <- c("KLF4","ELF4","PRDM1","ID2","JUN","TBX21","TCF7","CCR5","CST7","NKG7","CD28","CD27","GZMK","VSIR","TNFRSF18",'TNFRSF9','TNFRSF4','CD226','EOMES','ITGAE','ALCAM','BTLA','CD244','ENTPD1','HAVCR2','CTLA4','LAG3','TOX','TIGIT','PDCD1')
plot <- DotPlot(exhaustedT, features = Exhaustion, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("exh_markers.pdf", plot = plot, width = 5, height = 2.5, dpi = 300)

#Fig3C
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


#FigS3A-B
eff <- c("RUNX3","KLRG1","FAS","SPN","CD44","OSR2","CXCR1")
plot <- DotPlot(exhaustedT, features = eff, dot.scale = 10, scale = F, cols = c("lightgrey", "#bd7fef"),scale.min = 0, scale.max = 100)+coord_flip()+scale_color_gradientn(colors = c("lightgrey", "#bd7fef"), limits = c(0, 5))
ggsave("eff.pdf", plot = plot, width = 6, height = 5, dpi = 300)

Development <- c("CXCR5","SLAMF6","PRDM1",'CD69','CD101')
plot <- DotPlot(exhaustedT, features = Development, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("develop.pdf", plot = plot, width = 6, height = 5, dpi = 300)

#FigS3C
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


#FigS3D
Exhaustion <- c("KLF4","ELF4","PRDM1","ID2","JUN","TBX21","TCF7","CCR5","CST7","NKG7","CD28","CD27","GZMK","VSIR","TNFRSF18",'TNFRSF9','TNFRSF4','CD226','EOMES','ITGAE','ALCAM','BTLA','CD244','ENTPD1','HAVCR2','CTLA4','LAG3','PDCD1','TOX','TIGIT')
plot <- DotPlot(CD8T_human, features = Exhaustion, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("more_exh_markers_human.pdf", plot = plot, width = 8, height = 10, dpi = 300)

#FigS3E
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

#FigS3F
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
