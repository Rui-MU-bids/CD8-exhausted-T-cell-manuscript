library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(dplyr)

#FigS6A
Idents(CD8T_mice) <- CD8T_mice$subtype
Idents(CD8T_mice) <- factor(Idents(CD8T_mice),levels = c("C1","C2","C3","C5","C6","Mki67",'others'))
plot <- DimPlot(
  CD8T_mice,
  combine = TRUE, 
  label = F,
  pt.size = 0.3,
  cols = c("#E69F00", "#56B4E9", "#009E73", "#0072B2","#C26841","hotpink", "thistle")
)
ggsave("Mice_CD8_umap.pdf", plot = plot, width = 6, height = 5, dpi = 300)


#FigS6B
Exhaustion <- c("Klf4","Elf4","Prdm1","Id2","Jun","Tbx21","Tcf7","Ccr5","Cst7","Nkg7","Cd28","Cd27","Gzmk","Vsir","Tnfrsf18",'Tnfrsf9','Tnfrsf4','Cd226','Eomes','Itgae','Alcam','Btla','Cd244','Entpd1','Havcr2','Ctla4','Lag3','Pdcd1','Tox','Tigit')
plot <- DotPlot(CD8T_mice, features = Exhaustion, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("exh_markers_mice_CD8.pdf", plot = plot, width = 8, height = 10, dpi = 300)


#FigS6C
subtypes <- c('Ccl3',"Ccl4","Ccl5","Klrd1","Ifng","Cd200r4",
              "S100a4","S100a6","S100a13","Lgals1","Lgals3",
              "Il7r","Tcf7","Ccr7","Ccr8","Sell","Cd83","Il2ra","Slamf6",
              "Ifit1","Ifit2","Ifit3","Ifi208","Ifit3b","Oas3","Oasl2","Isg15","Stat2","Rsad2",
              "Ung","Cdca7","Zcchc6","Aldoc","Nrgn","Mgea5","Nhp2l1","Sssca1","Fam26f","Fdx1l")
DefaultAssay(CD8T_mice) <- 'SCT'
avg_exp <- AverageExpression(CD8T_mice, assays = "SCT", features = subtypes, group.by = 'subtype', return.seurat = T)
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


#FigS6D
Idents(CD8T_mice) <- CD8T_mice$subtype
CD8_labels <- as.character(Idents(CD8T_mice))
CD8_labels[CD8_labels %in% c("C1", "C2", "C3", "C5", "C6")] <- "exh"
CD8T_mice$CD8_labels <- CD8_labels
Idents(CD8T_mice) <- CD8T_mice$CD8_labels
Idents(CD8T_mice) <- factor(Idents(CD8T_mice),levels = c("exh","Mki67",'others'))

mice_inte_exhaustedT <- AddModuleScore(
  object = mice_inte_exhaustedT,
  features = Exhaustion,
  name = 'exhausted_score_mice'
)

vln <- VlnPlot(mice_inte_exhaustedT, 'exhausted_score_mice1', pt.size = 0) + coord_cartesian(ylim = c(-0.8, 1))
vln + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "red")

exhaustion_data <- data.frame(
  score = CD8T_mice$exhausted_score_mice1,
  group = CD8T_mice$CD8_labels
)

exhaustion_data$group <- factor(exhaustion_data$group, levels = c("exh", "Mki67", "others"))

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
  c("exh", "Mki67"),
  c("exh", "others"), 
  c("Mki67", "others")
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