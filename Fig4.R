library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

#Fig4A
Exhaustion <- list(c("KLF4","ELF4","PRDM1","ID2","JUN","TBX21","TCF7","CCR5","CST7","NKG7","CD28","CD27","GZMK","VSIR","TNFRSF18",'TNFRSF9','TNFRSF4','CD226','EOMES','ITGAE','ALCAM','BTLA','CD244','ENTPD1','HAVCR2','CTLA4','LAG3','TOX','TIGIT','PDCD1'))
Idents(human_inte_exhaustedT) <- human_inte_exhaustedT$humanC
Idents(human_inte_exhaustedT) <- factor(Idents(human_inte_exhaustedT),levels = c("C1","C2","C3","C4","C5"))

human_inte_exhaustedT <- AddModuleScore(
  object = human_inte_exhaustedT,
  features = Exhaustion,
  name = 'exhausted_score',
  seed = 123
)
vln <- VlnPlot(human_inte_exhaustedT, 'exhausted_score1', pt.size = 0) + coord_cartesian(ylim = c(-0.3, 1))
vln + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "red")

exhaustion_data <- data.frame(
  score = human_inte_exhaustedT$exhausted_score1,
  group = human_inte_exhaustedT$humanC
)

exhaustion_data$group <- factor(exhaustion_data$group, levels = c("C1","C2","C3","C4","C5"))

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
  c("C1", "C2"),
  c("C1", "C3"),
  c("C1", "C4"),
  c("C1", "C5"),
  c("C2", "C3"),
  c("C2", "C4"),
  c("C2", "C5"),
  c("C3", "C4"),
  c("C3", "C5"),
  c("C4", "C5")
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


#Fig4B
Exhaustion <- list(c("Klf4","Elf4","Prdm1","Id2","Jun","Tbx21","Tcf7","Ccr5","Cst7","Nkg7","Cd28","Cd27","Gzmk","Vsir","Tnfrsf18",'Tnfrsf9','Tnfrsf4','Cd226','Eomes','Itgae','Alcam','Btla','Cd244','Entpd1','Havcr2','Ctla4','Lag3','Tox','Tigit','Pdcd1'))
Idents(mice_inte_exhaustedT) <- mice_inte_exhaustedT$humanC
Idents(mice_inte_exhaustedT) <- factor(Idents(mice_inte_exhaustedT),levels = c("C1","C2","C3","C5","C6"))

mice_inte_exhaustedT <- AddModuleScore(
  object = mice_inte_exhaustedT,
  features = Exhaustion,
  name = 'exhausted_score_mice',
  seed = 123
)
vln <- VlnPlot(mice_inte_exhaustedT, 'exhausted_score_mice1', pt.size = 0) + coord_cartesian(ylim = c(-0.3, 1))
vln + stat_summary(fun = mean, geom = "crossbar", width = 0.5, color = "red")

exhaustion_data <- data.frame(
  score = mice_inte_exhaustedT$exhausted_score_mice1,
  group = mice_inte_exhaustedT$humanC
)

exhaustion_data$group <- factor(exhaustion_data$group, levels = c("C1","C2","C3","C5","C6"))

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
  c("C1", "C2"),
  c("C1", "C3"),
  c("C1", "C5"),
  c("C1", "C6"),
  c("C2", "C3"),
  c("C2", "C5"),
  c("C2", "C6"),
  c("C3", "C5"),
  c("C3", "C6"),
  c("C5", "C6")
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


#Fig4C
library(readxl)
genes <- read_xlsx('signature genes.xlsx')
colnames(genes) <- genes[2,]
genes <- genes[-c(1:2),]
genes <- as.vector(genes)
gene_vector <- as.character(unlist(genes))
gene_vector <- gene_vector[!is.na(gene_vector) & gene_vector != ""]

DefaultAssay(human_inte_exhaustedT) <- 'SCT'
Idents(human_inte_exhaustedT) <- human_inte_exhaustedT$humanC
avg_exp <- AverageExpression(human_inte_exhaustedT, assays = "SCT", features = gene_vector, group.by = 'humanC')
avg_exp <- avg_exp$SCT
avg_exp <- as.matrix(avg_exp)
scaled_exp <- t(scale(t(avg_exp)))
dist_matrix <- as.dist(1 - cor(avg_exp, method = "pearson"))
hclust_result <- hclust(dist_matrix, method = "ward.D2")
plot(hclust_result, 
     xlab = "Clusters",
     ylab = "Height",
     cex = 0.8)

similarity_matrix <- cor(scaled_exp, method = "pearson")

cor_with_pvalues <- function(mat) {
  n_clusters <- ncol(mat)
  correlation_matrix <- matrix(NA, nrow = n_clusters, ncol = n_clusters)
  pvalue_matrix <- matrix(NA, nrow = n_clusters, ncol = n_clusters)
  
  rownames(correlation_matrix) <- colnames(correlation_matrix) <- colnames(mat)
  rownames(pvalue_matrix) <- colnames(pvalue_matrix) <- colnames(mat)
  
  for(i in 1:n_clusters) {
    for(j in 1:n_clusters) {
      if(i == j) {
        correlation_matrix[i, j] <- 1
        pvalue_matrix[i, j] <- 0
      } else {
        cor_test <- cor.test(mat[, i], mat[, j], method = "pearson")
        correlation_matrix[i, j] <- cor_test$estimate
        pvalue_matrix[i, j] <- cor_test$p.value
      }
    }
  }
  
  return(list(correlation = correlation_matrix, pvalue = pvalue_matrix))
}

# Calculate correlations and p-values
cor_results <- cor_with_pvalues(scaled_exp)
similarity_matrix <- cor_results$correlation
pvalue_matrix <- cor_results$pvalue

# Function to get significance stars
get_significance_stars <- function(p_val) {
  if(p_val < 0.1) return("*")
  else return("")
}

significance_matrix <- matrix(
  sapply(pvalue_matrix, get_significance_stars),
  nrow = nrow(pvalue_matrix),
  ncol = ncol(pvalue_matrix)
)
rownames(significance_matrix) <- rownames(pvalue_matrix)
colnames(significance_matrix) <- colnames(pvalue_matrix)

similarity_heatmap_with_stars <- Heatmap(
  similarity_matrix,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = hclust_result,
  cluster_columns = hclust_result,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Cluster Similarity Matrix_Human",
  row_title = "Clusters",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", similarity_matrix[i, j]), 
              x, y + unit(1, "mm"), 
              gp = gpar(fontsize = 8, fontface = "bold"))
    
    if(significance_matrix[i, j] != "") {
      grid.text(significance_matrix[i, j], 
                x, y - unit(2, "mm"), 
                gp = gpar(fontsize = 10, col = "black", fontface = "bold"))
    }
  },
  heatmap_legend_param = list(
    title = "Pearson\nCorrelation",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(similarity_heatmap_with_stars)


#Fig4D
genes <- read_xlsx('signature genes.xlsx', sheet = 2)
colnames(genes) <- genes[2,]
genes <- genes[-c(1:2),]
genes <- as.vector(genes)
gene_vector <- as.character(unlist(genes))
gene_vector <- gene_vector[!is.na(gene_vector) & gene_vector != ""]


Idents(mice_inte_exhaustedT) <- mice_inte_exhaustedT$humanC
DefaultAssay(mice_inte_exhaustedT) <- 'SCT'
avg_exp <- AverageExpression(mice_inte_exhaustedT, assays = "SCT",  features = gene_vector, group.by = 'humanC')
avg_exp <- avg_exp$SCT
avg_exp <- as.matrix(avg_exp)
scaled_exp <- t(scale(t(avg_exp)))
dist_matrix <- as.dist(1 - cor(avg_exp, method = "pearson"))
hclust_result <- hclust(dist_matrix, method = "ward.D2")
plot(hclust_result, 
     xlab = "Clusters",
     ylab = "Height",
     cex = 0.8)

similarity_matrix <- cor(scaled_exp, method = "pearson")

cor_with_pvalues <- function(mat) {
  n_clusters <- ncol(mat)
  correlation_matrix <- matrix(NA, nrow = n_clusters, ncol = n_clusters)
  pvalue_matrix <- matrix(NA, nrow = n_clusters, ncol = n_clusters)
  
  rownames(correlation_matrix) <- colnames(correlation_matrix) <- colnames(mat)
  rownames(pvalue_matrix) <- colnames(pvalue_matrix) <- colnames(mat)
  
  for(i in 1:n_clusters) {
    for(j in 1:n_clusters) {
      if(i == j) {
        correlation_matrix[i, j] <- 1
        pvalue_matrix[i, j] <- 0
      } else {
        # Perform correlation test
        cor_test <- cor.test(mat[, i], mat[, j], method = "pearson")
        correlation_matrix[i, j] <- cor_test$estimate
        pvalue_matrix[i, j] <- cor_test$p.value
      }
    }
  }
  
  return(list(correlation = correlation_matrix, pvalue = pvalue_matrix))
}

# Calculate correlations and p-values
cor_results <- cor_with_pvalues(scaled_exp)
similarity_matrix <- cor_results$correlation
pvalue_matrix <- cor_results$pvalue

# Function to get significance stars
get_significance_stars <- function(p_val) {
  if(p_val < 0.05) return("*")
  else return("")
}

significance_matrix <- matrix(
  sapply(pvalue_matrix, get_significance_stars),
  nrow = nrow(pvalue_matrix),
  ncol = ncol(pvalue_matrix)
)
rownames(significance_matrix) <- rownames(pvalue_matrix)
colnames(significance_matrix) <- colnames(pvalue_matrix)

similarity_heatmap_with_stars <- Heatmap(
  similarity_matrix,
  name = "Correlation",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = hclust_result,
  cluster_columns = hclust_result,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "Cluster Similarity Matrix_Mice",
  row_title = "Clusters",
  cell_fun = function(j, i, x, y, width, height, fill) {

    grid.text(sprintf("%.2f", similarity_matrix[i, j]), 
              x, y + unit(1, "mm"), 
              gp = gpar(fontsize = 8, fontface = "bold"))
    
    if(significance_matrix[i, j] != "") {
      grid.text(significance_matrix[i, j], 
                x, y - unit(2, "mm"), 
                gp = gpar(fontsize = 10, col = "black", fontface = "bold"))
    }
  },
  heatmap_legend_param = list(
    title = "Pearson\nCorrelation",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(similarity_heatmap_with_stars)


