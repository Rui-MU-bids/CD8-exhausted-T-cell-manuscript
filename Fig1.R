#Fig1F & G
#run by command line: Rscript
options(future.globals.maxSize = 8000 * 1024^2) 
suppressPackageStartupMessages({
  library(Seurat)
  library(pheatmap)
  library(reshape2)
  library(ggplot2)
  library(dplyr)
})


merged <- readRDS('merged.RDS')
human_inte_exhaustedT <- readRDS('human_inte_exhaustedT.RDS')

Idents(merged) <- merged$meta.cluster
Tex <- subset(merged, idents = c('CD8.c11.Tex.PDCD1','CD8.c12.Tex.CXCL13','CD8.c13.Tex.myl12a','CD8.c14.Tex.TCF7'))
Tex <- SCTransform(Tex)


# Get count matrices
Tex_counts <- Tex@assays$SCT$counts
human_counts <- human_inte_exhaustedT@assays$SCT$counts

# Find common genes
common_genes <- intersect(rownames(Tex_counts), rownames(human_counts))

# Subset to common genes
Tex_counts <- Tex_counts[common_genes, ]
CD8T_subset <- human_counts[common_genes, ]

# Calculate correlation matrix
#cat('start calculating correlaiton')
start_time <- Sys.time()
cor_matrix <- cor(as.matrix(Tex_counts), as.matrix(CD8T_subset))


pdf("correlation_heatmap.pdf", width = 10, height = 8)
pheatmap(cor_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         treeheight_row = 0,
         treeheight_col = 0,
         main = "Cell-Cell Correlation")
dev.off()

end_time <- Sys.time()
runtime <- end_time - start_time
cat("Completed in", round(runtime, 2), attr(runtime, "units"), "\n")

# Save correlation matrix
saveRDS(cor_matrix, "cor_matrix.RDS")


# Calculate average expression per cell
Tex_avg <- rowMeans(Tex_counts)
human_avg <- rowMeans(CD8T_subset)

# Perform correlation test
overall_cor_test <- cor.test(Tex_avg, human_avg, method = "pearson")

# Create data frame for plotting
plot_data <- data.frame(
  Tex_expression = Tex_avg,
  Human_expression = human_avg,
  Gene = names(Tex_avg)
)


p <- ggplot(plot_data, aes(x = Tex_expression, y = Human_expression)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 1.2) +
  labs(
    x = "Average Expression in Zheng et al Tex Cells",
    y = "Average Expression in Human Exhausted T Cells",
    title = "Gene Expression Correlation: Tex vs Human Exhausted T Cells",
    subtitle = paste0("r = ", round(overall_cor_test$estimate, 3), 
                      ", p = ", format(overall_cor_test$p.value, scientific = TRUE))
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )






#Fig1E
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


