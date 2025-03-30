#Fig2A
subtypes <- c('CCL3',"CCL4","KLRD1","EGR2","IFNG",
              "COL1A1","COL6A1","CCN1","CCN2",
              "IL7R","TCF7","CCR7","SELL",
              "HSPA1A","HSPA1B","HSPA6","DNAJA1","DNAJB1",
              "IFIT1","IFIT2","IFIT3","OAS1","OAS3","ISG15")
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

#Fig2B
Exhaustion <- c("ENTPD1","HAVCR2","CTLA4","LAG3","PDCD1","TOX","TIGIT")
plot <- DotPlot(exhaustedT, features = Exhaustion, dot.scale = 8, scale = T, cols = c("lightgrey", "#bd7fef"))+coord_flip()
ggsave("exh_markers.pdf", plot = plot, width = 5, height = 2.5, dpi = 300)

#Fig2C
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
