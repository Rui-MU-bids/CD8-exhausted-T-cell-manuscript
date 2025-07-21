#Fig6A
ac_exT <- subset(LCMV_ac_T, subset = Cd3e >0 & Cd8a > 0 & Cd4 == 0 & Pdcd1 >0 & Mki67 == 0 & Aif1 == 0)
ac_exT <- SCTransform(ac_exT,vst.flavor = 'v2')
ac_exT <- RunPCA(ac_exT)
ac_exT <- FindNeighbors(ac_exT,  dims = 1:30, verbose = FALSE, reduction = 'pca')
ac_exT <- FindClusters(ac_exT, resolution = 0.5)
ac_exT <- RunUMAP(ac_exT, dims = 1:30, verbose = FALSE, reduction = 'pca') 

Idents(ac_exT) <- ac_exT$SCT_snn_res.0.5
rename<- c("6","7")
ac_exT@meta.data$humanC<- ifelse(ac_exT@meta.data$SCT_snn_res.0.5 %in% rename, "C1", NA)
rename<- c("1","2","9")
ac_exT@meta.data$humanC<- ifelse(ac_exT@meta.data$SCT_snn_res.0.5%in% rename, "C2", ac_exT@meta.data$humanC)
rename<- c("4","5")
ac_exT@meta.data$humanC<- ifelse(ac_exT@meta.data$SCT_snn_res.0.5%in% rename, "C3", ac_exT@meta.data$humanC)
rename<- c("3")
ac_exT@meta.data$humanC<- ifelse(ac_exT@meta.data$SCT_snn_res.0.5%in% rename, "C4", ac_exT@meta.data$humanC)
rename<- c("0")
ac_exT@meta.data$humanC<- ifelse(ac_exT@meta.data$SCT_snn_res.0.5%in% rename, "C6", ac_exT@meta.data$humanC)
rename <- c("8")
ac_exT@meta.data$humanC<- ifelse(ac_exT@meta.data$SCT_snn_res.0.5%in% rename, "unknown", ac_exT@meta.data$humanC)
table(ac_exT$type, ac_exT$humanC)
ggplot(data, aes(fill = Cluster, y = Count, x = type)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2","#C26841")) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



#Fig6B-D
LCMV_exT <- subset(LCMV_T, subset = Cd3e >0 & Cd8a > 0 & Pdcd1 >0 & Mki67 == 0) #Cd4 all 0
LCMV_exT <- SCTransform(LCMV_exT)
LCMV_exT <- RunPCA(LCMV_exT)
Inte_LCMV_exT <- IntegrateLayers(
  object = LCMV_exT, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
Inte_LCMV_exT <- FindNeighbors(Inte_LCMV_exT,  dims = 1:30, verbose = FALSE, reduction = 'harmony')
Inte_LCMV_exT <- FindClusters(Inte_LCMV_exT, resolution = 0.8)
Inte_LCMV_exT <- RunUMAP(Inte_LCMV_exT, dims = 1:30, verbose = FALSE,reduction = 'harmony') 
Inte_LCMV_exT <- PrepSCTFindMarkers(Inte_LCMV_exT)
Inte_LCMV_exT.markers <- FindAllMarkers(Inte_LCMV_exT, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
Inte_LCMV_exT.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write_xlsx(Inte_LCMV_exT.markers,"label_LCMV_exT_clusters.xlsx")

Inte_LCMV_exT <- subset(Inte_LCMV_exT, idents = c("0","1","2","3","4","6","7","8","9","10"))
rename<- c("1","4")
Inte_LCMV_exT@meta.data$humanC<- ifelse(Inte_LCMV_exT@meta.data$SCT_snn_res.0.8 %in% rename, "C1", NA)
rename<- c("0","3","9")
Inte_LCMV_exT@meta.data$humanC<- ifelse(Inte_LCMV_exT@meta.data$SCT_snn_res.0.8%in% rename, "C2", Inte_LCMV_exT@meta.data$humanC)
rename<- c("2")
Inte_LCMV_exT@meta.data$humanC<- ifelse(Inte_LCMV_exT@meta.data$SCT_snn_res.0.8%in% rename, "C3", Inte_LCMV_exT@meta.data$humanC)
rename<- c("6","10")
Inte_LCMV_exT@meta.data$humanC<- ifelse(Inte_LCMV_exT@meta.data$SCT_snn_res.0.8%in% rename, "C4", Inte_LCMV_exT@meta.data$humanC)
rename<- c("7")
Inte_LCMV_exT@meta.data$humanC<- ifelse(Inte_LCMV_exT@meta.data$SCT_snn_res.0.8%in% rename, "C5", Inte_LCMV_exT@meta.data$humanC)
rename<- c("8")
Inte_LCMV_exT@meta.data$humanC<- ifelse(Inte_LCMV_exT@meta.data$SCT_snn_res.0.8%in% rename, "C6", Inte_LCMV_exT@meta.data$humanC)
table(Inte_LCMV_exT$humanC, Inte_LCMV_exT$orig.ident)

Idents(Inte_LCMV_exT) <- Inte_LCMV_exT$humanC
Idents(Inte_LCMV_exT) <- Inte_LCMV_exT$orig.ident

LCMV_exT <- subset(Inte_LCMV_exT, idents = "CTL")
Idents(LCMV_exT) <- LCMV_exT$humanC
Idents(LCMV_exT) <- factor(Idents(LCMV_exT),levels = c("C1","C2","C3","C4","C5","C6"))
DimPlot(
  LCMV_exT,
  combine = T, 
  label = F,
  cols = c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#C26841")
)


PD_exT <- subset(Inte_LCMV_exT, idents = "PD1")
Idents(PD_exT) <- PD_exT$humanC
Idents(PD_exT) <- factor(Idents(PD_exT),levels = c("C1","C2","C3","C4","C5","C6"))
DimPlot(
  PD_exT,
  combine = T, 
  label = F,
  cols = c("#E69F00", "#56B4E9","#009E73", "#F0E442", "#0072B2","#C26841")
)

ggplot(data, aes(fill = Cluster, y = Count, x = sample)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442", "#0072B2","#C26841")) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




#Fig6E
exT <- merge(exT_CTL, c(exT_PD, aT_exT,aPT_exT))
exT <- SCTransform(exT)
exT <- RunPCA(exT)
exT <- IntegrateLayers(
  object = exT, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
exT <- FindNeighbors(exT,dims = 1:30, verbose = FALSE, reduction = 'harmony')
exT <- FindClusters(exT,resolution = 0.6)
FMC.exT.anchors <- FindTransferAnchors(reference = FMC_exT, query = exT, dims = 1:30,
                                       reference.reduction = "pca",normalization.method = 'SCT')
predictions <- TransferData(anchorset = FMC.exT.anchors, refdata = FMC_exT$humanC,dims = 1:30)
exT <- AddMetaData(exT, metadata = predictions)

ggplot(data, aes(fill = Cluster, y = Count, x = sample)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#0072B2","#C26841")) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#Fig6F
CRC_qc <- SCTransform(CRC_qc)
CRC_qc <- RunPCA(CRC_qc)
CRC_qc <- IntegrateLayers(
  object = CRC_qc, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
CRC_qc <- FindNeighbors(CRC_qc, dims = 1:30, reduction = 'harmony')
CRC_qc <- FindClusters(CRC_qc, resolution = 0.5)
Tcells <- subset(CRC_qc, idents = c('0','1','3','4','7','11','14','15','17'))
exh <- subset(Tcells, subset = CD3E > 0 & CD8A > 0 & PDCD1 > 0 & MKI67 == 0 & CD4 == 0, slot = 'counts')
exh <- FindNeighbors(exh, dims = 1:30, reduction = 'harmony')
exh <- FindClusters(exh, resolution = 0.5)
exh <- RunUMAP(exh, dims = 1:30,reduction = 'harmony')
exh.markers <- FindAllMarkers(exh, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1, recorrect_umi=FALSE)
exh.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write_xlsx(exh.markers,"CRC_exh_clusters_SCT_0.5.xlsx")

Idents(exh) <- exh$SCT_snn_res.0.5
rename<- c('4','5')
exh@meta.data$humanC<- ifelse(exh@meta.data$SCT_snn_res.0.5 %in% rename, "C1", NA)
rename<- c("0","2")
exh@meta.data$humanC <- ifelse(exh@meta.data$SCT_snn_res.0.5 %in% rename, "C2", exh@meta.data$humanC)
rename<- c("1","3")
exh@meta.data$humanC <- ifelse(exh@meta.data$SCT_snn_res.0.5 %in% rename, "C3", exh@meta.data$humanC)
rename<- c("7")
exh@meta.data$humanC <- ifelse(exh@meta.data$SCT_snn_res.0.5 %in% rename, "C4", exh@meta.data$humanC)

table(exh_filter$humanC, exh_filter$Treatment)
ggplot(data, aes(fill = Cluster, y = Count, x = sample)) + 
  geom_bar(position = "fill", stat = "identity") + 
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73","#F0E442")) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )