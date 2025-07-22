options(future.globals.maxSize = 8000 * 1024^2)
library(Seurat)
library(ggplot2)
CRC_Tcells <- CreateSeuratObject(CRC_Tcells_RNAcount)
Lung_Tcells <- CreateSeuratObject(Lung_Tcells_RNAcount)
Melanoma_Tcells <- CreateSeuratObject(Melanoma_Tcells_RNAcount)
MPNST_Tcells <- CreateSeuratObject(MPNST_Tcells_RNAcount)
Ovarian_Tcells <- CreateSeuratObject(Ovarian_Tcells_RNAcount)
PA_Tcells <- CreateSeuratObject(PA_Tcells_RNAcount)
Pancreatic_Tcells <- CreateSeuratObject(Pancreatic_Tcells_RNAcount)
Thyroid_Tcells <- CreateSeuratObject(Thyroid_Tcells_RNAcount)

Tcells <- merge(CRC_Tcells, 
                c(Lung_Tcells,Melanoma_Tcells,MPNST_Tcells,Ovarian_Tcells,PA_Tcells,Pancreatic_Tcells,Thyroid_Tcells),
                add.cell.ids = c('CRC','Lung','Melanoma','MPNST','Ovarian','PA','Pancreatic','Thyroid'))

Tcells$orig.idents <- sapply(strsplit(colnames(Tcells), "_"), `[`, 1)
table(Tcells$orig.idents)
Tcells <- JoinLayers(Tcells)
Tcells[['RNA']] <- split(Tcells[['RNA']], f = Tcells$orig.idents)

Idents(Tcells) <- Tcells$orig.idents
Tcells[["percent.mt"]] <- PercentageFeatureSet(Tcells, pattern = "MT-")
Tcells[["percent.rib"]] <- PercentageFeatureSet(Tcells, pattern = "RP-")
VlnPlot(Tcells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"),ncol=4)

Tcells <- SCTransform(Tcells, vars.to.regress = c('percent.mt','percent.rib'))
Tcells  <- RunPCA(Tcells)

Tcells <-  IntegrateLayers(
  object = Tcells, method = HarmonyIntegration,normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

Tcells <- FindNeighbors(Tcells, dims = 1:30, reduction = 'harmony')
Tcells <- FindClusters(Tcells, resolution = 0.2)

Tcells <- RunUMAP(Tcells, dims = 1:30, reduction = 'harmony')
DimPlot(Tcells, label = TRUE)

saveRDS(Tcells, 'Tcells_humancancer.RDS')

exhaustedT <- subset(Tcells, subset = CD3E > 0 & CD8A > 0 & PDCD1 > 0 & MKI67 == 0 & CD4 == 0, slot = 'counts')
saveRDS(exhaustedT, 'exhaustedT.RDS')





