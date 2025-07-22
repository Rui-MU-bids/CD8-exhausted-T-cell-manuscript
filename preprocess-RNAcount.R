library(Seurat)

###extract dhg and Till PA RNA count matrix
PA_Tcells <- readRDS('PA1/PA_T.RDS')
table(PA_Tcells$orig.ident)
DefaultAssay(PA_Tcells) <- 'RNA'
PA_Tcells[['RNA']] <- JoinLayers(PA_Tcells[['RNA']])
PA_Tcells_RNAcount <- PA_Tcells[["RNA"]]$counts


###extract lung cancer RNA count matrix
Lung_Tcells <- readRDS('Lung cancer/Lung_Tcells.RDS')
table(Lung_Tcells$orig.ident)
DefaultAssay(Lung_Tcells) <- 'RNA'
Lung_Tcells[['RNA']] <- JoinLayers(Lung_Tcells[['RNA']])
Lung_Tcells_RNAcount <- Lung_Tcells[['RNA']]$counts


###extract MPNST cancer RNA count matrix
MPNST_Tcells <- readRDS('MPNST/MPNST_T.RDS')
table(MPNST_Tcells$label)
DefaultAssay(MPNST_Tcells) <- 'RNA'
MPNST_Tcells[['RNA']] <- JoinLayers(MPNST_Tcells[['RNA']])
MPNST_Tcells_RNAcount <- MPNST_Tcells[['RNA']]$counts


###extract Thyroid cancer RNA count matrix
Thyroid_Tcells1 <- readRDS('Thyroid cancer/Tcells_PTC.RDS')
Thyroid_Tcells2 <- readRDS('Thyroid cancer/Tcells_ATC.RDS')
table(Thyroid_Tcells1$orig.ident)
table(Thyroid_Tcells2$orig.ident)
DefaultAssay(Thyroid_Tcells1) <- 'RNA'
DefaultAssay(Thyroid_Tcells2) <- 'RNA'
Thyroid_Tcells1[['RNA']] <- JoinLayers(Thyroid_Tcells1[['RNA']])
Thyroid_Tcells2[['RNA']] <- JoinLayers(Thyroid_Tcells2[['RNA']])
Thyroid_Tcells1[['RNA']] <- split(Thyroid_Tcells1[['RNA']], f = Thyroid_Tcells1$orig.ident)
Thyroid_Tcells2[['RNA']] <- split(Thyroid_Tcells2[['RNA']], f = Thyroid_Tcells2$orig.ident)
Thyroid_Tcells <- merge(Thyroid_Tcells1, Thyroid_Tcells2)
Thyroid_Tcells[['RNA']] <- JoinLayers(Thyroid_Tcells[['RNA']])
Thyroid_Tcells_RNAcount <- Thyroid_Tcells[['RNA']]$counts

###extract Pancreatic cancer RNA count matrix
Pancreatic_Tcells <- readRDS('Pancreatic cancer/Tcells_Panc.RDS')
DefaultAssay(Pancreatic_Tcells) <- 'RNA'
Pancreatic_Tcells[['RNA']] <- JoinLayers(Pancreatic_Tcells[['RNA']])
Pancreatic_Tcells_RNAcount <- Pancreatic_Tcells[['RNA']]$counts


###extract Ovarian cancer RNA count matrix
Ovarian_Tcells <- readRDS('Ovarian/Tcells_OC.RDS')
DefaultAssay(Ovarian_Tcells) <- 'RNA'
Ovarian_Tcells[['RNA']] <- JoinLayers(Ovarian_Tcells[['RNA']])
Ovarian_Tcells_RNAcount <- Ovarian_Tcells[['RNA']]$counts


###extract Melanoma cancer RNA count matrix
Melanoma_Tcells <- readRDS('Melanoma/Mela_T.RDS')
DefaultAssay(Melanoma_Tcells) <- 'RNA'
Melanoma_Tcells[['RNA']] <- split(Melanoma_Tcells[['RNA']], f = Melanoma_Tcells$orig.ident)
Melanoma_Tcells[['RNA']] <- JoinLayers(Melanoma_Tcells[['RNA']])
Melanoma_Tcells_RNAcount <- Melanoma_Tcells[['RNA']]$counts


###extract CRC cancer RNA count matrix
CRC_Tcells <- readRDS('CRC/cc_T.RDS')
DefaultAssay(CRC_Tcells) <- 'RNA'
table(CRC_Tcells$samples)
CRC_Tcells[['RNA']] <- split(CRC_Tcells[['RNA']], f = CRC_Tcells$samples)
CRC_Tcells[['RNA']] <- JoinLayers(CRC_Tcells[['RNA']])
CRC_Tcells_RNAcount <- CRC_Tcells[['RNA']]$counts

saveRDS(PA_Tcells_RNAcount, 'PA_Tcells_RNAcount.RDS')
saveRDS(Lung_Tcells_RNAcount, 'Lung_Tcells_RNAcount.RDS')
saveRDS(MPNST_Tcells_RNAcount, 'MPNST_Tcells_RNAcount.RDS')
saveRDS(Ovarian_Tcells_RNAcount, 'Ovarian_Tcells_RNAcount.RDS')
saveRDS(Pancreatic_Tcells_RNAcount, 'Pancreatic_Tcells_RNAcount.RDS')
saveRDS(Thyroid_Tcells_RNAcount, 'Thyroid_Tcells_RNAcount.RDS')
saveRDS(Melanoma_Tcells_RNAcount, 'Melanoma_Tcells_RNAcount.RDS')
saveRDS(CRC_Tcells_RNAcount, 'CRC_Tcells_RNAcount.RDS')