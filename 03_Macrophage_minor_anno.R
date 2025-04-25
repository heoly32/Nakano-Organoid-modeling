library(Seurat)
library(dplyr)
library(dittoSeq)
library(harmony)
library(ggplot2)

## Set directories
out_dir <- "Analysis/5K_2L_Merge_Macrophage/"
rds_file <- "Analysis/5K_2L_Merge/Merged_5K_2L.rds"
system(paste0("mkdir -p ", out_dir))

## Load Merged Seurat object with major cell type labels
merged <- readRDS(rds_file)
merged

## Extract Macrophage and mast cells only
seu <- subset(merged, `main.types` == "Macrophages" | `main.types` == "Mast")
seu
rm(merged)
gc()

## preprocessing
### SCTransform and PCA
seu <- SCTransform(seu,
                   assay = "RNA", 
                   new.assay.name = "SCT",
                   vst.flavor = "v2", 
                   verbose = T, 
                   vars.to.regress = c("percent.mt")) %>% 
  RunPCA(npcs = 50, verbose = FALSE)
### Harmony
seu <- RunHarmony(seu, group.by.vars = c("SampleID", "Groups", "Tissue"), 
                  reduction.use = "pca")
### ElbowPlot
ElbowPlot(seu, ndims = 50, reduction = "harmony")

### UMAP and TSNE
n_pcs <- 30
reduc_use <- "harmony"
seu <- RunUMAP(seu, reduction = reduc_use, dims = 1:n_pcs, verbose = FALSE) %>% 
  RunTSNE(reduction = reduc_use, dims = 1:n_pcs, verbose = FALSE)
## Neighbors and clustering
seu <- FindNeighbors(seu, reduction = reduc_use, dims = 1:n_pcs, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.4)

## Find all markers for each cluster
markers <- FindAllMarkers(seu, logfc.threshold = 0.25, min.pct = 0.1, only.pos = T)

markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) -> top10

## Heatmap with top10 markers 
DoHeatmap(seu, top10$gene)

seu$sub.types <- seu$seurat_clusters %>% as.character()

seu$sub.types[seu$seurat_clusters %in% c(0, 4, 8)] <- "NLRP3+ TAM"
seu$sub.types[seu$seurat_clusters %in% c(2, 3, 7)] <- "C1Q+ TAM"
seu$sub.types[seu$seurat_clusters %in% c(1)] <- "SPP1+ TAM"
seu$sub.types[seu$seurat_clusters %in% c(9)] <- "CXCL9+ TAM"
seu$sub.types[seu$seurat_clusters %in% c(5, 12)] <- "Monocyte"
seu$sub.types[seu$seurat_clusters %in% c(6, 10, 11)] <- "NLRP3+ TAM"

dim1 <- dittoDimPlot(seu, var = "sub.types", reduction.use = "umap", do.label = T, labels.size = 3)
dim1

saveRDS(seu, paste0(out_dir, "/Seurat_Obj_wSubtypes.rds"))


