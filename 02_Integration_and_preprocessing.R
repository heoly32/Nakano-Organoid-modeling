library(Seurat)
library(dplyr)
library(dittoSeq)
library(harmony)
library(ggplot2)

## Set directories
out_dir <- "Analysis/5K_2L_Merge/"

system(paste0("mkdir -p ", out_dir))

source("Custom_Vis_scRNAseq.R")

## Load Seurat Objects
rds_files <- list(`RCC_2` = "Analysis/RCC_2/Merged/SeuratObj_Merged.rds", 
                  `RCC_12` = "Analysis/RCC_12/Merged/Seurat_Obj_all.rds", 
                  `RCC_25` = "Analysis/RCC_25/Merged/Seurat_Obj_all.rds",
                  `RCC_20` = "Analysis/RCC_20/Merged/Seurat_Obj_all.rds", 
                  `RCC_23` = "Analysis/RCC_23/Merged_Organoid/Seurat_Obj_all.rds", 
                  `NSCLC_24` = "Analysis/NSCLC_24/Merged/SeuratObj_Merged.rds",
                  `NSCLC_16` = "Analysis/36656LK/NSCLC_16/Seurat_Obj_all.rds")

seu.list <- lapply(rds_files, readRDS)

## Add SampleID columns before integration
for(n in names(seu.list)) {
  seu.list[[n]]$SampleID <- n
}

seu.list <- lapply(seu.list, 
                   function(x) {
                     x$Groups <- stringr::str_split(colnames(x), "_", simplify = T)[, 2] %>% as.character
                     x$Groups <- ifelse(x$Groups == "Fresh", "Fresh", "Organoid") %>% as.character
                     x
                   })

## Merge all metadata
metadata <- do.call(rbind, 
                    lapply(names(seu.list), function(x) {
                      t_seu <- seu.list[[x]]
                      data.frame(row.names = colnames(t_seu), check.names = F,
                                 Treatment = t_seu$orig.ident, ## Fresh, Control, and CD47
                                 SampleID = t_seu$SampleID, ## Each sample ID
                                 Groups = t_seu$Groups, ## Fresh | Organoid
                                 ID_Treatment = paste0(t_seu$SampleID, "_", t_seu$orig.ident))
                    }))


## Find co-exist genes over all samples
gene_names <- unname(do.call(c, lapply(seu.list, function(x) { rownames(x) })))
coexist_genes <- names(table(gene_names))[table(gene_names) == 7]
coexist_genes %>% length


## Create merged Seurat object
merged <- CreateSeuratObject(counts = cbind(seu.list[[1]]@assays$RNA$counts[coexist_genes, ],
                                            seu.list[[2]]@assays$RNA$counts[coexist_genes, ], 
                                            seu.list[[3]]@assays$RNA$counts[coexist_genes, ], 
                                            seu.list[[4]]@assays$RNA$counts[coexist_genes, ], 
                                            seu.list[[5]]@assays$RNA$counts[coexist_genes, ], 
                                            seu.list[[6]]@assays$RNA$counts[coexist_genes, ], 
                                            seu.list[[7]]@assays$RNA$counts[coexist_genes, ]), 
                             project = "5K_2L", min.features = 200)

## Add Metadata
merged <- AddMetaData(merged, metadata)

## Add extra column "Tissue"
merged$Tissue <- "Kidney"
merged$Tissue[grepl("^NSCLC", merged$SampleID)] <- "Lung"

table(merged$Tissue, merged$Treatment)

## Add percent.mt values
merged[['percent.mt']] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged

# Pre-processing
## SCTransform and PCA
merged <- SCTransform(merged,
                      assay = "RNA", 
                      new.assay.name = "SCT",
                      vst.flavor = "v2", 
                      verbose = T, 
                      vars.to.regress = c("percent.mt")) %>% 
  RunPCA(npcs = 50, verbose = FALSE)

## PCA plots 
## Colored by SampleID, Treatment, Groups, and Tissue
options(repr.plot.width=8, repr.plot.height=6)
dittoDimPlot(merged, var = "SampleID", reduction.use = "pca")
dittoDimPlot(merged, var = "Treatment", reduction.use = "pca")
dittoDimPlot(merged, var = "Groups", reduction.use = "pca")
dittoDimPlot(merged, var = "Tissue", reduction.use = "pca")

options(repr.plot.width=12, repr.plot.height=12)
dittoDimPlot(merged, var = "SampleID", reduction.use = "pca", split.by = "SampleID")
dittoDimPlot(merged, var = "Treatment", reduction.use = "pca", split.by = "Treatment")
dittoDimPlot(merged, var = "Groups", reduction.use = "pca", split.by = "Groups")
dittoDimPlot(merged, var = "Tissue", reduction.use = "pca", split.by = "Tissue")

## Harmony
### Set SampleID, Groups, and Tissue as the names of covariates that hormony will remove its effect on the data
merged <- RunHarmony(merged, group.by.vars = c("SampleID", "Groups", "Tissue"), reduction.use = "pca")

## Harmony plots
options(repr.plot.width=8, repr.plot.height=6)
dittoDimPlot(merged, var = "SampleID", reduction.use = "harmony")
dittoDimPlot(merged, var = "Treatment", reduction.use = "harmony")
dittoDimPlot(merged, var = "Groups", reduction.use = "harmony")
dittoDimPlot(merged, var = "Tissue", reduction.use = "harmony")

options(repr.plot.width=12, repr.plot.height=12)
dittoDimPlot(merged, var = "SampleID", reduction.use = "harmony", split.by = "SampleID")
dittoDimPlot(merged, var = "Treatment", reduction.use = "harmony", split.by = "Treatment")
dittoDimPlot(merged, var = "Groups", reduction.use = "harmony", split.by = "Groups")
dittoDimPlot(merged, var = "Tissue", reduction.use = "harmony", split.by = "Tissue")

## ElbowPlot on Harmony
options(repr.plot.width=8, repr.plot.height=6)

ElbowPlot(merged, ndims = 50, reduction = "harmony")


## Set 20, as number of PCs
n_pcs <- 20
reduc_use <- "harmony"

# Run UMAP and TSNE
merged <- RunUMAP(merged, reduction = reduc_use, dims = 1:n_pcs, verbose = FALSE) %>% 
  RunTSNE(reduction = reduc_use, dims = 1:n_pcs, verbose = FALSE)

# Run FindNeighbors and Clustering
merged <- FindNeighbors(merged, reduction = reduc_use, dims = 1:n_pcs, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.4)

## UMAP and TSNE plots with cluster label
umap_and_tsne <- dittoDimPlot(merged, var = "seurat_clusters", reduction.use = "umap", do.label = T, legend.show = F) + 
  dittoDimPlot(merged, var = "seurat_clusters", reduction.use = "tsne", do.label = T, legend.show = F)
umap1 <- dittoDimPlot(merged, var = "seurat_clusters", split.by = "SampleID", reduction.use = "umap", do.label = T, legend.show = F)
tsne1 <- dittoDimPlot(merged, var = "seurat_clusters", split.by = "SampleID", reduction.use = "tsne", do.label = T, legend.show = F)
umap2 <- dittoDimPlot(merged, var = "seurat_clusters", split.by = "Groups", reduction.use = "umap", do.label = T, legend.show = F)
tsne2 <- dittoDimPlot(merged, var = "seurat_clusters", split.by = "Groups", reduction.use = "tsne", do.label = T, legend.show = F)

umap_and_tsne
umap1
tsne1
umap2
tsne2


## Feature plots with major markers for Immune cells, CD3 T, Macrophages, Mast, CD4 T, CD8 T, B, MAST, Endothelial, and Fibroblast
options(repr.plot.width=12, repr.plot.height=10)
used_markers <- c("PTPRC", "CD3D", "FTH1", "LYZ", "CD68", "NKG7", "CD4", "CD8A",
                  "CD79A", "KIT", "PECAM1", "COL1A1")

dim_markers <- vis_dimplot(merged, used_markers)
dot_markers1 <- dittoDotPlot(merged, used_markers, group.by = "seurat_clusters")
dot_markers2 <- dittoDotPlot(merged, used_markers, group.by = "seurat_clusters", split.by = "orig.ident")

dim_markers
dot_markers1
dot_markers2

## Find all markers for each clusters
markers <- FindAllMarkers(merged, logfc.threshold = 0.25, min.pct = 0.1, only.pos = T, verbose = F)

markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) -> top5

write.table(x = markers, file = paste0(out_dir, "Table_DE_byClusters.txt"), quote = F, sep = "\t", row.names = F, col.names = T)

## Cell type annotation
### Unknown, Endothelial and Fibroblasts will be excluded in the future analysis
merged$main.types <- merged$seurat_clusters %>% as.character()

merged$main.types.info[merged$`SCT_snn_res.0.4` %in% c(1, 3, 4, 7, 9, 13, 16, 17)] <- "Macrophage"
merged$main.types[merged$main.types %in% c(0, 8, 15, 20, 21)] <- "CD8 T"
merged$main.types[merged$main.types %in% c(2)] <- "CD4 T"
merged$main.types[merged$main.types %in% c(6)] <- "NK"
merged$main.types[merged$main.types %in% c(9)] <- "Monocytes"
merged$main.types[merged$main.types %in% c(10)] <- "B"
merged$main.types[merged$main.types %in% c(14)] <- "Mast"
merged$main.types[merged$main.types %in% c(19)] <- "DC"
merged$main.types[merged$main.types %in% c(5, 11, 12, 18)] <- "Exclude"

merged <- subset(merged, main.types != "Exclude")

dim1 <- dittoDimPlot(merged, "main.types", do.label = T, labels.size = 3, legend.show = F)
dim2 <- dittoDimPlot(merged, "main.types.info", do.label = T, labels.size = 3, legend.show = F)

dim1
dim2
## Heatmap with mean expression on major markers
heat_main <- vis_mean_heatmap(merged, used_markers, mean_by = "main.types")
heat_main

## Bar plot - Cell type distribution
bar1 <- dittoBarPlot(merged, var = "main.types", group.by = "Treatment", split.by = "SampleID")
bar1

## Save merged Seurat object to RDS
saveRDS(merged, paste0(out_dir, "/Merged_5K_2L.rds"))



