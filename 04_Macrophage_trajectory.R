suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(dittoSeq))
suppressMessages(library(harmony))
suppressMessages(library(SeuratWrappers))

suppressMessages(library(monocle3))
library(patchwork)

## Set directories
out_dir <- "Analysis/5K_2L_Merge_Macro_Trajectory/"

system(paste0("mkdir -p ", out_dir))

suppressMessages(source("/Users/lyongheo/Projects/RScripts/Custom_Vis_scRNAseq.R"))

## Load Seurat object
rds_file <- "Analysis/5K_2L_Merge_Macrophage/Seurat_Obj_wSubtypes.rds" 

seu <- readRDS(rds_file)

## Exclude Fresh samples
t_seu <- subset(seu, Treatment != "Fresh")
t_seu$ident <- t_seu$sub.types

## Seurat object to cell data set
cds <- as.cell_data_set(t_seu)
cds <- cluster_cells(cds)

## Set minor types to cluster names for trajectory
celltypes <- cds$sub.types %>% as.character
names(celltypes) <- rownames(cds@colData)
cds@clusters[['UMAP']]$clusters <- celltypes

## Learn graph
cds <- learn_graph(cds, )

### Order cells from Monocyte
cds <- order_cells(cds, root_cells = rownames(cds@colData)[cds@colData$sub.types == "Monocyte"][1])
p1_1 <- plot_cells(cds, graph_label_size = 1, 
                   alpha = 0.2, 
                   trajectory_graph_segment_size = 1.5, cell_size = 0.6, 
                   color_cells_by = "pseudotime", 
                   label_cell_groups = FALSE, label_leaves = F, 
                   label_branch_points = FALSE) + ggtitle("From Monocyte")


### Order cells from C1Q+ TAM
cds <- order_cells(cds, 
                   root_cells = rownames(cds@colData)[cds@colData$sub.types == "C1Q+ TAM"][1])
p2_1 <- plot_cells(cds, graph_label_size = 1, 
                   alpha = 0.2, 
                   trajectory_graph_segment_size = 1.5, cell_size = 0.6, 
                   color_cells_by = "pseudotime", 
                   label_cell_groups = FALSE, label_leaves = F, 
                   label_branch_points = FALSE) + ggtitle("From C1Q+ TAM")

### Add pseudotime from Monocyte and C1Q+ TAM to metadata in Seurat object
cds <- order_cells(cds, root_cells = rownames(cds@colData)[cds@colData$sub.types == "Monocyte"][1])

t_seu <- AddMetaData(
  object = t_seu, 
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Mono_PT"
)

cds <- order_cells(cds, 
                   root_cells = rownames(cds@colData)[cds@colData$sub.types == "C1Q+ TAM"][1])

t_seu <- AddMetaData(
  object = t_seu, 
  metadata = cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "C1Q_PT"
)