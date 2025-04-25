## Sample QC 
library(Seurat)
library(dplyr)

out_dir <- "Analysis/NSCLC_24/Control"

system(paste0("mkdir -p ", out_dir))

suppressMessages(source("./Custom_Vis_scRNAseq.R"))

## Load raw count in h5 file
raw_count_h5 <- "Data/NSCLC_24_Control/count/sample_filtered_feature_bc_matrix.h5"

cnt <- Read10X_h5(raw_count_h5)
seu <- CreateSeuratObject(cnt[['Gene Expression']], project = "Control", min.features = 300)
seu[['percent.mt']] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu


## General QC (Before)
Idents(seu) <- seu$orig.ident
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 

## Filter out low-qual cells
## Need to define different number sample by sample
seu <- subset(seu, subset = nFeature_RNA > 300 & 
                nFeature_RNA < 5000 & 
                nCount_RNA < 20000 &
                percent.mt < 5 
)
seu

## General QC (After)
Idents(seu) <- seu$orig.ident
options(repr.plot.width=14, repr.plot.height=6)
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

options(repr.plot.width=10, repr.plot.height=6)
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

## Save Seurat object to RDS file
saveRDS(seu, paste0(out_dir, "/SeuratObj_AfterQC.rds"))


