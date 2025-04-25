# Organoid modeling of functional tumor-associated macrophages reveals phagocytosis-induced immune cell dynamics
## Workflow - Seurat v5, Monocle 3 and CellChat v2
 1. Sample QC
 2. Integration and preprocessing
   2-1. Pre-processing: SCTransform v2, and Principal component analysis
   2-2. Integration of dataset using Harmony v2
   2-3. Dimensionality reduction: Uniform manifold approximation and projection (UMAP)
   2-4. Clustering using Louvain algorithm
   2-5. Manual cell type annotation based on well-know marker genes
 3. Identification of minor types of macrophages (TAMs, and Monocyte)
 4. Trajectory inference analysis using Monocle 3 
 5. Cell-cell interaction based on ligand-receptor database using CellChat v2
