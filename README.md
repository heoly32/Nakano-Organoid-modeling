# Organoid modeling of functional tumor-associated macrophages reveals phagocytosis-induced immune cell dynamics
## Workflow - Seurat v5, Monocle 3 and CellChat v2
 1. QC
 2. Pre-processing: SCTransform v2, and Principal component analysis
 3. Integration using Harnomy v2
 4. Dimensionality reduction: Uniform manifold approximation and projection (UMAP)
 5. Clustering using Louvain algorithm
 6. Manual cell type annotation based on well-know marker genes
 7. Cell-cell interaction based on ligand-receptor database using CellChat v2
 8. Trajectory inference analysis using Monocle 3 
