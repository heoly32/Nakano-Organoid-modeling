require(dplyr)
require(ggplot2)
require(Seurat)
require(ComplexHeatmap)
require(RColorBrewer)

vis_dimplot <- function(seu, markers, main = NULL, socket = c("data"), reduc = "umap", min.expr = 0, max.expr = NULL, n_col = 3, size = 0.5, split.by = NULL, font_size = 25, apply_log = FALSE) {
  act_assay <- DefaultAssay(seu)
  cnt <- seu[[act_assay]]@data
  if(socket == "scaled") {
    cnt <- seu[[act_assay]]@scale.data
  }
  if(apply_log) {
    cnt <- log(cnt)
  }
  dt <- NULL

  metadata <- seu@meta.data
  if(!is.null(split.by) & !any(colnames(metadata) %in% split.by)) {
    stop("Not found in metadata: ", split.by)
  }

  markers <- markers[markers %in% rownames(cnt)]
  for(marker in markers) {
    temp_dt <- data.frame(row.names = colnames(seu),
                          ident = seu$orig.ident,
                          UMAP_1 = seu[[reduc]]@cell.embeddings[, 1],
                          UMAP_2 = seu[[reduc]]@cell.embeddings[, 2],
                          gene = marker,
                          expr = cnt[marker, ]
                          )
    if(!is.null(split.by)) {
      temp_dt$split_by <- metadata[, split.by]
    }

    if(is.null(dt)) {
      dt <- temp_dt
    } else {
      dt <- rbind(dt, temp_dt)
    }
  }
  dt$gene <- factor(dt$gene, levels = markers)

  myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

  max.expr <- ifelse(is.null(max.expr), max(dt$expr), max.expr)
  sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(min.expr, max.expr))


  g1 <- ggplot(dt, aes(x = UMAP_1, y = UMAP_2, col = expr)) +
    geom_point(size = size) + xlab(paste0(reduc, "_1")) + ylab(paste0(reduc, "_2"))

  if(apply_log) {
    g1 <- g1 + labs(col = "log_expr")
  }
  if(is.null(split.by)) {
    g1 + facet_wrap(. ~ gene, strip.position = "top", ncol = n_col, scales = "free") +
    theme_bw() +
    sc +
    ggtitle(main) +
    theme(strip.text.y.left = element_text(angle = 90), plot.title = element_text(hjust = 0.5, size = font_size), strip.text.x = element_text(size = font_size), strip.text.y = element_text(size = font_size))
  } else {
    g1 + facet_grid(gene ~ split_by, switch = "y") +
    theme_bw() +
    sc +
    ggtitle(main) +
    theme(strip.text.y.left = element_text(angle = 90), plot.title = element_text(hjust = 0.5, size = font_size), strip.text.x = element_text(size = font_size), strip.text.y = element_text(size = font_size))
  }
}

vis_mean_heatmap <- function(seu, genes, assay = "SCT", slot = "data", main = "Heatmap by mean values", mean_by = "seurat_clusters", genes.by = NULL, use.log2 = FALSE) {
  meta <- seu@meta.data %>% as.data.frame()
  meta <- meta[order(meta[, mean_by]), ]
  genes_filt <- genes %in% rownames(seu)
  genes <- genes[genes_filt]
  if(!is.null(genes.by)) { genes.by <- genes.by[genes_filt] }
  ex <- seu[[assay]]@data %>% as.matrix()
  ex <- ex[genes, rownames(meta)]

  if(use.log2) {
    ex <- log2(ex + 1)
  }

  mean_dt <- do.call(cbind,
                     lapply(genes,
                            function(x) {
                              aggr_res <- aggregate(ex[x, ], list(meta[, mean_by]), FUN=mean)
                              t_dt <- data.frame(row.names = aggr_res[, 1], values = aggr_res[, 2])
                              colnames(t_dt) <- x
                              t_dt
                            })
                     )
  #mean_dt <- t(mean_dt)

  Heatmap(mean_dt, row_names_side = "left",
          cluster_rows = F, cluster_columns = F,
          name = "Mean",
          row_split = genes.by,
          row_title_rot = 0,
          col = rev(RColorBrewer::brewer.pal(11, "Spectral")))
}
