

cell.marker.list <- readRDS("cell.marker.list.rds")

ssgsea.score <- gsva(as.matrix(exp.vsd), cell.marker.list, method='ssgsea', 
                     kcdf='Gaussian', abs.ranking=TRUE)
scell.mat <- as.data.frame(ssgsea.score)


p <- pheatmap(as.matrix(scell.mat), scale = "row",
              color = colorRampPalette(c(rep("blue",2),"white",rep("red",2)))(501),
              cluster_row = F, cluster_col = F, border_color = NA,
              annotation_col = anno.col.fix, 
              annotation_colors = anno.col.color.fix,
              clustering_method = "ward.D",
              clustering_distance_rows = "manhattan",
              fontsize_col = 3,
              fontsize_row = 12)






