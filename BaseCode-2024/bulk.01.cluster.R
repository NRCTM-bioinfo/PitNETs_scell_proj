


library(ggplot2)
library(pheatmap)
library(DESeq2)
library(BiocGenerics)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(genefilter)
library(PoiClaClu)
library(gplots)
library(scatterplot3d)
library(ggpubr)
library(grid)
library(Rtsne)
library(openxlsx)
library(limma)
library(sva)
library(stringr)
library(ggthemes)
library(Rtsne)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggthemes)
library(enrichplot)
library(GSEABase)
library(scales)
library(GSVA)
library(estimate)
library(MCPcounter)
library(xCell)

library(glmnet)
library(dplyr)

library(enrichR)





color.group <- c(G1 = "#20854E", G2 = "#0027B5", G3 = "#7876B1",
                 G4 = "#E18727", G5 = "#BC3C29", G6 = "#6F99AD")


########## GSEA
gsea.path <- c("gsea/msigdb_v7.5.1_GMTs")
gsea.hallmark <- read.gmt(paste0(gsea.path, "/h.all.v7.5.1.symbols.gmt"))
gsea.kegg <- read.gmt(paste0(gsea.path, "/c2.cp.kegg.v7.5.1.symbols.gmt"))
gsea.go <- read.gmt(paste0(gsea.path, "/c5.go.bp.v7.5.1.symbols.gmt"))
gsea.all <- read.gmt(paste0(gsea.path, "/msigdb.v7.5.1.symbols.gmt"))


gene.id <- read.table("gene.id.v40.txt")
gene.id.coding <- gene.id[which(gene.id$V2 == "protein_coding"), ]





############################
vars <- rowVars(exp.data)
names(vars) <- rownames(exp.data)
vars <- vars[order(vars, decreasing = T)]

for (i in 2:10) {

cutoff <- i * 0.02

gene.list <- names(vars)[1:floor(length(vars) * cutoff )]
exp.sub <- exp.data[match(gene.list, rownames(exp.data)), ]

p <- pheatmap(as.matrix(exp.sub), scale = "row",
              color = colorRampPalette(c("#00599F","#00599F","#00599F","white","#D01910","#D01910","#D01910"))(501),
              cluster_row = T, cluster_col = T, border_color = NA,
              cutree_cols = 6, 
              annotation_col = anno.col.fix, 
              annotation_colors = anno.col.color.fix,
              clustering_method = "ward.D2",
              clustering_distance_rows = "manhattan",
              clustering_distance_cols = "manhattan",
              fontsize_col = 3,
              fontsize_row = 1)
ggsave(paste0("output/1.basic/cluster/Heatmap-", cutoff, ".pdf"), p, width = 15, height = 10)

}




















