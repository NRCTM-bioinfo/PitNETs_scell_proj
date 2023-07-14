

library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
library(matrixStats)
library(sva)
library(ggpubr)
library(openxlsx)
library(stringr)
library(scran)
library(ggthemes)
library(ggthemes)
library(destiny)
library(grDevices)
library(reticulate)
library(L1Graph)
library(Biobase)
library(scatterplot3d)
library(monocle)
library(pheatmap)
library(harmony)
library(SingleR)
library(Signac)
library(dplyr)

rm(list = ls())

color.lib <- c("#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
               "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#F4B3BE",
               "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
               "#F4A11D", "#8DC8ED", "#4C6CB0", "#8A1C1B", "#CBCC2B", "#EA644C",
               "#634795", "#005B1D", "#26418A", "#CB8A93", "#F1E404", "#E22826",
               "#50C1FF", "#F4D31D", "#F4A11D", "#82C800", "#8B5900", "#858ED1",
               "#FF72E1", "#CB50B2", "#007D9B", "#26418A", "#8B495F", "#FF394B",
               "#873636", "#053d0b")


color.anno <- c(Tumor = "#FDBF6F", Cycling_tumor = "#f7e547",
                M1 = "#f998bf", M2 = "#ff566f", MDSC = "#efb8ec", pDC = "#ff653a", 
                Neutrophil = "#a72ab7", PMN_MDSC = "#ff7ae2",
                Bcell = "#5969e0", CD4_T = "#78ceed", CD8_T = "#198be8", NKcell = "#49d1a8",
                Mast = "#8A1C1B", Endothelial = "#34a84d", Epithelial = "#32960a", Fibroblast = "#0bbfbc")

color.cell <- c(Tumor = "#FDBF6F", 
                MDC = "#cc6a92", Neutrophil = "#a72ab7", B_cell = "#81d7e8", T_NK_cell = "#2077bf", 
                Mast = "#07753e", Endothelial = "#34a84d", Epithelial = "#76d64d", Fibroblast = "#c2ed5e")

color.lineage <- c(Tumor = "#ffc47c", Immune = "#47b6f7", Stromal = "#89dd58")



hypo <- readRDS("obj/PitNETs.rds")
table(hypo@meta.data[, c("Sample","Subtype")])





#############################
pdf(paste0(out.path, "/2.filter.vlnplot.pdf"), width = 25, height = 8)
VlnPlot(hypo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0,
        group.by = "Sample", ncol = 3)
dev.off()

pdf(paste0(out.path, "/2.filter.geneplot.pdf"), width = 12, height = 10)
plot1 <- FeatureScatter(hypo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hypo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

hypo <- NormalizeData(hypo, normalization.method = "LogNormalize", scale.factor = ncol(hypo) )
hypo <- FindVariableFeatures(hypo, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(hypo), 10)

# plot variable features with and without labels
pdf(paste0(out.path, "/3.VariableFeaturePlot.pdf"), width = 12, height = 7)
plot1 <- VariableFeaturePlot(hypo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
dev.off()

all.genes <- rownames(hypo)
hypo <- ScaleData(hypo, features = VariableFeatures(hypo))
hypo <- RunPCA(hypo, features = VariableFeatures(object = hypo))

p <- DimPlot(hypo, reduction = "pca", group.by = "Sample", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/4.PCA.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "pca", group.by = "Subtype", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/4.PCA.Subtype.pdf"), p, width = 8, height = 7)






##################### Harmony
length(VariableFeatures(hypo))
hypo <- FindVariableFeatures(hypo, selection.method = "vst", nfeatures = 3000)
hypo <- ScaleData(hypo, features = VariableFeatures(hypo))
hypo <- RunPCA(hypo, features = VariableFeatures(object = hypo))


set.seed(123)
hypo <- RunHarmony(hypo, "Sample", max.iter.harmony = 5)
hypo <- FindNeighbors(hypo, reduction = "harmony", dims = 1:30)
hypo <- FindClusters(hypo, resolution = 0.5)
hypo <- RunUMAP(hypo, reduction = "harmony", dims = 1:30)
hypo <- RunTSNE(hypo, reduction = "harmony", dims = 1:30, perplexity = 50)







ncol = 6
pt.size = 0.3

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size, cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.UMAP.cluster.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size, label = TRUE, label.size = 8,
             cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.UMAP.cluster.label.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size,
             group.by = "Sample", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.UMAP.sample.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size, 
             group.by = "Subtype", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.UMAP.Subtype.pdf"), p, width = 8, height = 7)

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size, 
             group.by = "seurat_clusters", split.by = "Sample", ncol = ncol, cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.UMAP.sample.split.c.pdf"), p, width = 17, height = 12)



p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.tSNE.cluster.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, label = TRUE, label.size = 8,
             cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.tSNE.cluster.label.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, 
             group.by = "Sample", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.tSNE.sample.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, 
             group.by = "Subtype", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.tSNE.Subtype.pdf"), p, width = 8, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, 
             group.by = "seurat_clusters", split.by = "Sample", ncol = ncol, cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/5.tSNE.sample.split.c.pdf"), p, width = 17, height = 12)






p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, label = TRUE, label.size = 6,
             group.by = "CellType", cols = color.cell) + theme_base()
ggsave(paste0(out.path, "/6.tSNE.cell.label.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size, label = TRUE, label.size = 6,
             group.by = "CellType", cols = color.cell) + theme_base()
ggsave(paste0(out.path, "/6.UMAP.cell.label.pdf"), p, width = 9, height = 7)





p <- DimPlot(hypo, reduction = "tsne", pt.size = pt.size, label = FALSE, label.size = 6,
             group.by = "Lineage", cols = color.lineage) + theme_base()
ggsave(paste0(out.path, "/6.tSNE.lineage.label.pdf"), p, width = 8.5, height = 7)

p <- DimPlot(hypo, reduction = "umap", pt.size = pt.size, label = FALSE, label.size = 6,
             group.by = "Lineage", cols = color.lineage) + theme_base()
ggsave(paste0(out.path, "/6.UMAP.lineage.label.pdf"), p, width = 8.5, height = 7)


































