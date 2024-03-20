




############### Cell cycle
cc.genes

hypo <- CellCycleScoring(hypo, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

p <- DimPlot(hypo, reduction = "umap", pt.size = 1, label = FALSE, 
             group.by = "Phase", cols = c("#FDBF6F","#6A3D9A","#E7298A")) + theme_few()
ggsave(paste0(out.path, "/7.CC.umap.phase.pdf"), p, width = 9, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = 1, label = FALSE,
             group.by = "Phase", cols = c("#FDBF6F","#6A3D9A","#E7298A")) + theme_few()
ggsave(paste0(out.path, "/7.CC.tsne.phase.pdf"), p, width = 9, height = 7)

p <- VlnPlot(hypo, features = "S.Score", pt.size = 0, 
             group.by = "seurat_clusters", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/7.CC.score.S.Score.pdf"), p, width = 10, height = 5)

p <- VlnPlot(hypo, features = "G2M.Score", pt.size = 0, 
             group.by = "seurat_clusters", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/7.CC.score.G2M.Score.pdf"), p, width = 10, height = 5)






############## annoatation
ref <- readRDS(file = "SingleR/hs.BlueprintEncodeData.RDS")
pred.BlueprintEncodeData.main <- SingleR(test = hypo@assays$RNA@data, ref = ref, labels = ref$label.main)

ref <- readRDS(file = "SingleR/hs.BlueprintEncodeData.RDS")
pred.BlueprintEncodeData.fine <- SingleR(test = hypo@assays$RNA@data, ref = ref, labels = ref$label.fine)

ref <- readRDS(file = "SingleR/hs.HumanPrimaryCellAtlasData.RDS")
pred.HumanPrimaryCellAtlasData.main <- SingleR(test = hypo@assays$RNA@data, ref = ref, labels = ref$label.main)

ref <- readRDS(file = "SingleR/hs.HumanPrimaryCellAtlasData.RDS")
pred.HumanPrimaryCellAtlasData.fine <- SingleR(test = hypo@assays$RNA@data, ref = ref, labels = ref$label.fine)



hypo@meta.data$CellType.BlueprintEncodeData.main <- pred.BlueprintEncodeData.main$labels
hypo@meta.data$CellType.BlueprintEncodeData.fine <- pred.BlueprintEncodeData.fine$labels
hypo@meta.data$CellType.HumanPrimaryCellAtlasData.main <- pred.HumanPrimaryCellAtlasData.main$labels
hypo@meta.data$CellType.HumanPrimaryCellAtlasData.fine <- pred.HumanPrimaryCellAtlasData.fine$labels



p <- DimPlot(hypo, reduction = "tsne", pt.size = 0.5, label = TRUE, label.size = 4, raster = FALSE,
             group.by = "CellType.BlueprintEncodeData.main", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/8.SingleR.tsne.BlueprintEncodeData.main.pdf"), p, width = 10, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = 0.5, label = TRUE, label.size = 4, raster = FALSE,
             group.by = "CellType.BlueprintEncodeData.fine", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/8.SingleR.tsne.BlueprintEncodeData.fine.pdf"), p, width = 12, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = 0.5, label = TRUE, label.size = 4, raster = FALSE,
             group.by = "CellType.HumanPrimaryCellAtlasData.main", cols = color.lib) + theme_few()
ggsave(paste0(out.path, "/8.SingleR.tsne.HumanPrimaryCellAtlasData.main.pdf"), p, width = 12, height = 7)

p <- DimPlot(hypo, reduction = "tsne", pt.size = 0.5, label = TRUE, label.size = 4, raster = FALSE,
             group.by = "CellType.HumanPrimaryCellAtlasData.fine") + theme_few()
ggsave(paste0(out.path, "/8.SingleR.tsne.HumanPrimaryCellAtlasData.fine.pdf"), p, width = 20, height = 7)













################## UMAP
out.path.plot <- paste0(out.path, "/TSNE_marker_stromal")
cmd <- sprintf("mkdir %s", out.path.plot)
system(cmd)

genelist <- c("EPCAM","PFGRB","PTPRC","CD3E","GZMA","CSF3R","CD68","CD79B","KIT","PECAM1","CD4","CD8A","CD7","NCAM1","CD14","CD3G","CD79A","CA1","FUT4","ITGAM","KIT","ANPEP","CD33","S100A9","NGP","MBP","EGFL8","CMA","PRL","POMC","TRIM65","TRIM62","TBX19","POU1F1","NR5A1","GATA2","GH1","GH2","CD19","MS4A1","IGKC","IGHA1","PRLR","ESR1","DRD2","AVPR1B","GATA3","GATA1","IDH1","PTPRC")

genelist <- c("PTPRC","FGF7","MME","PECAM1","VWF","EPCAM","KRT18","PROM1","ALDH1A1","CST3","LUM","DCN","VIM","PDGFRA","COL1A2")

genelist <- unique(genelist)
total.genelist <- rownames(hypo)
genelist <- genelist[genelist %in% total.genelist]
length(genelist)

for (gene.name in genelist) {
  p <- FeaturePlot(object = hypo, features = gene.name, 
              cols = c("#CCCCCC", "red"), pt.size = 0.5, raster = FALSE,
              reduction = "tsne") + theme_base()
  ggsave(paste0(out.path.plot, "/", gene.name,".pdf"), p, width = 8, height = 7)
}








