


################ CNV
library(copykat)


out.path.cnv <- paste0(out.path, "/CNV")
cmd <- sprintf("mkdir %s", out.path.cnv)
system(cmd)

sub.mat <- hypo@assays$RNA@counts[rowSums(hypo@assays$RNA@counts) > 100, ]
sub.mat <- sub.mat[rowSums(sub.mat != 0) > 200 , ]

exp.rawdata <- as.matrix(sub.mat)
set.seed(1)
exp.rawdata <- exp.rawdata[, sample(1:ncol(exp.rawdata), ncol(exp.rawdata))]

sum(rownames(exp.rawdata) %in% c("PTPRC", "LYZ", "PECAM1"))

sum(rownames(sub.mat) %in% c("PTPRC", "LYZ", "PECAM1"))


for (i in 1:10) {
  message(i)
  sub.rawdata <- exp.rawdata[, floor((i-1)*6691.6+1):floor((i)*6691.6) ]
  copykat.result <- copykat::copykat(rawmat = sub.rawdata, 
      sam.name=paste0(out.path.cnv, "/hypo-", i), n.cores = 1)
  saveRDS(copykat.result, file = paste0(out.path.cnv, "/copykat.result.", i, ".rds"))
}


pred.cnv <- NULL
cna.mat <- NULL
cnv.pos <- ""
for (i in 1:10) {
  copykat.result <- readRDS(file = paste0("2.merge/20220707.ds23.merge/CNV/copykat.result.", i, ".rds"))
  pred.cnv <- rbind(pred.cnv, as.data.frame(copykat.result$prediction))
  cnv.pos <- intersect(cnv.pos, paste0(copykat.result$CNAmat$chrom, ":", copykat.result$CNAmat$chrompos))
}
pred.cnv$copykat <- pred.cnv$copykat.pred
#pred.cnv$copykat <- NA
#pred.cnv$copykat[which(pred.cnv$copykat.pred == "diploid")] = "aneuploid"
#pred.cnv$copykat[which(pred.cnv$copykat.pred == "aneuploid")] = "diploid"
write.xlsx(pred.cnv, paste0(out.path.cnv, "/cnv.results.txt"), overwrite = T)




hypo@meta.data$CNVpred <- pred.cnv$copykat[match(rownames(hypo@meta.data), pred.cnv$cell.names)]
hypo@meta.data$CNVpred[which(is.na(hypo@meta.data$CNVpred))] <- "diploid"



table(hypo@meta.data[, c("CNVpred","CNVpred16")])




gexp <- FetchData(hypo, vars = c("seurat_clusters", "Sample","UMAP_1", "UMAP_2","tSNE_1","tSNE_2","CNVpred","Subgroup"))
p <- ggscatter(gexp, x = "tSNE_1", y = "tSNE_2",
          color = "CNVpred", 
          ellipse = F, size = 0.5, 
          palette = c("#003366","#DDDDDD","#FFFFFF"),
          main = "CNVpred",
          mean.point = F,
          star.plot = F)
p <- p + theme_few()
ggsave(paste0(out.path.cnv, "/tSNE.CNVpred.pdf"), p, width = 9, height = 8)


p <- ggscatter(gexp, x = "tSNE_1", y = "tSNE_2",
          color = "CNVpred", 
          ellipse = F, size = 0.2, 
          palette = c("#003366","#DDDDDD","#FFFFFF"),
          main = "CNVpred",
          mean.point = F,
          star.plot = F)
p <- p + theme_few()
p <- p + facet_wrap( ~ Sample, ncol = 6)
ggsave(paste0(out.path.cnv, "/tSNE.CNVpred.sample.split.pdf"), p, width = 16, height = 11)








