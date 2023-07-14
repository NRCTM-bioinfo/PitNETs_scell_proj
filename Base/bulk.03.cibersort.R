######################## cibersort
plot.data <- read.xlsx('cibersort.tumor.xlsx')

plot.data$Subgroups <- meta.data$Subgroups[match(plot.data$Mixture, meta.data$SampleID)]
rownames(plot.data) <- plot.data$Mixture
plot.data[, 2:11] <- plot.data[, 2:11]/plot.data$`Absolute.score.(sig.score)`
p <- pheatmap(as.matrix(t(plot.data[, 2:11])), scale = "row",
              color = colorRampPalette(c("#00599F","#00599F","#00599F","white","#D01910","#D01910","#D01910"))(501),
              cluster_row = F, cluster_col = F, border_color = NA,
              cutree_cols = 6, 
              annotation_col = anno.col.fix, 
              annotation_colors = anno.col.color.fix,
              clustering_method = "ward.D2",
              clustering_distance_rows = "manhattan",
              cluster_cols = as.hclust(col_dend),
              fontsize_col = 10,
              fontsize_row = 10)
ggsave(paste0("cibersort.tam.Heatmap.pdf"), p, width = 14, height = 10)

for ( i in 2:11) {
  plot.data$Exp <- plot.data[, i]
  p1 <- ggboxplot(plot.data, 
                  x = "Subgroups", y = "Exp",
                  color = "black", fill = "white",
                  palette = color.group,
                  order = names(color.group),
                  add = "jitter", add.param = list(color = "Subgroups", size = 2),
                  xlab = "", ylab = paste0("Percentage"),
                  main = paste0(colnames(plot.data)[i], " Subgroups"),
                  legend = "bottom", outlier.shape = NA)
  p1 <- p1 + stat_compare_means(aes(label = ..p.format..), method = "wilcox.test",
                                ref.group = ".all.", hide.ns = F, size = 3) 
  p1 <- p1 + theme_base() + theme(plot.background = element_blank())
  ggsave(paste0("box.", colnames(plot.data)[i],".pdf"), p1, width = 8, height = 5)
  
}



plot.info <- NULL
sub <- data.frame(Patient = rownames(plot.data),
                  CellType = "Tumor_POU1F1",
                  Percent = as.numeric(plot.data[, "Tumor_POU1F1"]))
plot.info <- rbind(plot.info, sub)
sub <- data.frame(Patient = rownames(plot.data),
                  CellType = "Tumor_TBX19",
                  Percent = as.numeric(plot.data[, "Tumor_TBX19"]))
plot.info <- rbind(plot.info, sub)
sub <- data.frame(Patient = rownames(plot.data),
                  CellType = "Tumor_NR5A1",
                  Percent = as.numeric(plot.data[, "Tumor_NR5A1"]))
plot.info <- rbind(plot.info, sub)
sub <- data.frame(Patient = rownames(plot.data),
                  CellType = "Immune",
                  Percent = as.numeric(rowSums(plot.data[, c("Neutrophil","Mast","T_NK_cell","MDC")])))
plot.info <- rbind(plot.info, sub)
sub <- data.frame(Patient = rownames(plot.data),
                  CellType = "Stromal",
                  Percent = as.numeric(rowSums(plot.data[, c("Fibroblast","Epithelial","Endothelial")])))
plot.info <- rbind(plot.info, sub)


plot.info$CellType <- factor(as.character(plot.info$CellType), levels = rev( c("Tumor_POU1F1","Tumor_TBX19","Tumor_NR5A1","Immune","Stromal") ) )
p <- ggbarplot(plot.info, x = "Patient", y = "Percent", 
               fill = "CellType", color = "CellType", width = 0.9,
               palette = c(Tumor_POU1F1 = "#ffbfb2", Tumor_TBX19 = "#ffab44", Tumor_NR5A1 = "#f2dc71", Stromal = "#32960a", Immune = "#c95080"),
               size = 0, 
               order = meta.data$SampleID,
               legend = "right", xlab = "", ylab = "Cell Percentage (%)")
p <- p + theme_bw() 
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 3))
p <- p + scale_y_continuous(expand = c(0,0))
ggsave(paste0("cibersrot.percentage.pdf"), p, width = 20, height = 5)

plot.info$Subtype <- meta.data$Subtype[match(plot.info$Patient, meta.data$SampleID)]
plot.info$Subtype <- factor(as.character(plot.info$Subtype), levels = c("PRL_GH","TPIT","NFPA"))
p2 <- ggboxplot(plot.info, 
               x = "CellType", y = "Percent",
               color = "black", fill = "Subtype",
               palette = color.subtype, width = 0.6,
               bxp.errorbar = TRUE,
               order = c("Tumor_POU1F1","Tumor_TBX19","Tumor_NR5A1","Immune","Stromal"),
               xlab = "", ylab = paste0("Cell percentage"),
               main = paste0(" Subtype"),
               legend = "bottom", outlier.shape = NA)
p2 <- p2 + stat_compare_means(aes(group = Subtype), method = "anova", size = 2) 
p2 <- p2 + theme_base() + theme(plot.background = element_blank())
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5, size = 10))
ggsave(paste0("cibersrot.percentage.box.pdf"), p2 , width = 8, height = 5)
write.xlsx(plot.info, "cibersrot.percentage.xlsx")





