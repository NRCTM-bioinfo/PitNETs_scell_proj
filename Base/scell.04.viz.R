




##########################################################################################
############### Percentage
##########################################################################################
out.path.lineage <- paste0(out.path, "/Percentage.Lineage")
cmd <- sprintf("mkdir %s", out.path.lineage)
system(cmd)

ident.meta <- data.frame(hypo@meta.data)
ident.info.meta <- table(ident.meta[, c("Lineage", "Sample")])
ident.info.per <- ident.info.meta

plot.data <- NULL

for (i in 1:ncol(ident.info.meta)) {
  plot.data.sub <- as.data.frame(ident.info.meta[, i])
  colnames(plot.data.sub)[1] = "CellCount"
  plot.data.sub$Percentage <- plot.data.sub$CellCount / sum(plot.data.sub$CellCount) * 100
  plot.data.sub$CellType <- rownames(plot.data.sub)
  plot.data.sub$Sample <- colnames(ident.info.meta)[i]
  plot.data <- rbind(plot.data, plot.data.sub)
  ident.info.per[, i] <-  ident.info.per[, i]/sum(ident.info.per[, i])*100
}

write.xlsx(ident.info.meta, paste0(out.path.lineage, "/total.count.xlsx"), rowNames = T, overwrite = T)
write.xlsx(ident.info.per, paste0(out.path.lineage, "/total.percentage.xlsx"), rowNames = T, overwrite = T)


plot.data$CellType <- factor(plot.data$CellType, levels = names(color.lineage) )

order = colnames(ident.info.meta)
p <- ggbarplot(plot.data, x = "Sample", y = "Percentage", 
          fill = "CellType", color = "CellType", width = 0.7,
          palette = color.lineage,
          size = 0, 
          legend = "right", xlab = "", ylab = "Cell Percentage (%)",
          order = rev(order))
p <- p + theme_base() + coord_flip()
p <- p + theme(axis.text.x=element_text( angle=90 ))
ggsave(paste0(out.path.lineage, "/Percentage.pdf"), p, width = 7, height = 10)

p <- ggbarplot(plot.data, x = "Sample", y = "CellCount", 
          fill = "CellType", color = "CellType", width = 0.7,
          palette = color.lineage,
          size = 0, 
          legend = "right", xlab = "", ylab = "Cell CellCount",
          order = rev(order))
p <- p + theme_base() + coord_flip()
p <- p + theme(axis.text.x=element_text( angle=90 ))
ggsave(paste0(out.path.lineage, "/Count.pdf"), p, width = 7, height = 10)


plot.data$Subtype <- hypo$Subtype[match(plot.data$Sample, hypo$Sample)]


for (i in c("Subtype")) {
plot.data$Type <- plot.data[, i]
plot.data.sub <- plot.data[which(plot.data$Type != ""), ]
#if (i == "Risk") plot.data.sub$Type <- factor(as.character(plot.data.sub$Type), levels = c("Low","High"))
p <- ggboxplot(plot.data.sub, x = "CellType", y = "Percentage",
               fill = "white", color = "Type", size = 1, width = 0.5,
               palette = c("blue","red","green","orange"), main = paste0(i, " anova"),
               add = c("jitter"), add.param = list(color = "Type", size = 1, width = 0.3),
               ylab = "Percentage", xlab = "") + theme_base()
p <- p + stat_compare_means(aes(group = Type, label = paste0(..p.format..)), method = "anova" ) 
p <- p + theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))
ggsave(paste0(out.path.lineage, "/boxplot.", i,".pdf"), p, width = 12, height = 6)
}

write.xlsx(plot.data, paste0(out.path.lineage, "/total.boxplot.xlsx"), overwrite = T)








