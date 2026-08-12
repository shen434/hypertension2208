library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)

set.seed(123)

target_gene <- "ATP6V1A"

countMat <- GetAssayData(pbmc, layer = "counts")
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 10000)
hvgs <- VariableFeatures(pbmc)
data <- as.data.frame(countMat[unique(c(target_gene, hvgs)), ])

cell_mean <- colMeans(data, na.rm = TRUE)
top_cell_idx <- order(cell_mean, decreasing = TRUE)[1:min(10000, length(cell_mean))]
data <- data[, top_cell_idx]

result <- scTenifoldKnk(countMatrix = data,
                        gKO = target_gene,
                        qc = TRUE,
                        qc_mtThreshold = 0.1,
                        qc_minLSize = 1000,
                        nc_nNet = 5,
                        nc_nCells = 300)

df <- result$diffRegulation
df <- df[df$gene != target_gene, ]
outTab <- df[df$p.adj < 0.05, ]
write.table(outTab, file = "sigDiff.txt", sep = "\t", quote = FALSE, row.names = FALSE)

top_genes <- head(df[order(-df$FC), ], 20)
p1 <- ggplot(top_genes, aes(x = reorder(gene, FC), y = FC)) +
  geom_bar(stat = 'identity', fill = '#5A9BD4') +
  coord_flip() +
  labs(title = "Top 20 Differentially Regulated Genes", x = "Gene", y = "FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
pdf("barplot_top20.pdf", width = 6, height = 5)
print(p1)
dev.off()

df$log_p.adj <- -log10(df$p.adj)
df$significant <- ifelse(df$p.adj < 0.05, "Significant", "Not significant")
label_genes <- subset(df, p.adj < 0.05)
y_upper <- quantile(df$log_p.adj, 0.999, na.rm = TRUE)
p2 <- ggplot(df, aes(x = Z, y = log_p.adj, color = significant)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray50")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_text_repel(data = label_genes, aes(label = gene), size = 3, max.overlaps = 50) +
  labs(x = "Z-score", y = "-log10(p.adj)") +
  theme_classic() +
  coord_cartesian(ylim = c(0, y_upper)) + theme(legend.position = "none")
pdf("volcano.pdf", width = 6, height = 5)
print(p2)
dev.off()

sig_count <- table(df$significant)
sig_df <- as.data.frame(sig_count)
colnames(sig_df) <- c("category", "count")
sig_df$percentage <- paste0(round(sig_df$count / sum(sig_df$count) * 100, 1), "%")
p3 <- ggplot(sig_df, aes(x = "", y = count, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  geom_text(aes(label = percentage), position = position_stack(vjust = 0.5), size = 4) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = c("Significant" = "red", "Not significant" = "lightgray")) +
  labs(title = "Proportion of Significant Genes", fill = "") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text = element_blank(), panel.grid = element_blank(),
        legend.position = "right", plot.title = element_text(hjust = 0.5))
pdf("pie_significant.pdf", width = 6, height = 5)
print(p3)
dev.off()

pvalueFilter <- 0.05
adjPvalFilter <- 1
colorSel <- if (adjPvalFilter > 0.05) "pvalue" else "p.adjust"

rt <- read.table("sigDiff.txt", header = TRUE, sep = "\t", check.names = FALSE)
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
rt <- cbind(rt, entrezIDs)
rt <- rt[rt[, "entrezIDs"] != "NA", ]
gene <- rt$entrezID

kk_go <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1,
                  qvalueCutoff = 1, ont = "all", readable = TRUE)
GO <- as.data.frame(kk_go)
GO <- GO[GO$pvalue < pvalueFilter & GO$p.adjust < adjPvalFilter, ]
write.table(GO, file = "GO_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("GO_barplot.pdf", width = 9, height = 7)
bar <- barplot(kk_go, drop = TRUE, showCategory = 10, label_format = 100,
               split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bar)
dev.off()

pdf("GO_bubble.pdf", width = 9, height = 7)
bub <- dotplot(kk_go, showCategory = 10, orderBy = "GeneRatio", label_format = 100,
               split = "ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY ~ ., scale = 'free')
print(bub)
dev.off()

ontology.col <- c("#6A89C2FF", "#D6616BFF", "#67B88BFF")
data <- GO[order(GO$pvalue), ]
datasig <- data[data$pvalue < 0.05, ]
BP <- head(datasig[datasig$ONTOLOGY == "BP", , drop = FALSE], 6)
CC <- head(datasig[datasig$ONTOLOGY == "CC", , drop = FALSE], 6)
MF <- head(datasig[datasig$ONTOLOGY == "MF", , drop = FALSE], 6)
data <- rbind(BP, CC, MF)
main.col <- ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio <- Gene / BgGene
logpvalue <- -log(data$pvalue, 10)
logpvalue.col <- brewer.pal(n = 6, name = "Reds")
f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10), colors = logpvalue.col)
BgGene.col <- f(logpvalue)
df_go <- data.frame(GO = data$ID, start = 1, end = max(BgGene))
rownames(df_go) <- df_go$GO
bed2 <- data.frame(GO = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 <- data.frame(GO = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 <- data.frame(GO = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5

pdf("GO_circlize.pdf", width = 10, height = 10)
par(omi = c(0.1, 0.1, 0.1, 1.5))
circos.par(track.margin = c(0.01, 0.01))
circos.genomicInitialize(df_go, plotType = "none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8,
              facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
}
f2 <- colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                         col = value[, 2], border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1],
                                         adj = 0, cex = 0.8, ...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1,
                                         col = '#BA55D3', border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1],
                                         cex = 0.9, adj = 0, ...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white",
                    bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim <- get.cell.meta.data("cell.xlim")
                      cell.ylim <- get.cell.meta.data("cell.ylim")
                      for (j in 1:9) {
                        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[, 1],
                                         col = value[, 2], border = NA, ...)
                    })
circos.clear()

middle.legend <- Legend(
  labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
  type = "points", pch = c(15, 15, 17),
  legend_gp = gpar(col = c('pink', '#BA55D3', ontology.col[1])),
  title = "", nrow = 3, size = unit(3, "mm")
)
circle_size <- unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend <- Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function"),
  type = "points", pch = 15,
  legend_gp = gpar(col = ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 3, size = unit(3, "mm"),
  grid_height = unit(5, "mm"), grid_width = unit(5, "mm")
)
logp.legend <- Legend(
  labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '>=10'),
  type = "points", pch = 16,
  legend_gp = gpar(col = logpvalue.col), title = "-log10(Pvalue)",
  title_position = "topcenter",
  grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd <- packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
dev.off()

kk_kegg <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
kk_kegg@result$Description <- gsub(" - Homo sapiens \\(human\\)", "", kk_kegg@result$Description)
KEGG <- as.data.frame(kk_kegg)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) {
  paste(rt$gene[match(strsplit(x, "/")[[1]], as.character(rt$entrezID))], collapse = "/")
}))
KEGG <- KEGG[KEGG$pvalue < pvalueFilter & KEGG$p.adjust < adjPvalFilter, ]
write.table(KEGG, file = "KEGG_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

showNum <- min(30, nrow(KEGG))
pdf("KEGG_barplot.pdf", width = 8, height = 7)
barplot(kk_kegg, drop = TRUE, showCategory = showNum, label_format = 100, color = colorSel)
dev.off()

pdf("KEGG_bubble.pdf", width = 8, height = 7)
dotplot(kk_kegg, showCategory = showNum, orderBy = "GeneRatio", label_format = 100, color = colorSel)
dev.off()