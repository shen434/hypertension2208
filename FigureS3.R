library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

pvalueFilter <- 0.05
p.adjustFilter <- 0.05

rt <- read.table("module_blue.txt", header = FALSE, sep = "\t", check.names = FALSE)

genes <- unique(as.vector(rt[, 1]))
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezIDs)
gene <- entrezIDs[entrezIDs != "NA"]
gene <- gsub("c\\(\"(\\d+)\".*", "\\1", gene)

kk <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 1, qvalueCutoff = 1, ont = "all", readable = TRUE)
GO <- as.data.frame(kk)
GO <- GO[GO$pvalue < pvalueFilter & GO$p.adjust < p.adjustFilter, ]
write.table(GO, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

options(timeout = 600)
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x) paste(rt$genes[match(strsplit(x, "/")[[1]], as.character(rt$entrezID))], collapse = "/")))
KEGG <- na.omit(KEGG)
KEGG <- KEGG[KEGG$pvalue < pvalueFilter & KEGG$p.adjust < p.adjustFilter, ]
write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

KEGG <- KEGG[, -c(1, 2)]
KEGG <- cbind(ONTOLOGY = "KEGG", KEGG)
GK <- rbind(GO, KEGG)
GK$ONTOLOGY <- factor(GK$ONTOLOGY, levels = c("BP", "CC", "MF", "KEGG"))

data_selected <- GK %>% arrange(ONTOLOGY, pvalue) %>% group_by(ONTOLOGY) %>% slice_head(n = 10) %>% ungroup()
data_selected$log10p <- -log10(data_selected$pvalue)
data_selected <- data_selected %>% arrange(ONTOLOGY, desc(log10p)) %>% 
  mutate(Description = factor(Description, levels = rev(unique(Description))))

ontology_colors <- c("BP" = "#D73027", "CC" = "#1A9850", "MF" = "#4575B4", "KEGG" = "#D4731E")

group_positions <- data_selected %>% group_by(ONTOLOGY) %>% summarise(
  y_start = min(as.numeric(Description)), y_end = max(as.numeric(Description)),
  .groups = 'drop') %>%
  mutate(y_start = y_start - 0.4, y_end = y_end + 0.4)

p <- ggplot(data_selected, aes(x = log10p, y = Description, fill = ONTOLOGY)) +
  geom_col(width = 0.8) +
  geom_rect(data = group_positions,
            aes(xmin = -max(data_selected$log10p) * 0.15, xmax = -max(data_selected$log10p) * 0.07, ymin = y_start, ymax = y_end, fill = ONTOLOGY),
            inherit.aes = FALSE, alpha = 0.8) +
  geom_text(data = group_positions,
            aes(x = -max(data_selected$log10p) * 0.11, y = (y_start + y_end) / 2, label = ONTOLOGY, fill = ONTOLOGY),
            inherit.aes = FALSE, size = 4, fontface = "bold", color = "black") +
  geom_text(aes(label = Description, x = 0.05), hjust = 0, size = 3.5, 
            color = "black", fontface = "plain") + 
  scale_fill_manual(values = ontology_colors, name = "Ontology") + 
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.05)), limits = c(-max(data_selected$log10p) * 0.2, NA)) + 
  labs(title = "Enrichment Analysis", x = "-log10(pvalue)", y = "") +
  theme_void(base_size = 12) +  
  theme(axis.title.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right", legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 30, unit = "pt"),
        axis.line.x = element_line(color = "black", size = 0.3),
        axis.ticks.x = element_line(color = "black"), axis.text.x = element_text(size = 10))

pdf(file = "barplot.pdf", width = 7.5, height = 10)
print(p)
dev.off()

ontology.col <- c("#E178A8", "#A0E070", "#B5E1FF", "#FFB143")
data <- GK[order(GK$pvalue), ]
datasig <- data[data$pvalue < 0.05, , drop = FALSE]
data <- datasig %>% group_by(ONTOLOGY) %>% slice_head(n = 5)
main.col <- ontology.col[as.numeric(as.factor(data$ONTOLOGY))]

BgGene <- as.numeric(sapply(strsplit(data$BgRatio, "/"), '[', 1))
Gene <- as.numeric(sapply(strsplit(data$GeneRatio, '/'), '[', 1))
ratio <- Gene / BgGene
logpvalue <- -log(data$pvalue, 10)
logpvalue.col <- brewer.pal(n = 6, name = "Reds")
f <- colorRamp2(breaks = c(0, 2, 4, 6, 8, 10), colors = logpvalue.col)
BgGene.col <- f(logpvalue)
df <- data.frame(GK = data$ID, start = 1, end = max(BgGene))
rownames(df) <- df$GK
bed2 <- data.frame(GK = data$ID, start = 1, end = BgGene, BgGene = BgGene, BgGene.col = BgGene.col)
bed3 <- data.frame(GK = data$ID, start = 1, end = Gene, BgGene = Gene)
bed4 <- data.frame(GK = data$ID, start = 1, end = max(BgGene), ratio = ratio, col = main.col)
bed4$ratio <- bed4$ratio / max(bed4$ratio) * 9.5

pdf(file = "circlize.pdf", width = 10, height = 10)
par(omi = c(0.1, 0.1, 0.1, 1.5))
circos.par(track.margin = c(0.01, 0.01))
circos.genomicInitialize(df, plotType = "none")
circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), mean(ylim), sector.index, cex = 0.8, facing = "bending.inside", niceFacing = TRUE)
}, track.height = 0.08, bg.border = NA, bg.col = main.col)

for (si in get.all.sector.index()) {
  circos.axis(h = "top", labels.cex = 0.6, sector.index = si, track.index = 1,
              major.at = seq(0, max(BgGene), by = 100), labels.facing = "clockwise")
}
f2 <- colorRamp2(breaks = c(-1, 0, 1), colors = c("green", "black", "red"))
circos.genomicTrack(bed2, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = value[, 2], 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], adj = 0, cex = 0.8, ...)
                    })
circos.genomicTrack(bed3, ylim = c(0, 1), track.height = 0.1, bg.border = "white",
                    panel.fun = function(region, value, ...) {
                      i <- getI(...)
                      circos.genomicRect(region, value, ytop = 0, ybottom = 1, col = '#BA55D3', 
                                         border = NA, ...)
                      circos.genomicText(region, value, y = 0.4, labels = value[, 1], cex = 0.9, adj = 0, ...)
                    })
circos.genomicTrack(bed4, ylim = c(0, 10), track.height = 0.35, bg.border = "white", bg.col = "grey90",
                    panel.fun = function(region, value, ...) {
                      cell.xlim <- get.cell.meta.data("cell.xlim")
                      cell.ylim <- get.cell.meta.data("cell.ylim")
                      for (j in 1:9) {
                        y <- cell.ylim[1] + (cell.ylim[2] - cell.ylim[1]) / 10 * j
                        circos.lines(cell.xlim, c(y, y), col = "#FFFFFF", lwd = 0.3)
                      }
                      circos.genomicRect(region, value, ytop = 0, ybottom = value[, 1], col = value[, 2], 
                                         border = NA, ...)
                    })
circos.clear()

middle.legend <- Legend(
  labels = c('Number of Genes', 'Number of Select', 'Rich Factor(0-1)'),
  type = "points", pch = c(15, 15, 17), legend_gp = gpar(col = c('pink', '#BA55D3', ontology.col[1])),
  title = "", nrow = 3, size = unit(3, "mm")
)
circle_size <- unit(1, "snpc")
draw(middle.legend, x = circle_size * 0.42)

main.legend <- Legend(
  labels = c("Biological Process", "Cellular Component", "Molecular Function", "KEGG"), type = "points", pch = 15,
  legend_gp = gpar(col = ontology.col), title_position = "topcenter",
  title = "ONTOLOGY", nrow = 4, size = unit(3, "mm"), grid_height = unit(5, "mm"),
  grid_width = unit(5, "mm")
)
logp.legend <- Legend(
  labels = c('(0,2]', '(2,4]', '(4,6]', '(6,8]', '(8,10]', '>=10'),
  type = "points", pch = 16, legend_gp = gpar(col = logpvalue.col), title = "-log10(Pvalue)",
  title_position = "topcenter", grid_height = unit(5, "mm"), grid_width = unit(5, "mm"),
  size = unit(3, "mm")
)
lgd <- packLegend(main.legend, logp.legend)
draw(lgd, x = circle_size * 0.85, y = circle_size * 0.55, just = "left")
dev.off()
