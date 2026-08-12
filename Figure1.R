library(limma)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(scatterplot3d)

threshold_logFC <- 1
threshold_adjP <- 0.05
max_display_genes <- 50

file_expr <- "Sample Type Matrix.csv"
tmp_head <- readLines(file_expr, 1)
sep <- ifelse(grepl(",", tmp_head), ",", "\t")
expr_raw <- read.table(file_expr, header=TRUE, sep=sep, check.names=FALSE, stringsAsFactors=FALSE)
rownames(expr_raw) <- expr_raw[,1]
expr_mat <- expr_raw[,-1, drop=FALSE]
expr_mat <- as.matrix(expr_mat)
expr_mat <- apply(expr_mat, 2, as.numeric)
rownames(expr_mat) <- rownames(expr_raw)
colnames(expr_mat) <- colnames(expr_raw)[-1]
combined_expr <- expr_mat

sample_names <- colnames(expr_mat)
group_info <- ifelse(grepl("_Control$", sample_names, ignore.case=TRUE), "Control",
                     ifelse(grepl("_Treat$", sample_names, ignore.case=TRUE), "Treatment", "Unknown"))
num_ctrl <- sum(group_info == "Control")
num_treat <- sum(group_info == "Treatment")

group_labels <- factor(group_info, levels=c("Control", "Treatment"))
design_mat <- model.matrix(~0 + group_labels)
colnames(design_mat) <- c("Control", "Treatment")
fit <- lmFit(expr_mat, design_mat)
contrast_mat <- makeContrasts(Treatment - Control, levels=design_mat)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)
all_diff_results <- topTable(fit2, adjust.method="fdr", number=Inf)

write.table(cbind(Gene=rownames(all_diff_results), all_diff_results),
            file="DE_results.csv", sep=",", quote=FALSE, row.names=FALSE)

significant_DEGs <- all_diff_results[with(all_diff_results,
                                          (abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP)), ]
output_DEGs <- cbind(Gene=rownames(significant_DEGs), significant_DEGs)
SE <- ifelse(as.numeric(output_DEGs[, "t"]) != 0,
             abs(as.numeric(output_DEGs[, "logFC"]) / as.numeric(output_DEGs[, "t"])),
             NA)
output_DEGs <- cbind(output_DEGs, SE=SE)
output_DEGs <- as.data.frame(output_DEGs, stringsAsFactors=FALSE)
desired_order <- c("Gene", "logFC", "SE", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
existing_cols <- intersect(desired_order, colnames(output_DEGs))
output_DEGs <- output_DEGs[, c(existing_cols, setdiff(colnames(output_DEGs), existing_cols))]
write.csv(output_DEGs, file="DE_significant_genes.csv", row.names=FALSE)

significance_status <- ifelse(
  (all_diff_results$adj.P.Val < threshold_adjP & abs(all_diff_results$logFC) > threshold_logFC),
  ifelse(all_diff_results$logFC > threshold_logFC, "Up regulated", "Down regulated"),
  "Not Significant"
)
all_diff_results$Significance <- significance_status

volcano_plot <- ggplot(all_diff_results, aes(x=logFC, y=-log10(adj.P.Val))) +
  geom_point(aes(color=Significance), size=2.5, alpha=0.8) +
  scale_color_manual(values=c("Down regulated"="#1E90FF", "Not Significant"="#808080", "Up regulated"="#FF4500")) +
  geom_vline(xintercept=c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black", linewidth=0.5) +
  geom_hline(yintercept=-log10(threshold_adjP), linetype="dashed", color="black", linewidth=0.5) +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        panel.grid=element_blank(),
        axis.line=element_line(color="black", linewidth=0.5))
ggsave("DE_volcano.pdf", volcano_plot, width=6, height=6)

de_table <- all_diff_results
sig_degs <- de_table[with(de_table, abs(logFC) > threshold_logFC & adj.P.Val < threshold_adjP), ]
up_top <- rownames(sig_degs)[order(-sig_degs$logFC)][1:25]
down_top <- rownames(sig_degs)[order(sig_degs$logFC)][1:25]
label_genes <- unique(c(up_top, down_top))
de_table$Gene <- rownames(de_table)
de_table$label <- ifelse(de_table$Gene %in% label_genes, de_table$Gene, "")
de_table$Group <- "Not Significant"
de_table$Group[de_table$logFC > threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Up regulated"
de_table$Group[de_table$logFC < -threshold_logFC & de_table$adj.P.Val < threshold_adjP] <- "Down regulated"

group_counts <- table(de_table$Group)
up_count <- ifelse("Up regulated" %in% names(group_counts), group_counts["Up regulated"], 0)
down_count <- ifelse("Down regulated" %in% names(group_counts), group_counts["Down regulated"], 0)
not_count <- ifelse("Not Significant" %in% names(group_counts), group_counts["Not Significant"], 0)
legend_labels <- c(paste0("Up regulated (", up_count, ")"),
                   paste0("Down regulated (", down_count, ")"),
                   paste0("Not Significant (", not_count, ")"))

volcano_plot_labeled <- ggplot(de_table, aes(x=logFC, y=-log10(adj.P.Val), color=Group)) +
  geom_point(alpha=0.7, size=2) +
  geom_text_repel(aes(label=label), size=3, max.overlaps=50,
                  box.padding=0.3, point.padding=0.2, segment.color="grey50", show.legend=FALSE) +
  scale_color_manual(values=c("Up regulated"="#FF4500", "Down regulated"="#1E90FF", "Not Significant"="#808080"),
                     breaks=c("Up regulated", "Down regulated", "Not Significant"),
                     labels=legend_labels) +
  geom_vline(xintercept=c(-threshold_logFC, threshold_logFC), linetype="dashed", color="black") +
  geom_hline(yintercept=-log10(threshold_adjP), linetype="dashed", color="black") +
  labs(title="Volcano Plot with Labels", x="Log2 Fold Change", y="-Log10 Adjusted P-value") +
  theme_minimal(base_size=14) +
  theme(plot.title=element_text(face="bold", hjust=0.5),
        panel.grid=element_blank(),
        axis.line=element_line(color="black", linewidth=0.5))
ggsave("DE_volcano_with_labels.pdf", volcano_plot_labeled, width=8, height=8)

diff_gene_list <- data.frame(gene=rownames(significant_DEGs))
write.table(diff_gene_list, file="DEG_geneList.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

pca_result <- prcomp(t(combined_expr), scale.=TRUE)
pca_df <- data.frame(Sample=colnames(combined_expr),
                     PC1=pca_result$x[,1],
                     PC2=pca_result$x[,2],
                     Group=factor(c(rep("Control", num_ctrl), rep("Treatment", num_treat))))
pca_var <- pca_result$sdev^2
pca_var_perc <- round(100 * pca_var / sum(pca_var), 1)
pca_plot <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Group)) +
  stat_ellipse(level=0.95, linetype="dashed", linewidth=1) +
  geom_point(size=4, alpha=0.8, shape=16) +
  geom_text_repel(aes(label=Sample), size=4, fontface="bold", max.overlaps=15) +
  scale_color_manual(values=c("Control"="#0072B2", "Treatment"="#E69F00")) +
  labs(title="PCA Analysis", x=paste("PC1 (", pca_var_perc[1], "%)", sep=""),
       y=paste("PC2 (", pca_var_perc[2], "%)", sep=""), color="") +
  theme_classic(base_size=16) +
  theme(plot.title=element_text(face="bold", hjust=0.5, size=18),
        axis.title=element_text(face="bold", size=16),
        legend.position="top",
        panel.grid=element_blank(),
        axis.line=element_line(color="black", linewidth=0.5))
ggsave("DE_PCA.pdf", pca_plot, width=7, height=7)

pca_3d_df <- data.frame(Sample=colnames(combined_expr),
                        PC1=pca_result$x[,1],
                        PC2=pca_result$x[,2],
                        PC3=pca_result$x[,3],
                        Group=factor(c(rep("Control", num_ctrl), rep("Treatment", num_treat))))
pca_var_3d <- pca_result$sdev^2
pca_var_perc_3d <- round(100 * pca_var_3d / sum(pca_var_3d), 2)
color_bin <- c("#66C2A5", "#FC8D62")
group_colors <- color_bin[as.numeric(pca_3d_df$Group)]
pdf("DE_PCA_3D.pdf", width=10, height=10)
par(mar=c(4,4,6,2))
s3d <- scatterplot3d(x=pca_3d_df$PC2, y=pca_3d_df$PC1, z=pca_3d_df$PC3,
                     color=group_colors, pch=16, cex.symbols=1.8, scale.y=0.7, angle=45,
                     xlab=paste0("PC2 (", pca_var_perc_3d[2], "%)"),
                     ylab=paste0("PC1 (", pca_var_perc_3d[1], "%)"),
                     zlab=paste0("PC3 (", pca_var_perc_3d[3], "%)"),
                     main="", col.axis="#444444", col.grid="#CCCCCC", grid=TRUE, box=TRUE,
                     xlim=range(pca_3d_df$PC2)*1.1, ylim=range(pca_3d_df$PC1)*1.1, zlim=range(pca_3d_df$PC3)*1.1,
                     cex.axis=0.8, cex.lab=1.0)
par(xpd=TRUE)
legend("top", legend=levels(pca_3d_df$Group), col=color_bin, pch=16, title="",
       inset=c(0,-0.1), horiz=TRUE, cex=1.4, pt.cex=1.8, bg="white", box.col="gray60", bty="n")
dev.off()

if(nrow(significant_DEGs) > 0) {
  ordered_DEGs <- significant_DEGs[order(as.numeric(as.vector(significant_DEGs$logFC))), ]
  ordered_gene_names <- rownames(ordered_DEGs)
  total_DEG_count <- length(ordered_gene_names)
  if(total_DEG_count > (max_display_genes * 2)) {
    selected_gene_set <- ordered_gene_names[c(1:max_display_genes, (total_DEG_count - max_display_genes + 1):total_DEG_count)]
  } else {
    selected_gene_set <- ordered_gene_names
  }
  heatmap_expr <- combined_expr[selected_gene_set, ]
  
  sample_annotation <- data.frame(Group=factor(c(rep("Control", num_ctrl), rep("Treat", num_treat))))
  rownames(sample_annotation) <- colnames(combined_expr)
  annotation_colors <- list(Group=c("Control"="#66C2A5", "Treat"="#FC8D62"))
  color_palette <- colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(255)
  pdf("DE_heatmap.pdf", height=8, width=10)
  pheatmap(heatmap_expr, annotation_col=sample_annotation, annotation_colors=annotation_colors,
           color=color_palette, cluster_cols=FALSE, show_colnames=FALSE, scale="row",
           fontsize=12, fontsize_row=7, fontsize_col=10, border_color=NA,
           main=paste("Differential Expression Heatmap\nControl:", num_ctrl, "| Treat:", num_treat))
  dev.off()
  
  heatmap_expr_scaled <- t(scale(t(heatmap_expr)))
  ha_top <- HeatmapAnnotation(
    Distribution = anno_density(heatmap_expr_scaled, type="heatmap", which="column", height=unit(2,"cm")),
    Group = sample_annotation$Group,
    col = annotation_colors,
    annotation_name_side="left",
    annotation_name_gp=gpar(fontsize=12)
  )
  color_palette_comp <- colorRamp2(c(-2,0,2), c("#313695","white","#A50026"))
  ht <- Heatmap(heatmap_expr_scaled, name="Expression\n(Z-score)", col=color_palette_comp,
                top_annotation=ha_top, cluster_columns=FALSE, show_column_names=FALSE,
                row_names_gp=gpar(fontsize=8),
                column_title=paste("Differential Expression Heatmap with Distribution\nControl:", num_ctrl, "| Treat:", num_treat),
                column_title_gp=gpar(fontsize=14, fontface="bold"),
                heatmap_legend_param=list(title_gp=gpar(fontsize=12, fontface="bold"), labels_gp=gpar(fontsize=10)))
  pdf("DE_heatmap_with_distribution.pdf", width=12, height=10)
  draw(ht)
  dev.off()
}


library(limma)
library(WGCNA)

expFile <- "merge.normalizze.txt"
rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)

Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
conCount <- length(Type[Type == "Control"])
treatCount <- length(Type[Type == "Treat"])
datExpr0 <- t(data)

gsg <- goodSamplesGenes(datExpr0, verbose = 0)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(datExpr0), method = "average")
pdf("1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6, mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 12000, col = "red")
dev.off()

clust <- cutreeStatic(sampleTree, cutHeight = 12000, minSize = 10)
datExpr0 <- datExpr0[clust == 1, ]

traitData <- data.frame(Con = c(rep(1, conCount), rep(0, treatCount)),
                        Treat = c(rep(0, conCount), rep(1, treatCount)))
row.names(traitData) <- colnames(data)
sameSample <- intersect(rownames(datExpr0), rownames(traitData))
datExpr0 <- datExpr0[sameSample, ]
datTraits <- traitData[sameSample, ]

sampleTree2 <- hclust(dist(datExpr0), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
pdf("2_sample_heatmap.pdf", width = 12, height = 12)
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

enableWGCNAThreads()
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 0)
pdf("3_scale_independence.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

softPower <- sft$powerEstimate
adjacency <- adjacency(datExpr0, power = softPower)

pdf("3_softConnectivity.pdf", width = 9, height = 5)
k <- softConnectivity(datE = datExpr0, power = softPower)
par(mfrow = c(1, 2))
hist(k)
scaleFreePlot(k, main = "Check Scale free topology")
dev.off()

TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")
pdf("4_gene_clustering.pdf", width = 12, height = 9)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()

minModuleSize <- 50
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
pdf("5_Dynamic_Tree.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

MEList <- moduleEigengenes(datExpr0, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
pdf("6_Clustering_module.pdf", width = 7, height = 6)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.25
abline(h = MEDissThres, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 0)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
pdf("7_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors <- mergedColors
MEs <- mergedMEs

nGenes <- ncol(datExpr0)
nSamples <- nrow(datExpr0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
pdf("8_Module_trait.pdf", width = 6, height = 5.5)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")
traitNames <- names(datTraits)
geneTraitSignificance <- as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", traitNames, sep = "")
names(GSPvalue) <- paste("p.GS.", traitNames, sep = "")

y <- datTraits[, 1]
GS1 <- as.numeric(cor(y, datExpr0, use = "p"))
GeneSignificance <- abs(GS1)
pdf("9_GeneSignificance.pdf", width = 11, height = 7)
plotModuleSignificance(GeneSignificance, mergedColors)
dev.off()

trait <- "Treat"
traitColumn <- match(trait, traitNames)
for (module in modNames) {
  column <- match(module, modNames)
  moduleGenes <- moduleColors == module
  if (nrow(geneModuleMembership[moduleGenes, ]) > 1) {
    outPdf <- paste("10_", trait, "_", module, ".pdf", sep = "")
    pdf(outPdf, width = 7, height = 7)
    par(mfrow = c(1, 1))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, traitColumn]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", trait),
                       main = "Module membership vs. gene significance",
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
    abline(v = 0.8, h = 0.5, col = "red")
    dev.off()
  }
}

probes <- colnames(datExpr0)
geneInfo0 <- data.frame(probes = probes, moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneTraitSignificance[, Tra], GSPvalue[, Tra])
  names(geneInfo0) <- c(oldNames, names(geneTraitSignificance)[Tra], names(GSPvalue)[Tra])
}
for (mod in 1:ncol(geneModuleMembership)) {
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, mod], MMPvalue[, mod])
  names(geneInfo0) <- c(oldNames, names(geneModuleMembership)[mod], names(MMPvalue)[mod])
}
geneOrder <- order(geneInfo0$moduleColor)
geneInfo <- geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.xls", sep = "\t", row.names = FALSE)

for (mod in 1:nrow(table(moduleColors))) {
  modules <- names(table(moduleColors))[mod]
  inModule <- moduleColors == modules
  modGenes <- colnames(datExpr0)[inModule]
  write.table(modGenes, file = paste0("module_", modules, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

geneSigFilter <- 0.5
moduleSigFilter <- 0.5
datMM <- cbind(geneModuleMembership, geneTraitSignificance)
datMM <- datMM[abs(datMM[, ncol(datMM)]) > geneSigFilter, ]
for (mmi in colnames(datMM)[1:(ncol(datMM) - 2)]) {
  dataMM2 <- datMM[abs(datMM[, mmi]) > moduleSigFilter, ]
  write.table(row.names(dataMM2), file = paste0("hubGenes_", mmi, ".txt"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}