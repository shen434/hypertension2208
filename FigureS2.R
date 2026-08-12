
library(Seurat)
library(data.table)
library(stringr)
library(tibble)
library(ggplot2)
library(patchwork)
library(scDblFinder)
library(dplyr)
library(ROGUE)
library(harmony)
library(clustree)
library(SingleR)
library(ggpubr)
library(GSVA)
library(GSEABase)

samples <- list.files("seurat/")
seurat_list <- list()
for (sample in samples) {
  data.path <- paste0("seurat/", sample)
  seurat_data <- Read10X(data.dir = data.path)
  if (is.list(seurat_data)) {
    if ("Gene Expression" %in% names(seurat_data)) {
      counts <- seurat_data[["Gene Expression"]]
    } else {
      counts <- seurat_data[[1]]
    }
  } else {
    counts <- seurat_data
  }
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   project = sample,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list[[sample]] <- seurat_obj
}
gene_sets <- lapply(seurat_list, rownames)
common_genes <- Reduce(intersect, gene_sets)
if (length(common_genes) == 0) stop()
if (length(common_genes) < length(gene_sets[[1]])) {
  seurat_list <- lapply(seurat_list, function(obj) obj[common_genes, ])
}
seurat_combined <- merge(seurat_list[[1]],
                         y = seurat_list[-1],
                         add.cell.ids = names(seurat_list))
pbmc <- JoinLayers(seurat_combined)
print(pbmc)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]] <- PercentageFeatureSet(pbmc, features = HB.genes)
FeatureScatter(pbmc, "nCount_RNA", "percent.mt", group.by = "orig.ident")
FeatureScatter(pbmc, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident")
theme.set2 <- theme(axis.title.x = element_blank())
plot.featrures <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group <- "orig.ident"
plots <- list()
for(i in c(1:length(plot.featrures))){
  plots[[i]] <- VlnPlot(pbmc, group.by = group, pt.size = 0,
                        features = plot.featrures[i]) + theme.set2 + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow = 2)
violin
ggsave("01_vlnplot_before_qc.pdf", plot = violin, width = 14, height = 8)
dim(pbmc)
quantile(pbmc$nFeature_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nFeature_RNA, seq(0.9, 1, 0.01))
quantile(pbmc$nCount_RNA, seq(0.01, 0.1, 0.01))
quantile(pbmc$nCount_RNA, seq(0.9, 1, 0.01))
quantile(pbmc$percent.mt, seq(0.9, 1, 0.01))
quantile(pbmc$percent.HB, seq(0.9, 1, 0.01))
minGene <- 200
maxGene <- 3000
minUMI <- 600
maxUMI <- 20000
pctMT <- 10
pctHB <- 0.1
pbmc <- subset(pbmc, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene &
                 nCount_RNA > minUMI & percent.mt < pctMT & percent.HB < pctHB)
plots <- list()
for(i in seq_along(plot.featrures)){
  plots[[i]] <- VlnPlot(pbmc, group.by = group, pt.size = 0,
                        features = plot.featrures[i]) + theme.set2 + NoLegend()
}
violin <- wrap_plots(plots = plots, nrow = 2)
violin
dim(pbmc)

pbmc <- pbmc %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30) %>%
  RunTSNE(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = 0.8)
sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce,
                   dbr = 0.05,
                   dbr.sd = 0.01,
                   clusters = "ident",
                   nfeatures = 2000,
                   k = 50,
                   verbose = TRUE)
pbmc$Double_score <- sce$scDblFinder.score
pbmc$Is_Double <- ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet")
if (mean(pbmc$Is_Double == "Doublet") > 0.2) {
  pbmc$Is_Double <- ifelse(pbmc$Double_score > 0.6, "Doublet", "Singlet")
}
p1 <- DimPlot(pbmc, reduction = "umap", group.by = "Is_Double")
p2 <- VlnPlot(pbmc, group.by = "Is_Double",
              features = c("nCount_RNA", "nFeature_RNA"),
              pt.size = 0, ncol = 2)
p3 <- ggplot(pbmc@meta.data, aes(x = Double_score, fill = Is_Double)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  theme_classic() +
  scale_fill_manual(values = c("Doublet" = "red", "Singlet" = "blue"))
print(p1); print(p2); print(p3)

pbmc <- NormalizeData(pbmc)
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search = g2m_genes, match = rownames(pbmc))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(pbmc))
pbmc <- CellCycleScoring(pbmc, g2m.features = g2m_genes, s.features = s_genes)
DimPlot(pbmc, group.by = "Phase", reduction = "umap")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"))

DefaultAssay(pbmc) <- "RNA"
top10 <- head(VariableFeatures(pbmc), 10)
pdf("07_top10_variable_features.pdf", width = 7, height = 6)
VariableFeaturePlot(object = pbmc)
dev.off()
pdf("08_top10_labeled.pdf", width = 7, height = 6)
LabelPoints(plot = VariableFeaturePlot(object = pbmc), points = top10, repel = TRUE)
dev.off()
pbmc <- RunPCA(pbmc, verbose = F)
pdf("09_pca_plot.pdf", width = 7, height = 6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()
pdf("10_pca_loadings.pdf", width = 10, height = 9)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()
pdf("11_pca_heatmap.pdf", width = 10, height = 9)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()
pdf("12_elbow_plot.pdf", width = 7, height = 6)
ElbowPlot(pbmc, ndims = 50)
dev.off()
pcs <- 1:40
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident", max_iter = 20)
seq <- seq(0.1, 2, by = 0.1)
pbmc <- FindNeighbors(pbmc, dims = pcs)
for (res in seq) {
  pbmc <- FindClusters(pbmc, resolution = res)
}
p1 <- clustree(pbmc, prefix = "RNA_snn_res.") + coord_flip()
p <- p1 + plot_layout(widths = c(3, 1))
ggsave("13_clustree.png", p, width = 30, height = 14)
res_cols <- grep("RNA_snn_res.", colnames(pbmc@meta.data), value = TRUE)
cluster_counts <- sapply(res_cols, function(x) length(unique(pbmc@meta.data[[x]])))
plot_data <- data.frame(
  Resolution = as.numeric(gsub("RNA_snn_res.", "", res_cols)),
  Clusters = cluster_counts
)
pdf("14_resolution_vs_clusters.pdf", width = 10, height = 6)
ggplot(plot_data, aes(x = Resolution, y = Clusters)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  theme_classic(base_size = 14) +
  labs(title = "Resolution vs Number of Clusters",
       x = "Resolution", y = "Number of Clusters") +
  scale_x_continuous(breaks = seq(0, 2, by = 0.1))
dev.off()
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = pcs) %>% FindClusters(resolution = 0.5)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = pcs) %>% RunTSNE(dims = pcs, reduction = "harmony")
pdf("15_umap_clusters.pdf", width = 7, height = 6)
DimPlot(pbmc, reduction = "umap", label = T)
dev.off()
pdf("16_umap_samples.pdf", width = 7, height = 6)
DimPlot(pbmc, reduction = "umap", label = F, group.by = "orig.ident")
dev.off()
pdf("17_umap_doublets.pdf", width = 7, height = 6)
DimPlot(pbmc, reduction = "umap", label = F, group.by = "Is_Double")
dev.off()
pdf("18_tsne_clusters.pdf", width = 7, height = 6)
DimPlot(pbmc, reduction = "tsne", label = T)
dev.off()
pdf("19_tsne_samples.pdf", width = 7, height = 6)
DimPlot(pbmc, reduction = "tsne", label = F, group.by = "orig.ident")
dev.off()
pdf("20_tsne_doublets.pdf", width = 7, height = 6)
DimPlot(pbmc, reduction = "tsne", label = F, group.by = "Is_Double")
dev.off()

DefaultAssay(pbmc) <- "RNA"
pbmc.markers2 <- FindAllMarkers(pbmc, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(pbmc.markers2, file = "markers.2.RNA.csv")
markers <- c(
  "CD3D", "CD3E", "CD4",
  "CD8A", "CD8B",
  "CD79A", "MS4A1", "CD19",
  "CD14", "FCGR3A", "LYZ", "FCN1",
  "KLRB1", "NKG7", "KLRD1", "GNLY",
  "MZB1", "JCHAIN", "XBP1"
)
p <- DotPlot(pbmc, features = markers) + RotatedAxis()
ggsave("21_dotplot_markers.pdf", p, width = 8, height = 5)

pbmc$celltype.1 <- recode(pbmc@meta.data$seurat_clusters,
                          "0" = "CD4+ T cell",
                          "1" = "CD4+ T cell",
                          "2" = "CD4+ T cell",
                          "3" = "CD8+ T cell",
                          "4" = "CD8+ T cell",
                          "5" = "CD4+ T cell",
                          "6" = "NK cell",
                          "7" = "Monocyte",
                          "8" = "B cell",
                          "9" = "Monocyte",
                          "10" = "CD8+ T cell",
                          "11" = "Plasma cell",
                          "12" = "Unknown")
table(pbmc@meta.data$celltype.1)
Biocols <- c('#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
             '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
             '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
             '#968175')
p <- DimPlot(pbmc, reduction = "umap", label = T, group.by = "celltype.1", cols = Biocols)
ggsave("22_umap_annotated.pdf", p, width = 8, height = 6)
pdf("23_cellType_umap.pdf", width = 8, height = 6)
clusterCornerAxes(object = pbmc,
                  reduction = 'umap',
                  clusterCol = "celltype.1",
                  noSplit = T, cellLabel = T,
                  addCircle = F, cicDelta = 0.1, cicAlpha = 0.2, cicLineSize = 1.5, cicLineColor = "grey50",
                  themebg = 'bwCorner', nbin = 200,
                  cellLabelSize = 5,
                  show.legend = F)
dev.off()
pdf("24_celltype_barplot.pdf", width = 6, height = 6)
cellRatioPlot(
  object = pbmc,
  sample.name = "orig.ident",
  celltype.name = "celltype.1",
  col.width = 0.7
)
dev.off()
celltype_proportion <- pbmc@meta.data %>%
  filter(!is.na(celltype.1)) %>%
  group_by(orig.ident, celltype.1) %>%
  summarise(count = n(), .groups = "drop_last") %>%
  mutate(proportion = count / sum(count) * 100) %>%
  dplyr::select(-count) %>%
  spread(key = celltype.1, value = proportion, fill = 0)
write.csv(celltype_proportion, "celltype_proportion_by_group.csv", row.names = FALSE)

target_genes <- c("MAP2K2", "AKT2", "MAPK3", "PLEC", "MAP2K1")
expr <- GetAssayData(pbmc, layer = "data")
plot_data <- data.frame(
  celltype = pbmc$celltype.1,
  sample = pbmc$orig.ident
)
for (g in target_genes) {
  if (g %in% rownames(expr)) {
    plot_data[[g]] <- as.numeric(expr[g, ])
  }
}
plot_long <- plot_data %>%
  pivot_longer(cols = all_of(target_genes), names_to = "gene", values_to = "expression") %>%
  filter(!is.na(expression))
dot_data <- plot_long %>%
  group_by(celltype, sample, gene) %>%
  summarise(
    avg_expr = mean(expression),
    pct_expr = sum(expression > 0) / n() * 100,
    .groups = "drop"
  )
dot_data$celltype <- factor(dot_data$celltype)
dot_data$gene <- factor(dot_data$gene, levels = target_genes)
pdf("25_schwann_markers.pdf", width = 10, height = 7)
ggplot(dot_data, aes(x = gene, y = celltype)) +
  geom_point(aes(size = pct_expr, color = avg_expr)) +
  facet_grid(~ sample, scales = "free_x", space = "free_x") +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey",
    high = "#B2182B",
    midpoint = 0,
    name = "Avg Exp"
  ) +
  scale_size_continuous(
    range = c(2, 8),
    breaks = c(25, 50, 75, 100),
    labels = c("25%", "50%", "75%", "100%"),
    name = "Pct Exp"
  ) +
  labs(x = NULL, y = "Cell Type", title = "Expression of Schwann cell markers") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(0, 0, 15, 0)),
    strip.background = element_rect(fill = "#E8E8E8", color = "black", size = 1),
    strip.text = element_text(size = 12, face = "bold", color = "black", margin = margin(8, 8, 8, 8)),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    panel.grid.major = element_line(color = "grey92", size = 0.3),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(15, 15, 15, 15)
  ) +
  guides(
    color = guide_colorbar(
      barwidth = unit(1.5, "lines"),
      barheight = unit(8, "lines"),
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "black",
      ticks.colour = "black"
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(color = "grey40")
    )
  )
dev.off()

markers <- c(
  "CD3D", "CD3E", "CD4",
  "CD8A", "CD8B",
  "CD79A", "MS4A1", "CD19",
  "CD14", "FCGR3A", "LYZ", "FCN1",
  "KLRB1", "NKG7", "KLRD1", "GNLY",
  "MZB1", "JCHAIN", "XBP1"
)
expr <- GetAssayData(pbmc, layer = "data")
plot_data <- data.frame(celltype = pbmc$celltype.1)
for (g in markers) {
  if (g %in% rownames(expr)) plot_data[[g]] <- as.numeric(expr[g, ])
}
plot_long <- plot_data %>%
  pivot_longer(cols = all_of(markers), names_to = "gene", values_to = "expression") %>%
  filter(!is.na(expression))
dot_data <- plot_long %>%
  group_by(celltype, gene) %>%
  summarise(
    med_expr = median(expression),
    pct_expr = sum(expression > 0) / n() * 100,
    .groups = "drop"
  )
dot_data$gene <- factor(dot_data$gene, levels = rev(markers))
dot_data$celltype <- factor(dot_data$celltype)
p <- ggplot(dot_data, aes(x = celltype, y = gene)) +
  geom_point(aes(size = pct_expr, color = med_expr)) +
  scale_color_gradient(low = "grey50", high = "red3", name = "Exp") +
  scale_size_continuous(range = c(2, 8), breaks = c(25, 50, 75, 100),
                        labels = c("25%", "50%", "75%", "100%"), name = "Pct Exp") +
  labs(x = "Cell Type", y = NULL, title = "Marker gene expression") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey92"),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "right"
  )
ggsave("26_marker_dotplot_sorted.pdf", plot = p, width = 10, height = 7)

expr <- as.matrix(GetAssayData(pbmc, layer = "data"))
expr <- expr[rowSums(expr) > 1, ]
gmt_file <- "my_genesets.gmt"
geneSets <- getGmt(gmt_file, geneIdType = SymbolIdentifier())
gene_list <- geneIds(geneSets)
gsva_result <- GSVA::gsva(
  expr = expr,
  gset.idx.list = gene_list,
  method = "ssgsea",
  kcdf = "Gaussian",
  min.sz = 1,
  max.sz = Inf,
  parallel.sz = 1,
  verbose = TRUE
)
gsva_score <- t(gsva_result)
for (col_name in colnames(gsva_score)) {
  pbmc <- AddMetaData(pbmc, metadata = gsva_score[, col_name], col.name = col_name)
}
p <- FeaturePlot(pbmc, features = colnames(gsva_score), reduction = "umap", ncol = 1)
p <- p & scale_color_gradient(low = "grey", high = "red")
ggsave("27_umap_ssGSEA.pdf", p, width = 6, height = 5)
pdf("28_vln_ssGSEA.pdf", width = 6, height = 5)
VlnPlot(object = pbmc, features = colnames(gsva_score), ncol = 1)
dev.off()
score_col <- colnames(gsva_score)[1]
plot_df <- data.frame(
  group = as.character(pbmc@meta.data$orig.ident),
  cell_type = as.character(pbmc@meta.data$celltype.1),
  score = as.numeric(pbmc@meta.data[[score_col]]),
  stringsAsFactors = FALSE
)
plot_df <- plot_df[!is.na(plot_df$cell_type) & !is.na(plot_df$score) & !is.na(plot_df$group), ]
plot_df$group <- factor(plot_df$group, levels = c("GSM6564434", "GSM6564435"))
cell_types <- unique(plot_df$cell_type)
stats_list <- lapply(cell_types, function(ct) {
  sub <- plot_df[plot_df$cell_type == ct, ]
  hc <- sub$score[sub$group == "GSM6564434"]
  les <- sub$score[sub$group == "GSM6564435"]
  if(length(hc) >= 3 && length(les) >= 3) {
    test <- wilcox.test(score ~ group, data = sub)
    data.frame(
      cell_type = ct,
      HC_n = length(hc),
      HC_median = round(median(hc), 4),
      lesional_n = length(les),
      lesional_median = round(median(les), 4),
      p_value = test$p.value,
      significance = ifelse(test$p.value < 0.001, "***",
                     ifelse(test$p.value < 0.01, "**",
                     ifelse(test$p.value < 0.05, "*", "ns")))
    )
  } else {
    data.frame(cell_type = ct, HC_n = length(hc), HC_median = NA,
               lesional_n = length(les), lesional_median = NA,
               p_value = NA, significance = "NA")
  }
})
stats_celltype <- do.call(rbind, stats_list)
write.csv(stats_celltype, "Table1_CellType_Wilcoxon_Statistics_ssGSEA.csv", row.names = FALSE)
ct_means <- aggregate(score ~ cell_type, plot_df, mean)
ct_order <- ct_means$cell_type[order(ct_means$score)]
plot_df$cell_type <- factor(plot_df$cell_type, levels = ct_order)
color_palette <- c("GSM6564434" = "#3B6790", "GSM6564435" = "#D66D6D")
p1 <- ggboxplot(
  plot_df,
  x = "cell_type",
  y = "score",
  fill = "group",
  palette = color_palette,
  width = 0.7,
  outlier.shape = NA,
  bxp.errorbar = TRUE,
  bxp.errorbar.width = 0.4
) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "ns")),
    bracket.size = 0.5,
    size = 4.5,
    hide.ns = FALSE
  ) +
  labs(x = "", y = "Schwann cell development and differentiation") +
  rotate_x_text(45) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(size = 11, face = "bold", color = "#2C3E50"),
    axis.text.y = element_text(size = 11, color = "#2C3E50"),
    axis.title.y = element_text(size = 13, face = "bold", color = "#2C3E50"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.3)
  )
overall_test <- wilcox.test(score ~ group, data = plot_df)
stats_overall <- data.frame(
  comparison = "",
  HC_n = sum(plot_df$group == "GSM6564434"),
  HC_median = round(median(plot_df$score[plot_df$group == "GSM6564434"]), 4),
  lesional_n = sum(plot_df$group == "GSM6564435"),
  lesional_median = round(median(plot_df$score[plot_df$group == "GSM6564435"]), 4),
  p_value = overall_test$p.value,
  significance = ifelse(overall_test$p.value < 0.001, "***",
                 ifelse(overall_test$p.value < 0.01, "**",
                 ifelse(overall_test$p.value < 0.05, "*", "ns")))
)
write.csv(stats_overall, "Table2_Overall_Wilcoxon_Statistics_ssGSEA.csv", row.names = FALSE)
p2 <- ggboxplot(
  plot_df,
  x = "group",
  y = "score",
  fill = "group",
  palette = color_palette,
  width = 0.5,
  outlier.shape = NA
) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols = c("***", "**", "*", "ns")),
    comparisons = list(c("GSM6564434", "GSM6564435")),
    bracket.size = 0.6,
    size = 5
  ) +
  labs(x = "", y = "Schwann cell development and differentiation") +
  theme(
    axis.text = element_text(size = 12, face = "bold", color = "#2C3E50"),
    axis.title.y = element_text(size = 13, face = "bold", color = "#2C3E50"),
    legend.position = "none",
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    panel.grid.major.y = element_line(color = "#ECECEC", linewidth = 0.3)
  )
pdf("29_boxplot_celltype_ssGSEA.pdf", width = 12, height = 6)
print(p1)
dev.off()
pdf("30_boxplot_overall_ssGSEA.pdf", width = 5, height = 5)
print(p2)
dev.off()
