library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)
library(assertthat)
library(SCpubr)

logFCfilter <- 1
adjPvalFilter <- 0.05

dirs <- list.dirs(".")
dirs_sample <- dirs[-1]
names(dirs_sample) <- gsub(".+\\/(.+)", "\\1", dirs_sample)
counts <- Read10X(data.dir = dirs_sample)
pbmc <- CreateSeuratObject(counts, min.cells = 3, min.features = 100)

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^mt-")

pdf("01_QC_vlnplot.pdf", width = 10, height = 6.5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

pbmc <- subset(pbmc,
               subset = nFeature_RNA > 200 &
                        nFeature_RNA < 10000 &
                        percent.mt < 15 &
                        nCount_RNA > 200)

pdf("01_QC_scatter.pdf", width = 13, height = 7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf("01_variable_features.pdf", width = 10, height = 6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(object = pbmc, npcs = 20, pc.genes = VariableFeatures(object = pbmc))
pbmc <- RunHarmony(pbmc, "orig.ident")

pdf("02_PCA_loadings.pdf", width = 10, height = 8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

pdf("02_PCA_plot.pdf", width = 7.5, height = 5)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

pdf("02_PCA_heatmap.pdf", width = 10, height = 8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf("02_PC_jackstraw.pdf", width = 8, height = 6)
JackStrawPlot(object = pbmc, dims = 1:20)
dev.off()

pcSelect <- 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)
pbmc <- FindClusters(pbmc, resolution = seq(0.5, 1.2, by = 0.1))
pbmc <- FindClusters(object = pbmc, resolution = 0.6)

res_cols <- grep("RNA_snn_res", colnames(pbmc@meta.data), value = TRUE)
pdf("03_clustree.pdf", width = 10, height = 8)
clustree::clustree(pbmc, prefix = "RNA_snn_res.",
                   node_colour = "sc3_stability",
                   node_size_range = c(4, 10),
                   edge_width = 0.8,
                   show_axis = TRUE) +
  ggtitle("Clustering Resolution Hierarchy (0.5 to 1.2)") +
  theme(legend.position = "right")
dev.off()

pdf("03_umap.pdf", width = 7, height = 6)
pbmc <- RunUMAP(object = pbmc, dims = 1:pcSelect)
DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)
dev.off()

write.table(pbmc$seurat_clusters, file = "03_clusters.txt", quote = FALSE, sep = "\t", col.names = FALSE)

pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers <- pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC))) > logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj)) < adjPvalFilter), ]
write.table(sig.markers, file = "03_cluster_markers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf("03_marker_heatmap.pdf", width = 15, height = 15)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

af <- pbmc
Idents(pbmc) <- "seurat_clusters"
Idents(af) <- "seurat_clusters"

refMarker <- read.table("cellmarker.txt", header = FALSE, sep = "\t", check.names = FALSE)
genes <- list()
for(i in 1:nrow(refMarker)){
  genes[[refMarker[i,1]]] <- unlist(strsplit(refMarker[i,2], "\\,"))
}

pdf("03_marker_dotplot.pdf", width = 15, height = 8)
do_DotPlot(sample = af, features = genes, dot.scale = 12,
           legend.length = 50, legend.framewidth = 2, font.size = 12)
dev.off()

annCells <- c("NK Cell", "Monocyte", "NK Cell", "T cell", "T cell",
              "T cell", "T cell", "NK Cell", "B Cell", "Monocyte",
              "Monocyte", "T cell", "Platelet", "B Cell", "Monocyte",
              "Platelet", "T cell", "NK Cell", "Dendritic cell", "NK Cell",
              "Plasma cell", "Monocyte")

clusterAnn <- data.frame(id = levels(Idents(af)), labels = annCells)
clusterAnn <- clusterAnn[, c("id", "labels")]
write.table(clusterAnn, file = "04_cluster_annotation.txt", quote = FALSE, sep = "\t", row.names = FALSE)

cellAnn <- c()
for(i in 1:length(pbmc$seurat_clusters)){
  index <- pbmc$seurat_clusters[i]
  cellAnn <- c(cellAnn, clusterAnn[index, 2])
}
cellAnnOut <- cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut) <- c("id", "labels")
write.table(cellAnnOut, file = "04_cell_annotation.txt", quote = FALSE, sep = "\t", row.names = FALSE)

newLabels <- annCells
names(newLabels) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, newLabels)
pdf("04_umap_annotated.pdf", width = 7.5, height = 6)
DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)
dev.off()

new_cell_type_mapping <- c(
  "0" = "NK Cell", "1" = "Monocyte", "2" = "NK Cell",
  "3" = "T cell", "4" = "T cell", "5" = "T cell",
  "6" = "T cell", "7" = "NK Cell", "8" = "B Cell",
  "9" = "Monocyte", "10" = "Monocyte", "11" = "T cell",
  "12" = "Platelet", "13" = "B Cell", "14" = "Monocyte",
  "15" = "Platelet", "16" = "T cell", "17" = "NK Cell",
  "18" = "Dendritic cell", "19" = "NK Cell", "20" = "Plasma cell",
  "21" = "Monocyte"
)
metadata <- as.data.frame(pbmc@meta.data)
metadata <- metadata %>%
  mutate(cell_type = recode(as.character(seurat_clusters), !!!new_cell_type_mapping))
pbmc@meta.data <- metadata