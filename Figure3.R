library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(wordcloud)
library(ComplexHeatmap)

pvalThreshold <- 0.05
padjThreshold <- 0.05
colorParameter <- if (padjThreshold > 0.05) "pvalue" else "p.adjust"

geneFile <- "genes.csv"
geneData <- read.csv(geneFile, header = TRUE, check.names = FALSE)
geneSymbols <- unique(as.vector(geneData[, 1]))
entrezMapping <- mget(geneSymbols, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrezIDs <- as.character(entrezMapping)
validGenes <- entrezIDs[entrezIDs != "NA"]

goAnalysis <- enrichGO(gene = validGenes,
                       OrgDb = org.Hs.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       ont = "all",
                       readable = TRUE)
goResult <- as.data.frame(goAnalysis)
filteredGO <- goResult[goResult$pvalue < pvalThreshold & goResult$p.adjust < padjThreshold, ]
write.table(filteredGO, file = "GO_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)

pdf("GO_barplot.pdf", width = 8, height = 7)
barPlot <- barplot(goAnalysis, drop = TRUE, showCategory = 10,
                   label_format = 100, split = "ONTOLOGY", color = colorParameter) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_fill_gradientn(colors = c("#FF6666","#FFB266","#FFFF99","#99FF99","#6666FF","#7F52A0","#B266FF"))
print(barPlot)
dev.off()

pdf("GO_bubble.pdf", width = 8, height = 7)
bubblePlot <- dotplot(goAnalysis, showCategory = 10, orderBy = "GeneRatio",
                      label_format = 100, split = "ONTOLOGY", color = colorParameter) +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  scale_color_gradientn(colors = c("#FFB266","#FFFF99","#99FF99","#6666FF","#7F52A0","#B266FF"))
print(bubblePlot)
dev.off()

topGO <- filteredGO %>% group_by(ONTOLOGY) %>% slice_head(n = 10)
pdf("GO_grouped_barplot.pdf", width = 14, height = 10)
groupBarPlot <- ggbarplot(topGO, x = "Description", y = "Count", fill = "ONTOLOGY",
                          color = "white", xlab = "", palette = "aaas",
                          legend = "right", sort.val = "desc", sort.by.groups = TRUE,
                          position = position_dodge(0.9)) +
  rotate_x_text(75) +
  theme(panel.background = element_blank(),
        axis.text.x = element_text(size = 10, color = "black")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  geom_text(aes(label = Count), position = position_dodge(0.9), vjust = -0.3, size = 3) +
  scale_y_continuous(limits = c(0, max(topGO$Count) + 5), expand = c(0, 0))
print(groupBarPlot)
dev.off()

# KEGG enrichment
gene_file_path <- "genes.csv"
raw_gene_data <- read.csv(gene_file_path, header = TRUE, check.names = FALSE)
unique_genes <- unique(as.vector(raw_gene_data[, 1]))
entrez_list <- mget(unique_genes, org.Hs.egSYMBOL2EG, ifnotfound = NA)
entrez_ids <- as.character(entrez_list)
gene_mapping_df <- data.frame(GeneSymbol = unique_genes, EntrezID = entrez_ids, stringsAsFactors = FALSE)
valid_entrez <- gene_mapping_df$EntrezID[gene_mapping_df$EntrezID != "NA"]

kegg_enrich_result <- enrichKEGG(gene = valid_entrez, organism = "hsa",
                                 pvalueCutoff = 1, qvalueCutoff = 1)
kegg_df <- as.data.frame(kegg_enrich_result)
kegg_df$geneID <- as.character(sapply(kegg_df$geneID, function(x) {
  id_vector <- strsplit(x, "/")[[1]]
  paste(gene_mapping_df$GeneSymbol[match(id_vector, gene_mapping_df$EntrezID)], collapse = "/")
}))
filtered_kegg <- kegg_df[kegg_df$pvalue < 1 & kegg_df$p.adjust < 1, ]
write.csv(filtered_kegg, file = "KEGG.csv", quote = FALSE, row.names = FALSE)

top_n <- 30
if (nrow(filtered_kegg) < top_n) top_n <- nrow(filtered_kegg)
top_kegg <- head(filtered_kegg[order(filtered_kegg$pvalue), ], top_n)

pdf("Lollipop_Modified.pdf", width = 9.5, height = 7.5)
lollipop_plot <- ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, color = p.adjust)) +
  geom_segment(aes(xend = Description, y = 0, yend = Count), size = 1.2) +
  geom_point(size = 6) +
  coord_flip() +
  scale_color_gradient(low = "#2c7bb6", high = "#d7191c") +
  labs(title = "Top 30 KEGG Pathways Enriched",
       x = "", y = "Gene Count", color = "Adj P-value") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(lollipop_plot)
dev.off()

pdf("Horizontal_Barplot_KEGG_Annotated_Modified.pdf", width = 9.5, height = 6)
barplot_plot <- ggplot(top_kegg, aes(x = reorder(Description, Count), y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  geom_text(aes(label = Count), hjust = -0.1, size = 5) +
  scale_fill_gradient(low = "#2c7bb6", high = "#d7191c") +
  labs(title = "Top 30 KEGG Pathways",
       x = "", y = "Gene Count", fill = "Adj P-value") +
  theme_bw(base_size = 16) +
  theme(panel.grid.major = element_line(color = "grey85"),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
print(barplot_plot)
dev.off()

pdf("WordCloud_KEGG_Modified.pdf", width = 10, height = 8)
wordcloud(words = top_kegg$Description, freq = top_kegg$Count,
          scale = c(5, 0.5), colors = brewer.pal(11, "Spectral"),
          random.order = FALSE, rot.per = 0.3, use.r.layout = FALSE)
title("KEGG Pathways Word Cloud", cex.main = 2, font.main = 2)
dev.off()

wrap_text <- function(text, max_chars = 25) {
  if (nchar(text) > max_chars) {
    words <- strsplit(text, " ")[[1]]
    lines <- c()
    current_line <- ""
    for (word in words) {
      if (nchar(current_line) + nchar(word) + 1 <= max_chars) {
        if (current_line == "") current_line <- word
        else current_line <- paste(current_line, word, sep = " ")
      } else {
        if (current_line != "") lines <- c(lines, current_line)
        current_line <- word
      }
    }
    if (current_line != "") lines <- c(lines, current_line)
    return(paste(lines, collapse = "\n"))
  }
  return(text)
}

calculate_pathway_association <- function(kegg_data, n_pathways = 10) {
  selected_pathways <- head(kegg_data[order(kegg_data$pvalue), ], n_pathways)
  assoc_matrix <- matrix(0, nrow = n_pathways, ncol = n_pathways)
  rownames(assoc_matrix) <- selected_pathways$Description
  colnames(assoc_matrix) <- selected_pathways$Description
  for (i in 1:n_pathways) {
    for (j in 1:n_pathways) {
      if (i != j) {
        genes_i <- strsplit(selected_pathways$geneID[i], "/")[[1]]
        genes_j <- strsplit(selected_pathways$geneID[j], "/")[[1]]
        overlap_count <- length(intersect(genes_i, genes_j))
        union_count <- length(union(genes_i, genes_j))
        jaccard_similarity <- ifelse(union_count > 0, overlap_count / union_count, 0)
        significance_weight <- (1 - selected_pathways$p.adjust[i]) * (1 - selected_pathways$p.adjust[j])
        association_score <- (jaccard_similarity * significance_weight * 10)
        assoc_matrix[i, j] <- round(association_score, 2)
      } else {
        assoc_matrix[i, j] <- 10
      }
    }
  }
  return(assoc_matrix)
}

if (nrow(top_kegg) >= 6) {
  real_assoc_matrix <- calculate_pathway_association(top_kegg, 6)
} else {
  real_assoc_matrix <- calculate_pathway_association(top_kegg, nrow(top_kegg))
}
wrapped_names <- sapply(rownames(real_assoc_matrix), wrap_text)
rownames(real_assoc_matrix) <- wrapped_names
colnames(real_assoc_matrix) <- wrapped_names
write.csv(real_assoc_matrix, file = "KEGG_Pathway_Association_Matrix.csv", row.names = TRUE)

col_fun <- colorRamp2(c(0, 5, 10), c("#2166AC", "#F7F7F7", "#B2182B"))
pdf("Association_Heatmap_Modified.pdf", width = 10, height = 8)
heatmap_obj <- Heatmap(
  real_assoc_matrix,
  name = "Association Score",
  col = col_fun,
  column_title = paste("KEGG Pathway Association Analysis (Top", matrix_size, "Pathways)"),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  row_title = "KEGG Pathways",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", real_assoc_matrix[i, j]), x, y,
              gp = gpar(fontsize = 9, fontface = "bold",
                       col = ifelse(real_assoc_matrix[i, j] > 5, "white", "#333333")))
  },
  heatmap_legend_param = list(
    title = "Association Score",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10),
    legend_height = unit(4, "cm"),
    legend_width = unit(1.2, "cm"),
    at = c(0, 2.5, 5, 7.5, 10),
    labels = c("0", "2.5", "5.0", "7.5", "10"),
    border = "#333333"
  ),
  border = TRUE,
  border_gp = gpar(col = "#333333", lwd = 0.8),
  rect_gp = gpar(col = "#CCCCCC", lwd = 0.3)
)
draw(heatmap_obj)
dev.off()