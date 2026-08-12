
library(limma)

input_file <- "geneMatrix.txt"
expression_data <- read.table(input_file, header = TRUE, sep = "\t", check.names = FALSE)

if (any(is.na(expression_data[, 1]))) {
  stop("NA in first column")
}

rownames(expression_data) <- make.names(expression_data[, 1], unique = TRUE)
expression_data <- expression_data[, -1]

pdf("1.rawBox.pdf", width = 30, height = 24)
par(mar = c(12, 6, 4, 2))
boxplot(expression_data,
        col = "blue",
        xaxt = "n",
        outline = FALSE,
        cex.axis = 3,
        cex.lab = 3,
        main = "Raw Expression Data")
axis(1, at = 1:ncol(expression_data), labels = colnames(expression_data),
     las = 2, cex.axis = 1.5)
dev.off()

expression_matrix <- as.matrix(expression_data)

quantiles <- quantile(expression_matrix, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
log_condition <- (quantiles[5] > 100) || ((quantiles[6] - quantiles[1]) > 50 && quantiles[2] > 0)

if (log_condition) {
  expression_matrix[expression_matrix < 0] <- 0
  expression_matrix <- log2(expression_matrix + 1)
}

normalized_data <- normalizeBetweenArrays(expression_matrix)
colnames(normalized_data) <- colnames(expression_matrix)

pdf("2.normalBox.pdf", width = 30, height = 24)
par(mar = c(12, 6, 4, 2))
boxplot(normalized_data,
        col = "red",
        xaxt = "n",
        outline = FALSE,
        cex.axis = 3,
        cex.lab = 3,
        main = "Normalized Expression Data")
axis(1, at = 1:ncol(normalized_data), labels = colnames(normalized_data),
     las = 2, cex.axis = 1.5)
dev.off()

normalized_df <- data.frame(geneNames = rownames(normalized_data), normalized_data, check.rows = FALSE)
write.table(normalized_df,
            file = "normalize.txt",
            quote = FALSE,
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)
