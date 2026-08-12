library(glmnet)
library(pROC)
library(limma)
library(reshape2)
library(ggpubr)

expFile <- "GSE234085.normalize.txt"
geneFile <- "gene.txt"
diseaseName <- "HTN"

rt <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
y <- gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y <- ifelse(y == "Control", 0, 1)

geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)

bioCol <- rainbow(nrow(geneRT), s = 0.9, v = 0.9)
aucText <- c()
pdf("ROC_genes.pdf", width = 5, height = 4.75)
for (k in 1:nrow(geneRT)) {
  x <- as.vector(geneRT[k, 1])
  roc1 <- roc(y, as.numeric(rt[x, ]))
  if (k == 1) {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "")
  } else {
    plot(roc1, print.auc = FALSE, col = bioCol[k], legacy.axes = TRUE, main = "", add = TRUE)
  }
  aucText <- c(aucText, paste0(x, ", AUC=", sprintf("%.3f", roc1$auc[1])))
}
legend("bottomright", aucText, lwd = 2, bty = "n", col = bioCol[1:nrow(geneRT)])
dev.off()

rt <- rt[as.vector(geneRT[, 1]), ]
rt <- as.data.frame(t(rt))
logit <- glm(y ~ ., family = binomial(link = 'logit'), data = rt)
pred <- predict(logit, newx = rt)
roc1 <- roc(y, as.numeric(pred))
ci1 <- ci.auc(roc1, method = "bootstrap")
ciVec <- as.numeric(ci1)
pdf("ROC_model.pdf", width = 5, height = 4.75)
plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = "Model")
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")
dev.off()

rt2 <- read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt2 <- as.matrix(rt2)
rownames(rt2) <- rt2[, 1]
exp <- rt2[, 2:ncol(rt2)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data <- avereps(data)
geneRT2 <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data <- t(data[as.vector(geneRT2[, 1]), ])
Type <- gsub("(.*)\\_(.*)", "\\2", row.names(data))
Type <- ifelse(Type == "Control", "Control", diseaseName)
Type <- factor(Type, levels = c("Control", diseaseName))
rt_box <- cbind(as.data.frame(data), Type)
data_melt <- melt(rt_box, id.vars = c("Type"))
colnames(data_melt) <- c("Type", "Gene", "Expression")

p <- ggboxplot(data_melt, x = "Gene", y = "Expression", fill = "Type",
               xlab = "", ylab = "Gene expression", legend.title = "Type",
               palette = c("#3B6790", "#D66D6D"), width = 0.6)
p <- p + rotate_x_text(45)
p1 <- p + stat_compare_means(aes(group = Type),
                             method = "wilcox.test",
                             symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                                symbols = c("***", "**", "*", " ")),
                             label = "p.signif")
pdf("boxplot.pdf", width = 5, height = 4.5)
print(p1)
dev.off()