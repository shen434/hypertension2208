library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidytable)
library(ggplot2)

geneFile <- "gene.txt"
ctdFile <- "CTD_chem_gene_ixns.csv"
drugFile <- "refer.Drug.txt"

ctdRT <- read.csv(ctdFile, header = FALSE, sep = ",", check.names = FALSE)
geneRT <- read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)

outTab <- ctdRT[ctdRT[, 4] %in% as.vector(geneRT[, 1]), ]
outTab <- outTab[outTab[, 7] == "Homo sapiens", ]
write.csv(outTab, file = "ingredient.csv", row.names = FALSE)

rt <- read.csv("ingredient.csv", header = TRUE, sep = ",", check.names = FALSE)
ctdRT <- rt[, c(1, 4)]
ingredients <- unique(as.vector(rt[, 1]))

drugRT <- read.delim(drugFile, header = TRUE, sep = "\t",
                     check.names = FALSE, fill = TRUE, quote = "\"")
drugRT <- drugRT[, c(2, 3)]

pvalueFilter <- 0.05
adjPvalFilter <- 0.05

kk <- enricher(ingredients,
               pvalueCutoff = 1, qvalueCutoff = 1,
               minGSSize = 5, maxGSSize = 500,
               TERM2GENE = drugRT)

DRUG <- as.data.frame(kk)
DRUG <- DRUG[DRUG$pvalue < pvalueFilter & DRUG$p.adjust < adjPvalFilter, ]
write.csv(DRUG, file = "drug_enrich.csv", row.names = FALSE)

enrichRT <- read.csv("drug_enrich.csv", header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
outTab2 <- data.frame()
for (drug in rownames(enrichRT)) {
  ingred <- enrichRT[drug, "geneID"]
  ingreds <- unlist(strsplit(ingred, "\\/"))
  ingredRT <- unique(ctdRT[ctdRT[, 1] %in% ingreds, ])
  ingredOut <- cbind(drug, ingredRT)
  outTab2 <- rbind(outTab2, ingredOut)
}
colnames(outTab2) <- c("Drug", "Ingredient", "Target")
write.table(outTab2, file = "drug_ingredient_target.txt", sep = "\t", quote = FALSE, row.names = FALSE)