library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

source("refer.ML.R")

Train_data <- read.table("data.train.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)
Train_expr <- Train_data[,1:(ncol(Train_data)-1),drop=F]
Train_class <- Train_data[,ncol(Train_data),drop=F]

Test_data <- read.table("data.test.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)
Test_expr <- Test_data[,1:(ncol(Test_data)-1),drop=F]
Test_class <- Test_data[,ncol(Test_data),drop=F]
Test_class$Cohort <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", row.names(Test_class))
Test_class <- Test_class[,c("Cohort", "Type")]

comgene <- intersect(colnames(Train_expr), colnames(Test_expr))
Train_expr <- as.matrix(Train_expr[,comgene])
Test_expr <- as.matrix(Test_expr[,comgene])
Train_set <- scaleData(data=Train_expr, centerFlags=T, scaleFlags=T)
Test_set <- scaleData(data=Test_expr, cohort=Test_class$Cohort, centerFlags=T, scaleFlags=T)

methodRT <- read.table("refer.methodLists.txt", header=T, sep="\t", check.names=F)
methods <- gsub("-| ", "", methodRT$Model)

classVar <- "Type"
min.selected.var <- 2
Variable <- colnames(Train_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))

preTrain.var <- list()
set.seed(123)
for (method in preTrain.method) {
  var_result <- RunML(method=method, Train_set=Train_set, Train_label=Train_class,
                      mode="Variable", classVar=classVar)
  preTrain.var[[method]] <- var_result
}
preTrain.var[["simple"]] <- colnames(Train_set)

empty_methods <- names(preTrain.var)[sapply(preTrain.var, length)==0]
if (length(empty_methods)>0) {
  methods_to_remove <- sapply(methods, function(m) {
    parts <- strsplit(m, "\\+")[[1]]
    if (length(parts)==1) parts <- c("simple", parts)
    parts[1] %in% empty_methods
  })
  methods <- methods[!methods_to_remove]
}

model <- list()
set.seed(123)
Train_set_bk <- Train_set
for (method in methods) {
  method_name <- method
  method <- strsplit(method, "\\+")[[1]]
  if (length(method)==1) method <- c("simple", method)
  Variable <- preTrain.var[[method[1]]]
  if (length(Variable)==0) next
  Train_set <- Train_set_bk[, Variable, drop=F]
  fit <- RunML(method=method[2], Train_set=Train_set, Train_label=Train_class,
               mode="Model", classVar=classVar)
  selected_var <- ExtractVar(fit)
  if (length(selected_var) > min.selected.var) {
    model[[method_name]] <- fit
  }
}
Train_set <- Train_set_bk
saveRDS(model, "model.MLmodel.rds")

FinalModel <- "multiLogistic"
if (FinalModel == "multiLogistic") {
  logisticmodel <- lapply(model, function(fit) {
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family="binomial",
               data=as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })
}
saveRDS(logisticmodel, "model.logisticmodel.rds")

model <- readRDS("model.MLmodel.rds")
methodsValid <- names(model)
if (length(methodsValid)==0) stop("No valid models")

RS_list <- list()
for (method in methodsValid) {
  RS_list[[method]] <- CalPredictScore(fit=model[[method]],
                                       new_data=rbind.data.frame(Train_set,Test_set))
}
riskTab <- as.data.frame(t(do.call(rbind, RS_list)))
riskTab <- cbind(id=row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep="\t", row.names=F, quote=F)

Class_list <- list()
for (method in methodsValid) {
  Class_list[[method]] <- PredictClass(fit=model[[method]],
                                       new_data=rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
classTab <- cbind(id=row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep="\t", row.names=F, quote=F)

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}
fea_df <- lapply(model, function(fit) data.frame(ExtractVar(fit)))
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, "model.genes.txt", sep="\t", row.names=F, col.names=T, quote=F)

AUC_list <- list()
for (method in methodsValid) {
  auc_result <- RunEval(fit=model[[method]], Test_set=Test_set, Test_label=Test_class,
                        Train_set=Train_set, Train_label=Train_class,
                        Train_name="Train", cohortVar="Cohort", classVar=classVar)
  AUC_list[[method]] <- auc_result
}
AUC_mat <- do.call(rbind, AUC_list)
aucTab <- cbind(Method=row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep="\t", row.names=F, quote=F)

library(ComplexHeatmap)
library(RColorBrewer)

AUC_mat <- read.table("model.AUCmatrix.txt", header=T, sep="\t", check.names=F, row.names=1, stringsAsFactors=F)
avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing=T)
AUC_mat <- AUC_mat[names(avg_AUC),]
avg_AUC <- as.numeric(format(avg_AUC, digits=3, nsmall=3))
CohortCol <- brewer.pal(n=ncol(AUC_mat), name="Paired")
names(CohortCol) <- colnames(AUC_mat)
cellwidth <- 1
cellheight <- 0.5

hm <- SimpleHeatmap(Cindex_mat=AUC_mat, avg_Cindex=avg_AUC,
                    CohortCol=CohortCol, barCol="steelblue",
                    cellwidth=cellwidth, cellheight=cellheight,
                    cluster_columns=F, cluster_rows=F)
pdf("model.AUCheatmap.pdf", width=cellwidth*ncol(AUC_mat)+6, height=cellheight*nrow(AUC_mat)*0.45)
draw(hm, heatmap_legend_side="right", annotation_legend_side="right")
dev.off()

library(pROC)

rsFile <- "model.riskMatrix.txt"
method <- "RF+Enet[alpha=0.8]"
riskRT <- read.table(rsFile, header=T, sep="\t", check.names=F, row.names=1)
CohortID <- gsub("(.*)\\_(.*)\\_(.*)", "\\1", rownames(riskRT))
CohortID <- gsub("(.*)\\.(.*)", "\\1", CohortID)
riskRT$Cohort <- CohortID

for (Cohort in unique(riskRT$Cohort)) {
  rt <- riskRT[riskRT$Cohort==Cohort,]
  y <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
  y <- ifelse(y=="Control", 0, 1)
  roc1 <- roc(y, as.numeric(rt[,method]))
  ci1 <- ci.auc(roc1, method="bootstrap")
  ciVec <- as.numeric(ci1)
  pdf(paste0("ROC.", Cohort, ".pdf"), width=5, height=4.75)
  plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=Cohort)
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
  dev.off()
}
library(limma)
library(pheatmap)

expFile <- "merge.normalzie.txt"
geneFile <- "hubGenes.csv"
rt <- read.table(expFile, header=T, sep="\t", check.names=F)
rt <- as.matrix(rt)
rownames(rt) <- rt[,1]
exp <- rt[,2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
data <- matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data <- avereps(data)
geneRT <- read.csv(geneFile, header=T, check.names=F)
data <- data[as.vector(geneRT[,1]),]
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
data <- data[,order(Type)]
Project <- gsub("(.+)\\_(.+)\\_(.+)", "\\1", colnames(data))
Type <- gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(data))
colnames(data) <- gsub("(.+)\\_(.+)\\_(.+)", "\\2", colnames(data))
names(Type) <- colnames(data)
Type <- as.data.frame(Type)
Type <- cbind(Project, Type)
pdf("hubGene.heatmap.pdf", width=7.5, height=3.5)
pheatmap(data, annotation=Type,
         color=colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(50),
         cluster_cols=F, show_colnames=F, scale="row",
         fontsize=8, fontsize_row=7, fontsize_col=8)
dev.off()

library(dplyr)
library(data.table)
library(coloc)
library(VariantAnnotation)
library(gwasglue)
library(locuscomparer)

drugFile <- "IEU.ref.txt"
outcomeFile <- "finngen_R12_I9_HYPTENS.gz"
finnAnnFile <- "finngen_R12_manifest.tsv"

drugRT <- read.table(drugFile, header=T, sep="\t", check.names=F)
row.names(drugRT) <- drugRT[,"id"]
allVcf <- list.files(pattern="*.vcf.gz$")

data0 <- data.table::fread(outcomeFile, header=T, sep="\t", check.names=F)
finn_info <- fread(finnAnnFile, data.table=F)
trait_row <- finn_info[grepl(outcomeFile, finn_info$path_https),]
data0$ncase.outcome <- trait_row$num_cases
data0$ncontrol.outcome <- trait_row$num_controls
data0$samplesize.outcome <- trait_row$num_cases + trait_row$num_controls

data1 <- data0 %>% dplyr::select("rsids","#chrom","pos","alt","ref","af_alt","beta","sebeta","pval","samplesize.outcome","ncase.outcome")
colnames(data1) <- c('SNP','chrom',"pos",'effect_allele','other_allele',"eaf","beta","se","P","samplesize","number_cases")
data2 <- as.data.frame(data1)
data2$varbeta <- data2$se^2
data2$MAF <- ifelse(data2$eaf<0.5, data2$eaf, 1-data2$eaf)
data3 <- subset(data2, !duplicated(SNP))
data3$s <- data3$number_cases/data3$samplesize
data3$z <- data3$beta/data3$se
GWASdataAll <- data3 %>% na.omit()

outTab <- data.frame()
for (eqtlFile in allVcf) {
  eqtlID <- gsub(".vcf.gz", "", eqtlFile)
  geneChr <- drugRT[eqtlID, "Chr"]
  geneStart <- drugRT[eqtlID, "Start"]
  geneEnd <- drugRT[eqtlID, "End"]
  geneName <- drugRT[eqtlID, "Symbol"]
  vcfRT <- readVcf(eqtlFile)
  data1 <- gwasvcf_to_TwoSampleMR(vcf=vcfRT, type="exposure")
  data2 <- data1 %>% dplyr::select("SNP","chr.exposure","pos.exposure","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","se.exposure","pval.exposure","samplesize.exposure")
  colnames(data2) <- c("SNP","chrom","Pos","effect_allele","other_allele","MAF","beta","se","P","samplesize")
  data3 <- as.data.frame(data2)
  data3$varbeta <- data3$se^2
  data3$z <- data3$beta/data3$se
  data3 <- subset(data3, !duplicated(SNP))
  geneData <- data3 %>% filter(chrom==geneChr, Pos>geneStart-100000, Pos<geneEnd+100000) %>% na.omit()
  lead <- geneData %>% dplyr::arrange(P)
  leadPos <- lead$Pos[1]
  QTLdata <- geneData %>% filter(Pos>leadPos-50000, Pos<leadPos+50000) %>% na.omit()
  sameSNP <- intersect(QTLdata$SNP, GWASdataAll$SNP)
  QTLdata <- QTLdata[QTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
  GWASdata <- GWASdataAll[GWASdataAll$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
  coloc_data <- list(dataset1=list(snp=QTLdata$SNP,beta=QTLdata$beta,varbeta=QTLdata$varbeta,N=QTLdata$samplesize,MAF=QTLdata$MAF,z=QTLdata$z,pvalues=QTLdata$P,type="quant"),
                     dataset2=list(snp=GWASdata$SNP,beta=GWASdata$beta,varbeta=GWASdata$varbeta,N=GWASdata$samplesize,MAF=GWASdata$MAF,z=GWASdata$z,pvalues=GWASdata$P,type="cc"))
  result <- coloc.abf(dataset1=coloc_data$dataset1, dataset2=coloc_data$dataset2)
  outTab <- rbind(outTab, data.frame(ID=eqtlID, Symbol=geneName,
                                     PP.H0=result$summary[2],
                                     PP.H1=result$summary[3],
                                     PP.H2=result$summary[4],
                                     PP.H3=result$summary[5],
                                     PP.H4=result$summary[6]))
  GWAS_fn <- GWASdata[,c('SNP','P')] %>% dplyr::rename(rsid=SNP, pval=P)
  eQTL_fn <- QTLdata[,c('SNP','P')] %>% dplyr::rename(rsid=SNP, pval=P)
  pdf(paste0(eqtlID, ".", geneName, ".pdf"), width=8, height=6)
  print(locuscompare(in_fn1=GWAS_fn, in_fn2=eQTL_fn, snp=NULL, title1='GWAS', title2='eQTL'))
  dev.off()
}
write.csv(outTab, file="coloc.result.csv", row.names=F)