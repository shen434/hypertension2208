smr-1.3.1-win.exe --bfile ./1000G.EUR/1000G_EUR --gwas-summary outcome.ma --beqtl-summary ./cis-eQTL-SMR_20191212/eQTLGen --maf 0.01 --out SMR.result --thread-num 10
smr-1.3.1-win.exe --bfile ./1000G.EUR/1000G_EUR --gwas-summary outcome.ma --beqtl-summary ./cis-eQTL-SMR_20191212/eQTLGen --out myplot --plot --probe ENSG00000107798 --probe-wind 500 --gene-list glist_hg19_refseq.txt


library(magick)
library(TeachingDemos)

inputFile="myplot.ENSG00000114573.txt"        
source("plot_SMR.r") 

SMRData = ReadSMRData(inputFile)
pdf(file="ATP6V1A.SMRLocusPlot.pdf", width=10, height=7)
SMRLocusPlot(data=SMRData, smr_thresh=0.001, heidi_thresh=0.05, plotWindow=1000, max_anno_probe=16)
dev.off()
pdf(file="ATP6V1A.SMREffectPlot.pdf", width=7, height=5.5)
SMREffectPlot(data=SMRData, trait_name="") 
dev.off()