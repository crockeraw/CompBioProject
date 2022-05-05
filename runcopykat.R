library(copykat)


gliob_data<-read.table("data/GSE84465_GBM_All_data.csv",header = T)

# cell.line sets the analysis to only have diploid or aneuploid
# n.cores should be set to 1 unless explicitly needing multiple cores, I use 2 on the longleaf server but only one on my computer

copykat.test <- copykat(rawmat=gliob_data,cell.line = "no",n.cores = 2,,sam.name = "GSE84465",genome = "hg20")

pred.test <- data.frame(copykat.test$prediction)
CNA.test <- data.frame(copykat.test$CNAmat)

CNVs<- read.table("_copykat_CNA_raw_results_gene_by_cell.txt")
counts<- CNVs[-1,c(8:ncol(CNVs))]
counts<- sapply(counts,as.numeric )
colnames(counts)<- CNVs[1,c(8:ncol(CNVs))]
rownames(counts)<- CNVs[-1,6]

colnames(counts)<- CNVs[1,c(6,8:ncol(CNVs))]
rownames(counts)<- CNVs[-1,6]

write.csv(counts,"CNA_matrix.csv")

means<- rowMeans(counts)
raw<- as.matrix(as.numeric(as.matrix(counts[-1])))

max(counts)
fulldf<- read.table("_copykat_CNA_results.txt")
