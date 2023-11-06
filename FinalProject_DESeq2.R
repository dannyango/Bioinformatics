setwd("E:/School UF/2023 Spring - MCB 6937 - Computational Genomics and Epigenomics/Final Project")

# Loading the package in R
library("DESeq2")
library("ggplot2")

# Read in counts data
control_rep1 <- read.delim("control_rep1_gene_summary", header=FALSE)
control_rep2 <- read.delim("control_rep2_gene_summary", header=FALSE)
cef_rep1 <- read.delim("cef_rep1_gene_summary", header=FALSE)
cef_rep2 <- read.delim("cef_rep2_gene_summary", header=FALSE)
chl_rep1 <- read.delim("chl_rep1_gene_summary", header=FALSE)
chl_rep2 <- read.delim("chl_rep2_gene_summary", header=FALSE)
#cip_rep1 <- read.delim("cip_rep1_gene_summary", header=FALSE)
#cip_rep2 <- read.delim("cip_rep2_gene_summary", header=FALSE)
#ery_rep1 <- read.delim("ery_rep1_gene_summary", header=FALSE)
#ery_rep2 <- read.delim("ery_rep2_gene_summary", header=FALSE)
#imi_rep1 <- read.delim("imi_rep1_gene_summary", header=FALSE)
#imi_rep2 <- read.delim("imi_rep2_gene_summary", header=FALSE)
kan_rep1 <- read.delim("kan_rep1_gene_summary", header=FALSE)
kan_rep2 <- read.delim("kan_rep2_gene_summary", header=FALSE)
#mit_rep1 <- read.delim("mit_rep1_gene_summary", header=FALSE)
#mit_rep2 <- read.delim("mit_rep2_gene_summary", header=FALSE)
#pol_rep1 <- read.delim("pol_rep1_gene_summary", header=FALSE)
#pol_rep2 <- read.delim("pol_rep2_gene_summary", header=FALSE)
tet_rep1 <- read.delim("tet_rep1_gene_summary", header=FALSE)
tet_rep2 <- read.delim("tet_rep2_gene_summary", header=FALSE)

counts <- read.delim("control_rep1_gene_summary", header=FALSE)
names(counts)[2] = "control_1"
counts$control_2 <- as.numeric(paste(control_rep2$V2))
counts$cef_1 <- as.numeric(paste(cef_rep1$V2))
counts$cef_2 <- as.numeric(paste(cef_rep2$V2))
counts$chl_1 <- as.numeric(paste(chl_rep1$V2))
counts$chl_2 <- as.numeric(paste(chl_rep2$V2))
#counts$cip_1 <- as.numeric(paste(cip_rep1$V2))
#counts$cip_2 <- as.numeric(paste(cip_rep2$V2))
#counts$ery_1 <- as.numeric(paste(ery_rep1$V2))
#counts$ery_2 <- as.numeric(paste(ery_rep2$V2))
#counts$imi_1 <- as.numeric(paste(imi_rep1$V2))
#counts$imi_2 <- as.numeric(paste(imi_rep2$V2))
counts$kan_1 <- as.numeric(paste(kan_rep1$V2))
counts$kan_2 <- as.numeric(paste(kan_rep2$V2))
#counts$mit_1 <- as.numeric(paste(mit_rep1$V2))
#counts$mit_2 <- as.numeric(paste(mit_rep2$V2))
#counts$pol_1 <- as.numeric(paste(pol_rep1$V2))
#counts$pol_2 <- as.numeric(paste(pol_rep2$V2))
counts$tet_1 <- as.numeric(paste(tet_rep1$V2))
counts$tet_2 <- as.numeric(paste(tet_rep2$V2))

rownames(counts) <- counts[,1]
gene_names <- counts[,1]
counts[,1] <- NULL

# Read in sample info
coldata <- read.delim("coldata_cef.txt", row.names=1)
coldata <- read.delim("coldata_chl.txt", row.names=1)
coldata <- read.delim("coldata_kan.txt", row.names=1)
coldata <- read.delim("coldata_tet.txt", row.names=1)

# Make sure the row names in colData match column name in counts
all(colnames(counts) %in% rownames(coldata))

# Make sure they are in the same order
all(colnames(counts) == rownames(coldata))

# Construct a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = ~ condition)

# Set the factor level
dds$condition <- relevel(dds$condition, ref = "untreated")

# Run DESeq
dds <- DESeq(dds)
#res <- results(dds)

# Summarize the results
res05 <- results(dds, alpha=0.05)
summary(res05)

# MA-plot
#plotMA(res05)

#Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef=2, type="apeglm")

cef_ma <- plotMA(resLFC, main = "CEF")
chl_ma <- plotMA(resLFC, main = "CHL")
kan_ma <- plotMA(resLFC, main = "KAN")
tet_ma <- plotMA(resLFC, main = "TET")

# this gives log2(n + 1)
ntd <- normTransform(dds)

# PCA analysis
vsd <- vst(dds, blind=FALSE)

cef_PCA <- plotPCA(vsd, intgroup=c("condition")) + ggtitle("CEF") + theme(plot.title = element_text(hjust = 0.5))
chl_PCA <- plotPCA(vsd, intgroup=c("condition")) + ggtitle("CHL") + theme(plot.title = element_text(hjust = 0.5))
kan_PCA <- plotPCA(vsd, intgroup=c("condition")) + ggtitle("KAN") + theme(plot.title = element_text(hjust = 0.5))
tet_PCA <- plotPCA(vsd, intgroup=c("condition")) + ggtitle("TET") + theme(plot.title = element_text(hjust = 0.5))

# Run heatmap
library("pheatmap")

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds))
df <- df[-c(2)]

cef_ph <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
chl_ph <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
kan_ph <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
tet_ph <- pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)

# Write the output
#write.csv(as.data.frame(res05), file="condition_treated_results.csv")

# Output genes with p < 0.05
#resSig <- subset(res05, padj < 0.05)
#write.csv(as.data.frame(resSig), file="condition_treated_results05.csv")

# Extract the normalized counts for expression levels
#dds.countsdata <- estimateSizeFactors(dds)
#NormalizedCount <- counts(dds.countsdata, normalized = T)
#write.csv(NormalizedCount, file="condition_treated_normcounts.csv", row.names = TRUE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
organism = "org.EcK12.eg.db"
library(organism, character.only = TRUE)
require(DOSE)

original_gene_list <- res05$log2FoldChange
names(original_gene_list) <- gene_names
gene_list <- na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
#keytypes(org.EcK12.eg.db)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ALIAS", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")

#====================================================================================

cef_dp <- dotplot(gse, showCategory=10, split=".sign") + ggtitle("CEF") + theme(plot.title = element_text(hjust = 0.5))
x2 <- pairwise_termsim(gse)
cef_em <- emapplot(x2) + ggtitle("CEF") + theme(plot.title = element_text(hjust = 0.5))
cef_rp <- ridgeplot(gse) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size = 8)) + ggtitle("CEF") + theme(plot.title = element_text(hjust = 0.5))

chl_dp <- dotplot(gse, showCategory=10, split=".sign") + ggtitle("CHL") + theme(plot.title = element_text(hjust = 0.5))
x2 <- pairwise_termsim(gse)
chl_em <- emapplot(x2) + ggtitle("CHL") + theme(plot.title = element_text(hjust = 0.5))
chl_rp <- ridgeplot(gse) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size = 8)) + ggtitle("CHL") + theme(plot.title = element_text(hjust = 0.5))

kan_dp <- dotplot(gse, showCategory=10, split=".sign") + ggtitle("KAN") + theme(plot.title = element_text(hjust = 0.5))
x2 <- pairwise_termsim(gse)
kan_em <- emapplot(x2) + ggtitle("KAN") + theme(plot.title = element_text(hjust = 0.5))
kan_rp <- ridgeplot(gse) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size = 8)) + ggtitle("KAN") + theme(plot.title = element_text(hjust = 0.5))

tet_dp <- dotplot(gse, showCategory=10, split=".sign") + ggtitle("TET") + theme(plot.title = element_text(hjust = 0.5))
x2 <- pairwise_termsim(gse)
tet_em <- emapplot(x2) + ggtitle("TET") + theme(plot.title = element_text(hjust = 0.5))
tet_rp <- ridgeplot(gse) + labs(x = "enrichment distribution") + theme(axis.text.y = element_text(size = 8)) + ggtitle("TET") + theme(plot.title = element_text(hjust = 0.5))

cef_PCA + chl_PCA + kan_PCA + tet_PCA
cef_ph + chl_ph + kan_ph + tet_ph
cef_dp + chl_dp + kan_dp + tet_dp
cef_em + chl_em + kan_em + tet_em
cef_rp + chl_rp + kan_rp + tet_rp
