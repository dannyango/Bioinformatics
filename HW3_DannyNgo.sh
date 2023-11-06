#1
module load fastqc
fastqc -t 16 control_rep1_1.fq.gz
fastqc -t 16 control_rep1_2.fq.gz

#2
module load trimmomatic

trimmomatic PE -phred33 -threads 16 control_rep1_1.fq.gz control_rep1_2.fq.gz \
control_rep1_1_paired_trimmed.fq.gz control_rep1_1_unpaired_trimmed.fq.gz \
control_rep1_2_paired_trimmed.fq.gz control_rep1_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 control_rep2_1.fq.gz control_rep2_2.fq.gz \
control_rep2_1_paired_trimmed.fq.gz control_rep2_1_unpaired_trimmed.fq.gz \
control_rep2_2_paired_trimmed.fq.gz control_rep2_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 control_rep3_1.fq.gz control_rep3_2.fq.gz \
control_rep3_1_paired_trimmed.fq.gz control_rep3_1_unpaired_trimmed.fq.gz \
control_rep3_2_paired_trimmed.fq.gz control_rep3_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 treatment_rep1_1.fq.gz treatment_rep1_2.fq.gz \
treatment_rep1_1_paired_trimmed.fq.gz treatment_rep1_1_unpaired_trimmed.fq.gz \
treatment_rep1_2_paired_trimmed.fq.gz treatment_rep1_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 treatment_rep2_1.fq.gz treatment_rep2_2.fq.gz \
treatment_rep2_1_paired_trimmed.fq.gz treatment_rep2_1_unpaired_trimmed.fq.gz \
treatment_rep2_2_paired_trimmed.fq.gz treatment_rep2_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 16 treatment_rep3_1.fq.gz treatment_rep3_2.fq.gz \
treatment_rep3_1_paired_trimmed.fq.gz treatment_rep3_1_unpaired_trimmed.fq.gz \
treatment_rep3_2_paired_trimmed.fq.gz treatment_rep3_2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

#3
module load hisat2
hisat2-build -f Mus_musculus.dna.fa Mus_musculus.dna

#4
hisat2 -p 16 -x Mus_musculus.dna -1 control_rep1_1_paired_trimmed.fq.gz -2 control_rep1_2_paired_trimmed.fq.gz -S control_rep1.sam
hisat2 -p 16 -x Mus_musculus.dna -1 control_rep2_1_paired_trimmed.fq.gz -2 control_rep2_2_paired_trimmed.fq.gz -S control_rep2.sam
hisat2 -p 16 -x Mus_musculus.dna -1 control_rep3_1_paired_trimmed.fq.gz -2 control_rep3_2_paired_trimmed.fq.gz -S control_rep3.sam
hisat2 -p 16 -x Mus_musculus.dna -1 treatment_rep1_1_paired_trimmed.fq.gz -2 treatment_rep1_2_paired_trimmed.fq.gz -S treatment_rep1.sam
hisat2 -p 16 -x Mus_musculus.dna -1 treatment_rep2_1_paired_trimmed.fq.gz -2 treatment_rep2_2_paired_trimmed.fq.gz -S treatment_rep2.sam
hisat2 -p 16 -x Mus_musculus.dna -1 treatment_rep3_1_paired_trimmed.fq.gz -2 treatment_rep3_2_paired_trimmed.fq.gz -S treatment_rep3.sam

#5
module load samtools
samtools view -bS control_rep1.sam > control_rep1.bam
samtools view -bS control_rep2.sam > control_rep2.bam
samtools view -bS control_rep3.sam > control_rep3.bam
samtools view -bS treatment_rep1.sam > treatment_rep1.bam
samtools view -bS treatment_rep2.sam > treatment_rep2.bam
samtools view -bS treatment_rep3.sam > treatment_rep3.bam

#6
samtools sort -n -o control_rep1_sorted.bam control_rep1.bam 
samtools sort -n -o control_rep2_sorted.bam control_rep2.bam
samtools sort -n -o control_rep3_sorted.bam control_rep3.bam 
samtools sort -n -o treatment_rep1_sorted.bam treatment_rep1.bam 
samtools sort -n -o treatment_rep2_sorted.bam treatment_rep2.bam 
samtools sort -n -o treatment_rep3_sorted.bam treatment_rep3.bam 

#7
module load htseq
htseq-count -f bam -m intersection-nonempty -s no -t gene -i ID -o control_rep1_counts control_rep1_sorted.bam Mus_musculus.genes.gff3 > control_rep1_gene_summary
htseq-count -f bam -m intersection-nonempty -s no -t gene -i ID -o control_rep2_counts control_rep2_sorted.bam Mus_musculus.genes.gff3 > control_rep2_gene_summary
htseq-count -f bam -m intersection-nonempty -s no -t gene -i ID -o control_rep3_counts control_rep3_sorted.bam Mus_musculus.genes.gff3 > control_rep3_gene_summary
htseq-count -f bam -m intersection-nonempty -s no -t gene -i ID -o treatment_rep1_counts treatment_rep1_sorted.bam Mus_musculus.genes.gff3 > treatment_rep1_gene_summary
htseq-count -f bam -m intersection-nonempty -s no -t gene -i ID -o treatment_rep2_counts treatment_rep2_sorted.bam Mus_musculus.genes.gff3 > treatment_rep2_gene_summary
htseq-count -f bam -m intersection-nonempty -s no -t gene -i ID -o treatment_rep3_counts treatment_rep3_sorted.bam Mus_musculus.genes.gff3 > treatment_rep3_gene_summary

#8
rm control_rep1_counts
rm control_rep2_counts
rm control_rep3_counts
rm treatment_rep1_counts
rm treatment_rep2_counts
rm treatment_rep3_counts

rm control_rep1.sam
rm control_rep2.sam
rm control_rep3.sam
rm treatment_rep1.sam
rm treatment_rep2.sam
rm treatment_rep3.sam

rm control_rep1.bam
rm control_rep2.bam
rm control_rep3.bam
rm treatment_rep1.bam
rm treatment_rep2.bam
rm treatment_rep3.bam

rm control_rep1_sorted.bam
rm control_rep2_sorted.bam
rm control_rep3_sorted.bam
rm treatment_rep1_sorted.bam
rm treatment_rep2_sorted.bam
rm treatment_rep3_sorted.bam

#9
setwd("C:/Users/Danny/Desktop/Assignment3")

#10
library(DESeq2)
library(pheatmap)

#11
countsdata <- read.delim("mice_counts.txt", row.names="gene_id")
coldata <- read.delim("sample_info.txt", row.names=1)

#12
all(colnames(countsdata) %in% rownames(coldata))
all(colnames(countsdata) == rownames(coldata))

#13
dds <- DESeqDataSetFromMatrix(
  countData = countsdata,
  colData = coldata,
  design = ~ condition)

#14
dds$condition <- relevel(dds$condition, ref = "untreated")

#15
keep <- rowSums(counts(dds)) >= 10
dds_keep <- dds[keep,]

#16
dds_run <- DESeq(dds_keep)
res <- results(dds_run)

#17
res1 <- results(dds_run, alpha=0.1)
summary(res1)

res05 <- results(dds_run, alpha=0.05)
summary(res05)

#18
plotMA(res)
plotMA(res1)
plotMA(res05)

#19
ntd <- normTransform(dds_run)
library(pheatmap)
select <- order(rowMeans(counts(dds_run,normalized=TRUE)), decreasing=TRUE)[1:30]
df <- as.data.frame(colData(dds)[,c("species", "condition")])
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=df)

#20
resSig <- subset(res, padj < 0.05)
write.csv(as.data.frame(resSig), file="condition_treated_results05.csv")

dds_run.countsdata <- estimateSizeFactors(dds)
NormalizedCount <- counts(dds_run.countsdata, normalized = T)
write.csv(NormalizedCount, file="condition_treated_normcounts.csv", row.names = TRUE)

#21
# term_id	term_name
GO:0005515	protein binding
GO:0015318	inorganic molecular entity transmembrane transporter activity
GO:0016917	GABA receptor activity
GO:0099106	ion channel regulator activity
GO:0015467	G-protein activated inward rectifier potassium channel activity
GO:0048731	system development
GO:0007267	cell-cell signaling
GO:0099601	regulation of neurotransmitter receptor activity
GO:0044057	regulation of system process
GO:0051128	regulation of cellular component organization
GO:0099173	postsynapse organization
GO:0035418	protein localization to synapse
GO:0007612	learning
GO:0098662	inorganic cation transmembrane transport
GO:0045202	synapse
GO:0034702	ion channel complex
GO:0098850	extrinsic component of synaptic vesicle membrane
GO:0043235	receptor complex

#genes queried
ENSMUSG00000003273
ENSMUSG00000003279
ENSMUSG00000005583
ENSMUSG00000005958
ENSMUSG00000010797
ENSMUSG00000010803
ENSMUSG00000019889
ENSMUSG00000019960
ENSMUSG00000020532
ENSMUSG00000020635
ENSMUSG00000021182
ENSMUSG00000021314
ENSMUSG00000021587
ENSMUSG00000021702
ENSMUSG00000021820
ENSMUSG00000022103
ENSMUSG00000022231
ENSMUSG00000022307
ENSMUSG00000022523
ENSMUSG00000025582
ENSMUSG00000026163
ENSMUSG00000026463
ENSMUSG00000026824
ENSMUSG00000026833
ENSMUSG00000027004
ENSMUSG00000027273
ENSMUSG00000027489
ENSMUSG00000027965
ENSMUSG00000028132
ENSMUSG00000028528
ENSMUSG00000028546
ENSMUSG00000028634
ENSMUSG00000028926
ENSMUSG00000029189
ENSMUSG00000029335
ENSMUSG00000029608
ENSMUSG00000029659
ENSMUSG00000031144
ENSMUSG00000031343
ENSMUSG00000032224
ENSMUSG00000032452
ENSMUSG00000032946
ENSMUSG00000033316
ENSMUSG00000034098
ENSMUSG00000034295
ENSMUSG00000034891
ENSMUSG00000036578
ENSMUSG00000036617
ENSMUSG00000037217
ENSMUSG00000037984
ENSMUSG00000038244
ENSMUSG00000038255
ENSMUSG00000038473
ENSMUSG00000038738
ENSMUSG00000038807
ENSMUSG00000039809
ENSMUSG00000041670
ENSMUSG00000041773
ENSMUSG00000043004
ENSMUSG00000043301
ENSMUSG00000044847
ENSMUSG00000045763
ENSMUSG00000046321
ENSMUSG00000046442
ENSMUSG00000046480
ENSMUSG00000050272
ENSMUSG00000050663
ENSMUSG00000051359
ENSMUSG00000052459
ENSMUSG00000052726
ENSMUSG00000052981
ENSMUSG00000054459
ENSMUSG00000054640
ENSMUSG00000058153
ENSMUSG00000059003
ENSMUSG00000059173
ENSMUSG00000059187
ENSMUSG00000060636
ENSMUSG00000063415
ENSMUSG00000067879
ENSMUSG00000072825
ENSMUSG00000087141
ENSMUSG00000094840
