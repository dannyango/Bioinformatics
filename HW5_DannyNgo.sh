#1
module load fastqc
fastqc -t 8 /blue/mcb4934/share/meixiazhao/Assignment5/*.fq.gz -o /blue/mcb4934/share/danny.ngo/Assignment5

module load trimmomatic
trimmomatic SE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment5/B73_H3K27ac_rep1.fq.gz \
B73_H3K27ac_rep1_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30

trimmomatic SE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment5/B73_H3K27ac_rep2.fq.gz \
B73_H3K27ac_rep2_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30

trimmomatic SE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment5/B73_Input_rep1.fq.gz \
B73_Input_rep1_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30

trimmomatic SE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment5/B73_Input_rep2.fq.gz \
B73_Input_rep2_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:30

#2
module load bowtie2
bowtie2-build Zm-B73_chr7.fa Zm-B73_chr7

#3
bowtie2 -x Zm-B73_chr7 -q B73_H3K27ac_rep1_trimmed.fq.gz -p 8 -k 5 -S B73_H3K27ac_rep1_trimmed.sam
bowtie2 -x Zm-B73_chr7 -q B73_H3K27ac_rep2_trimmed.fq.gz -p 8 -k 5 -S B73_H3K27ac_rep2_trimmed.sam
bowtie2 -x Zm-B73_chr7 -q B73_Input_rep1_trimmed.fq.gz -p 8 -k 5 -S B73_Input_rep1_trimmed.sam
bowtie2 -x Zm-B73_chr7 -q B73_Input_rep2_trimmed.fq.gz -p 8 -k 5 -S B73_Input_rep2_trimmed.sam

#4
module load samtools
samtools view -bq 1 B73_H3K27ac_rep1_trimmed.sam > B73_H3K27ac_rep1_trimmed_unique.bam
samtools view -bq 1 B73_H3K27ac_rep2_trimmed.sam > B73_H3K27ac_rep2_trimmed_unique.bam
samtools view -bq 1 B73_Input_rep1_trimmed.sam > B73_Input_rep1_trimmed_unique.bam
samtools view -bq 1 B73_Input_rep2_trimmed.sam > B73_Input_rep2_trimmed_unique.bam

#5
samtools sort -n -o B73_H3K27ac_rep1_trimmed_unique_sorted.bam B73_H3K27ac_rep1_trimmed_unique.bam
samtools sort -n -o B73_H3K27ac_rep2_trimmed_unique_sorted.bam B73_H3K27ac_rep2_trimmed_unique.bam
samtools sort -n -o B73_Input_rep1_trimmed_unique_sorted.bam B73_Input_rep1_trimmed_unique.bam
samtools sort -n -o B73_Input_rep2_trimmed_unique_sorted.bam B73_Input_rep2_trimmed_unique.bam

#6
module load macs
macs2 filterdup -i B73_H3K27ac_rep1_trimmed_unique_sorted.bam --keep-dup 1 -o B73_H3K27ac_rep1_trimmed_unique_sorted.bed
macs2 filterdup -i B73_H3K27ac_rep2_trimmed_unique_sorted.bam --keep-dup 1 -o B73_H3K27ac_rep2_trimmed_unique_sorted.bed
macs2 filterdup -i B73_Input_rep1_trimmed_unique_sorted.bam --keep-dup 1 -o B73_Input_rep1_trimmed_unique_sorted.bed
macs2 filterdup -i B73_Input_rep2_trimmed_unique_sorted.bam --keep-dup 1 -o B73_Input_rep2_trimmed_unique_sorted.bed

#7
macs2 callpeak -t B73_H3K27ac_rep1_trimmed_unique_sorted.bed -c B73_Input_rep1_trimmed_unique_sorted.bed -f AUTO -g 1.5e8 -q 0.05 --broad --outdir macs2_broad -n B73_H3K27ac_rep1
macs2 callpeak -t B73_H3K27ac_rep2_trimmed_unique_sorted.bed -c B73_Input_rep2_trimmed_unique_sorted.bed -f AUTO -g 1.5e8 -q 0.05 --broad --outdir macs2_broad -n B73_H3K27ac_rep2

#8
module load gcc/5.2.0
module load python
module load idr

mkdir macs2_broad
cd macs2_broad

sort -k8,8nr B73_H3K27ac_rep1_peaks.broadPeak > B73_H3K27ac_rep1_peaks_sorted.broadPeak
sort -k8,8nr B73_H3K27ac_rep2_peaks.broadPeak > B73_H3K27ac_rep2_peaks_sorted.broadPeak

#9
idr --samples B73_H3K27ac_rep1_peaks_sorted.broadPeak B73_H3K27ac_rep2_peaks_sorted.broadPeak \
--input-file-type broadPeak \
--rank p.value \
--output-file peaks-idr \
--plot \
--log-output-file peaks.idr.log

#10
awk -v OFS="\t" '{if ($7=="+") {print $1, $4-15000, $5, $7, $9} else {print $1, $4, $5+15000, $7, $9}}' Zm-B73_genes_chr7.gff3 > Zm-B73_genes_chr7.bed

module load bedtools
bedtools intersect -a /blue/mcb4934/share/danny.ngo/Assignment5/macs2_broad/B73_H3K27ac_rep1_peaks_sorted.broadPeak -b Zm-B73_genes_chr7.bed -wa -wb > Zm-B73_genes_chr7_peaks

#professor's answer
awk -v OFS="\t" '{if ($7 == "+") {print $1, $4-15000, $4-1, $7, $9} else {print $1, $5+1, $5+15000, $7, $9}}' Zm-B73_genes_chr7.gff3 > Zm-B73_genes_chr7.bed
bedtools intersect -a peaks-idr -b ../Zm-B73_genes_chr7.bed -wa -wb > peaks-idr_15kb_of_genes wc -l peaks-idr_15kb_of_genes

#11
rm *.sam
rm *.bam
