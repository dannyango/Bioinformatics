#1
module load bismark
mkdir genome
cd genome
cp ../Zm-B73_chr7.fa Zm-B73_chr7.fa
cd ..
bismark_genome_preparation --parallel 8 ./genome

#2
module load fastqc
module load trimmomatic

fastqc -t 8 /blue/mcb4934/share/meixiazhao/Assignment4/*.fq.gz -o /blue/mcb4934/share/danny.ngo/Assignment4

trimmomatic PE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_WT_rep1_R1.fq.gz \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_WT_rep1_R2.fq.gz \
WGBS_WT_rep1_R1_paired_trimmed.fq.gz WGBS_WT_rep1_R1_unpaired_trimmed.fq.gz \
WGBS_WT_rep1_R2_paired_trimmed.fq.gz WGBS_WT_rep1_R2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_WT_rep2_R1.fq.gz \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_WT_rep2_R2.fq.gz \
WGBS_WT_rep2_R1_paired_trimmed.fq.gz WGBS_WT_rep2_R1_unpaired_trimmed.fq.gz \
WGBS_WT_rep2_R2_paired_trimmed.fq.gz WGBS_WT_rep2_R2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_mutant_rep1_R1.fq.gz \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_mutant_rep1_R2.fq.gz \
WGBS_mutant_rep1_R1_paired_trimmed.fq.gz WGBS_mutant_rep1_R1_unpaired_trimmed.fq.gz \
WGBS_mutant_rep1_R2_paired_trimmed.fq.gz WGBS_mutant_rep1_R2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

trimmomatic PE -phred33 -threads 8 \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_mutant_rep2_R1.fq.gz \
/blue/mcb4934/share/meixiazhao/Assignment4/WGBS_mutant_rep2_R2.fq.gz \
WGBS_mutant_rep2_R1_paired_trimmed.fq.gz WGBS_mutant_rep2_R1_unpaired_trimmed.fq.gz \
WGBS_mutant_rep2_R2_paired_trimmed.fq.gz WGBS_mutant_rep2_R2_unpaired_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:75

#3
module load bismark

bismark ./genome -I 50 -N 1 --multicore 8 -q -1 WGBS_WT_rep1_R1_paired_trimmed.fq.gz -2 WGBS_WT_rep1_R2_paired_trimmed.fq.gz
bismark ./genome -I 50 -N 1 --multicore 8 -q -1 WGBS_WT_rep2_R1_paired_trimmed.fq.gz -2 WGBS_WT_rep2_R2_paired_trimmed.fq.gz
bismark ./genome -I 50 -N 1 --multicore 8 -q -1 WGBS_mutant_rep1_R1_paired_trimmed.fq.gz -2 WGBS_mutant_rep1_R2_paired_trimmed.fq.gz
bismark ./genome -I 50 -N 1 --multicore 8 -q -1 WGBS_mutant_rep2_R1_paired_trimmed.fq.gz -2 WGBS_mutant_rep2_R2_paired_trimmed.fq.gz

#4
deduplicate_bismark --bam -p -o WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.bam WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.bam
deduplicate_bismark --bam -p -o WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.bam WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.bam
deduplicate_bismark --bam -p -o WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.bam WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.bam
deduplicate_bismark --bam -p -o WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.bam WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.bam

#5
bismark2report WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam
bismark2report WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam
bismark2report WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam
bismark2report WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam

#6
bismark_methylation_extractor -p --multicore 8 --gzip WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam
bismark_methylation_extractor -p --multicore 8 --gzip WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam
bismark_methylation_extractor -p --multicore 8 --gzip WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam
bismark_methylation_extractor -p --multicore 8 --gzip WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.bam

#7
bismark2bedGraph --CX -o WGBS_WT_rep1 \
CHG_OB_WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHG_OT_WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OB_WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OT_WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OB_WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OT_WGBS_WT_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz

bismark2bedGraph --CX -o WGBS_WT_rep2 \
CHG_OB_WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHG_OT_WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OB_WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OT_WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OB_WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OT_WGBS_WT_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz

bismark2bedGraph --CX -o WGBS_mutant_rep1 \
CHG_OB_WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHG_OT_WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OB_WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OT_WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OB_WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OT_WGBS_mutant_rep1_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz

bismark2bedGraph --CX -o WGBS_mutant_rep2 \
CHG_OB_WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHG_OT_WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OB_WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CHH_OT_WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OB_WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz \
CpG_OT_WGBS_mutant_rep2_R1_paired_trimmed_bismark_bt2_pe.deduplicated.txt.gz

#8
coverage2cytosine --CX --genome_folder ./genome -o WGBS_WT_rep1_cytocine WGBS_WT_rep1.gz.bismark.cov.gz
coverage2cytosine --CX --genome_folder ./genome -o WGBS_WT_rep2_cytocine WGBS_WT_rep2.gz.bismark.cov.gz
coverage2cytosine --CX --genome_folder ./genome -o WGBS_mutant_rep1_cytocine WGBS_mutant_rep1.gz.bismark.cov.gz
coverage2cytosine --CX --genome_folder ./genome -o WGBS_mutant_rep2_cytocine WGBS_mutant_rep2.gz.bismark.cov.gz

#9
python metilene_calculate_proportion.py WGBS_WT_rep1_cytocine.CX_report.txt WGBS_WT_rep1_cytocine.CX_report_pro
python metilene_calculate_proportion.py WGBS_WT_rep2_cytocine.CX_report.txt WGBS_WT_rep2_cytocine.CX_report_pro
python metilene_calculate_proportion.py WGBS_mutant_rep1_cytocine.CX_report.txt WGBS_mutant_rep1_cytocine.CX_report_pro
python metilene_calculate_proportion.py WGBS_mutant_rep2_cytocine.CX_report.txt WGBS_mutant_rep2_cytocine.CX_report_pro

#10
python metilene_merge_samples.py \
WGBS_WT_rep1_cytocine.CX_report_pro \
WGBS_WT_rep2_cytocine.CX_report_pro \
WGBS_mutant_rep1_cytocine.CX_report_pro \
WGBS_mutant_rep2_cytocine.CX_report_pro \
All_4_samples_cytosines

python metilene_split_cytosines.py \
All_4_samples_cytosines \
All_4_samples_cytosines_CG \
All_4_samples_cytosines_CHG \
All_4_samples_cytosines_CHH

#11
module load metilene
metilene -a g1 -b g2 -t 8 All_4_samples_cytosines_CG > All_4_samples_cytosines_CG_DMRs
metilene -a g1 -b g2 -t 8 All_4_samples_cytosines_CHG > All_4_samples_cytosines_CHG_DMRs
metilene -a g1 -b g2 -t 8 All_4_samples_cytosines_CHH > All_4_samples_cytosines_CHH_DMRs

#12
awk '{if ($5>0) print $5}' All_4_samples_cytosines_CG_DMRs | wc -l
awk '{if ($5<0) print $5}' All_4_samples_cytosines_CG_DMRs | wc -l

#13
module load bedtools
grep -v "#" Zm-B73_genes_chr7.gff3 | awk '{print $1"\t"$4-1"\t"$5-1"\t"$9}' > Zm-B73_genes_chr7.bed
sort -k1,1 -k2,2n Zm-B73_genes_chr7.bed > Zm-B73_genes_chr7_sorted.bed
bedtools intersect -a All_4_samples_cytosines_CG_DMRs.bed -b Zm-B73_genes_chr7_sorted.bed -wa -wb > CG_DMRs_in_gene_bodies
awk -v OFS="\t" '{print $1, $2, $3}' All_4_samples_cytosines_CG_DMRs > All_4_samples_cytosines_CG_DMRs.bed
grep "chr7" CG_DMRs_in_gene_bodies | wc -l


#14
module load bowtie

bowtie-build -f Zm-B73_chr7.fa Zm-B73_chr7

#15
module load trimmomatic

trimmomatic SE -phred33 -threads 16 /blue/mcb4934/share/meixiazhao/Assignment4/sRNA-seq_WT_rep1.fq.gz \
sRNA-seq_WT_rep1_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15

trimmomatic SE -phred33 -threads 16 /blue/mcb4934/share/meixiazhao/Assignment4/sRNA-seq_WT_rep2.fq.gz \
sRNA-seq_WT_rep2_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15

trimmomatic SE -phred33 -threads 16 /blue/mcb4934/share/meixiazhao/Assignment4/sRNA-seq_mutant_rep1.fq.gz \
sRNA-seq_mutant_rep1_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15

trimmomatic SE -phred33 -threads 16 /blue/mcb4934/share/meixiazhao/Assignment4/sRNA-seq_mutant_rep2.fq.gz \
sRNA-seq_mutant_rep2_trimmed.fq.gz \
ILLUMINACLIP:${HPC_TRIMMOMATIC_ADAPTER}/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:15

#16
bowtie -v 0 Zm-B73_chr7 -p 4 -S -k 20 -q sRNA-seq_WT_rep1_trimmed.fq.gz > sRNA-seq_WT_rep1_trimmed.sam
bowtie -v 0 Zm-B73_chr7 -p 4 -S -k 20 -q sRNA-seq_WT_rep2_trimmed.fq.gz > sRNA-seq_WT_rep2_trimmed.sam
bowtie -v 0 Zm-B73_chr7 -p 4 -S -k 20 -q sRNA-seq_mutant_rep1_trimmed.fq.gz > sRNA-seq_mutant_rep1_trimmed.sam
bowtie -v 0 Zm-B73_chr7 -p 4 -S -k 20 -q sRNA-seq_mutant_rep2_trimmed.fq.gz > sRNA-seq_mutant_rep2_trimmed.sam

#17
module load samtools

samtools view -bS sRNA-seq_WT_rep1_trimmed.sam > sRNA-seq_WT_rep1_trimmed.bam
samtools view -bS sRNA-seq_WT_rep2_trimmed.sam > sRNA-seq_WT_rep2_trimmed.bam
samtools view -bS sRNA-seq_mutant_rep1_trimmed.sam > sRNA-seq_mutant_rep1_trimmed.bam
samtools view -bS sRNA-seq_mutant_rep2_trimmed.sam > sRNA-seq_mutant_rep2_trimmed.bam

#18
samtools sort -n -o sRNA-seq_WT_rep1_trimmed_sorted.bam sRNA-seq_WT_rep1_trimmed.bam
samtools sort -n -o sRNA-seq_WT_rep2_trimmed_sorted.bam sRNA-seq_WT_rep2_trimmed.bam
samtools sort -n -o sRNA-seq_mutant_rep1_trimmed_sorted.bam sRNA-seq_WT_rep1_trimmed.bam
samtools sort -n -o sRNA-seq_mutant_rep2_trimmed_sorted.bam sRNA-seq_WT_rep2_trimmed.bam

#19
awk -v OFS="\t" '{print $1, ".", "DMR", $2, $3, ".", ".", ".", "ID="$1"_"$2"_"$3}' All_4_samples_cytosines_CG_DMRs > All_4_samples_cytosines_CG_DMRs.gff3

module load htseq

htseq-count -f bam -m intersection-nonempty -s no -t DMR -i ID -o WT_sRNA-seq_WT_rep1_counts \
sRNA-seq_WT_rep1_trimmed_sorted.bam All_4_samples_cytosines_CG_DMRs.gff3 > sRNA-seq_WT_rep1_summary.txt

#20
awk -v OFS="\t" '{if ($2>10) {print $1, $2}}' sRNA-seq_WT_rep1_summary.txt | grep "chr7" | wc -l

#21
rm *_cytocine.CX_report.txt
rm *.sam

#22
We compared two replicates of the mutant sample to two replicates of the wild type sample. First, we found the DMRs are regions that have either increased or decreased methylation patterns. When genes are methylated, generally those genes would be inactivated and lead to decreased RNA transcription. In this assignment, we found more genes that were hypomethylated than were hypermethylated. We specifically analyzed only CG sequences, as opposed to CHG and CHH sites.

We also compared the differentially methylated regions (DMRs) that were mapped to the genes encoding small RNA. We just looked at the small RNAs for one of the wild type replicates due to time constraints. The above analyses ultimately found around 47 small RNAs that were differentially methylated. While we did find more DMRs that were hypomethylated, it would be interesting to see if the DMRs with the higher small RNA counts were hypomethylated as well. 