cd /blue/mcb4934/share/danny.ngo #change directory
mkdir Assignment2 #make new folder
cd Assignment2 #change directory to newly created folder

module load bwa
module load picard
module load samtools
module load gatk
module load bcftools

#Question 2
bwa index Zm-B73_chr6.fa #create bwa index of reference genome
picard CreateSequenceDictionary -R Zm-B73_chr6.fa -O Zm-B73_chr6.dict #picard dictionary
samtools faidx Zm-B73_chr6.fa #samtools index

#Question 3, p1
bwa mem -M -t 16 -R "@RG\tID:4\tSM:P1\tLB:mop1\tPL:ILLUMINA" Zm-B73_chr6.fa p1_1.fq p1_2.fq > p1.sam #align fastq reads
samtools view -bS -q 10 p1.sam > p1_q10.bam #convert from .bam to .sam
picard SortSam -I p1_q10.bam -O p1_q10_sorted.bam -SO coordinate #sort alignment
picard MarkDuplicates --REMOVE_DUPLICATES true --CREATE_INDEX true -I p1_q10_sorted.bam -O p1_q10_sorted_dedup.bam --METRICS_FILE p1_marked_dedup_metrics.txt #remove duplicates from PCR
gatk HaplotypeCaller -I p1_q10_sorted_dedup.bam -R Zm-B73_chr6.fa -ERC GVCF -O p1.g.vcf #variant calling

#Question 3, w1
bwa mem -M -t 16 -R "@RG\tID:4\tSM:W1\tLB:mop1\tPL:ILLUMINA" Zm-B73_chr6.fa w1_1.fq w1_2.fq > w1.sam #align fastq reads
samtools view -bS -q 10 w1.sam > w1_q10.bam #convert from .bam to .sam
picard SortSam -I w1_q10.bam -O w1_q10_sorted.bam -SO coordinate #sort alignment
picard MarkDuplicates --REMOVE_DUPLICATES true --CREATE_INDEX true -I w1_q10_sorted.bam -O w1_q10_sorted_dedup.bam --METRICS_FILE w1_marked_dedup_metrics.txt #remove duplicates from PCR
gatk HaplotypeCaller -I w1_q10_sorted_dedup.bam -R Zm-B73_chr6.fa -ERC GVCF -O w1.g.vcf #variant calling

#Question 4
gatk GenomicsDBImport -V p1.g.vcf -V w1.g.vcf --genomicsdb-workspace-path Samples_db_chr6 --intervals chr6 #import GVCF files into GenomicsDB
gatk GenotypeGVCFs -R Zm-B73_chr6.fa -V gendb://Samples_db_chr6 -O Samples_variants.vcf #combine into single vcf file

#Question 5
bcftools filter -i 'QUAL>30' Samples_variants.vcf > Samples_variants_QUAL30.vcf #filter out quality <30
gatk VariantsToTable -V Samples_variants_QUAL30.vcf -F CHROM -F POS -F TYPE -GF AD -O Samples_chr6.table #data summary
grep -w -c "SNP" Samples_chr6.table #counts number of SNPs
grep -w -c "INDEL" Samples_chr6.table #counts number of INDELs

#Question 6
gatk SelectVariants -R Zm-B73_chr6.fa -V Samples_variants_QUAL30.vcf --select-type-to-include SNP -O Samples_variants_QUAL30_SNPs.vcf #filter for only SNPs
bcftools view Samples_variants_QUAL30_SNPs.vcf -o Samples_variants_QUAL30_SNPs.bcf #convert .vcf to .bcf
bcftools index Samples_variants_QUAL30_SNPs.bcf
bcftools view -r chr6:50000000-60000000 Samples_variants_QUAL30_SNPs.bcf > Samples_variants_QUAL30_SNPs_subset.bcf  #filter only in region
bcftools view Samples_variants_QUAL30_SNPs_subset.bcf -H | grep -w -c "chr6" #counts number of SNPs
# professor's solution
# gatk VariantsToTable -V Samples_variants_QUAL30_SNPs_subset.bcf -F CHROM -F POS -F TYPE -GF AD -O Samples_chr6_subset.table
# grep -w -c 'SNP' Samples_chr6_subset.table