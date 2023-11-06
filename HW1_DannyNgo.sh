cd /blue/mcb4934/share/danny.ngo #change directory
mkdir Homework01 #make new folder
cd Homework01 #change directory to newly created folder

wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz #download gff3 file
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz #download fasta file

gunzip Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz #decompress file 1
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz #decompress file 2

grep -c ">" Zm-B73-REFERENCE-NAM-5.0.fa #counts number of sequences
# OR grep ">" Zm-B73-REFERENCE-NAM-5.0.fa | wc -l
grep -w -c "gene" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 #counts number of genes
# OR grep -w "gene" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | wc -l
grep -w -c "exon" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 #counts number of exons
# OR grep -w "gene" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | wc -l
grep -w "exon" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 > Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3_exon #create new file containing only exons
cut -f 1,3,4,5,7,9 Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3_exon > Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3_exon_part #create new file containing only specified colummns

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR22910370/SRR22910370 #download file
module load sra
fasterq-dump --split-files SRR22910370 #convert SRA to fastq; split into 2 files
module load fastqc
fastqc SRR22910370_1.fastq #check quality of file 1
fastqc SRR22910370_2.fastq #check quality of file 2
module load fastx_toolkit
fastq_to_fasta -i SRR22910370_1.fastq -o SRR22910370_1.fa #convert file 1 from fastq to fasta
fastq_to_fasta -i SRR22910370_2.fastq -o SRR22910370_2.fa #convert file 2 from fastq to fasta