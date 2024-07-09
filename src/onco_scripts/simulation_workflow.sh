#!/usr/bin/env sh

# Debugging
set -e
set -o pipefail

# Input parameters
divergence_rate=$1 # establish the divergence rate (0-1) in control sample
reference_genome=$2
mutation_rate=$3 # establish the mutation rate (0-1) in tumor after LOH
read_length=$4
paired_reads_number=$5
ID_control=$6
ID_tumor=$7
threads=$8

# Create haplotypes for diploid genome with divergence included
jloh sim --fasta $reference_genome --divergence $divergence_rate --out-fasta mut_A.fa --out-haplotypes mut.haplotypes_1.tsv --chrom-name-replace chr chr_A
jloh sim --fasta $reference_genome --divergence $divergence_rate --out-fasta mut_B.fa --out-haplotypes mut.haplotypes_2.tsv --chrom-name-replace chr chr_B


# Perform simulation for each haplotype without newly divergence
wgsim -N $paired_reads_number -1 $read_length -2 $read_length -r 0.0 -e 0.001 -R 0.05 mut_A.fa ${ID_control}_R1_A.fastq ${ID_control}_R2_A.fastq
wgsim -N $paired_reads_number -1 $read_length -2 $read_length -r 0.0 -e 0.001 -R 0.05 mut_B.fa ${ID_control}_R1_B.fastq ${ID_control}_R2_B.fastq

# Join haplotypes for each read
cat ${ID_control}_R1_A.fastq ${ID_control}_R1_B.fastq > ${ID_control}_R1.fastq
cat ${ID_control}_R2_A.fastq ${ID_control}_R2_B.fastq > ${ID_control}_R2.fastq

# Run BWA-MEM with specified index directory for control sample
bwa mem -R "@RG\tID:$ID_control\tSM:$ID_control\tLB:library1\tPL:ILLUMINA" -t ${threads} $reference_genome ${ID_control}_R1.fastq ${ID_control}_R2.fastq > aligned.$ID_control.sam 2> aligned.$ID_control.log

# Filter and compress the SAM file
samtools view -F 256 -F 0x4 -b -h -@${threads} aligned.$ID_control.sam > aligned.${ID_control}_f.bam 2> output_error_$ID_control.txt
samtools sort -@${threads} aligned.${ID_control}_f.bam > aligned.${ID_control}_fs.bam
samtools index aligned.${ID_control}_fs.bam

gatk HaplotypeCaller -R $reference_genome -I aligned.${ID_control}_fs.bam -O ${ID_control}_HC.vcf -A ExcessHet -A CountNs \
-A BaseQuality -A VariantType -A RMSMappingQuality -A IndelClassify -A IndelLength -A GcContent -A AlleleFraction -A Coverage \
-A FisherStrand -A MappingQuality -A MappingQualityZero > output_HC_${ID_control}.log 2>&1

# Filtering VCF files
bcftools view -i 'QUAL>=20 && VARIANT_TYPE=="snp" && FORMAT/DP>=4' ${ID_control}_HC.vcf > ${ID_control}_HC_filtered.vcf
bcftools view -i 'FORMAT/AF>=0.25 && FORMAT/AF<=0.75' -Oz -o ${ID_control}_HC_hetero.vcf.gz ${ID_control}_HC_filtered.vcf
bcftools index -t ${ID_control}_HC_hetero.vcf.gz

# Create truth sets for JLOH testing
# Bed from hetero SNPs in control
awk -v OFS='\t' '!/^#/ {print $1, $2-1, $2, $3}' ${ID_control}_HC_hetero.vcf.gz > ${ID_control}_HC_hetero.bed
bedtools merge -i ${ID_control}_HC_hetero.bed -d 500 -c 4 -o count > merged_intervals_${ID_control}.bed
# Maintain regions where at least 2 SNPs are present in the selected distance
awk '$4 >= 2' merged_intervals_${ID_control}.bed > intervals_with_2plus_snps_${ID_control}.bed

# Positive truth set
shuf -n 100 intervals_with_2plus_snps_${ID_control}.bed > true_positives_${ID_control}.bed
sort -k1,1 -k2,2n true_positives_${ID_control}.bed > true_positives_s_${ID_control}.bed

# Negative truth set
# Keep variants not selected to be positive truth set
bedtools intersect -v -a intervals_with_2plus_snps_${ID_control}.bed -b true_positives_${ID_control}.bed > candidate_negatives_${ID_control}.bed
# Select and sort 100 random negatives
shuf -n 100 candidate_negatives_${ID_control}.bed > true_negatives_${ID_control}.bed
sort -k1,1 -k2,2n true_negatives_${ID_control}.bed > true_negatives_${ID_control}_s.bed

# Keep only a VCF file from variants from the positive truth set (true positives)
bedtools intersect -header -a ${ID_control}_HC_hetero.vcf.gz -b true_positives_s_${ID_control}.bed > LOH_${ID_control}_HC.vcf

# Compress and index
bcftools view -Oz -o LOH_${ID_control}_HC.vcf.gz LOH_${ID_control}_HC.vcf
bcftools index -t LOH_${ID_control}_HC.vcf.gz

# Provoke LOH in control sample to simulate tumor genome: add variants as REF
cat mut_A.fa | bcftools consensus LOH_${ID_control}_HC.vcf.gz --haplotype R > mut_A_${ID_tumor}.fa
cat mut_B.fa | bcftools consensus LOH_${ID_control}_HC.vcf.gz --haplotype R > mut_B_${ID_tumor}.fa

# Perform simulation for each haplotype. Include divergence here to simulate mutations from cancer stage
wgsim -N $paired_reads_number -1 $read_length -2 $read_length -r $mutation_rate -e 0.001 -R 0.05 mut_A_${ID_tumor}.fa ${ID_tumor}_R1_A.fastq ${ID_tumor}_R2_A.fastq
wgsim -N $paired_reads_number -1 $read_length -2 $read_length -r $mutation_rate -e 0.001 -R 0.05 mut_B_${ID_tumor}.fa ${ID_tumor}_R1_B.fastq ${ID_tumor}_R2_B.fastq

# Join haplotypes for each read
cat ${ID_tumor}_R1_A.fastq ${ID_tumor}_R1_B.fastq > ${ID_tumor}_R1.fastq
cat ${ID_tumor}_R2_A.fastq ${ID_tumor}_R2_B.fastq > ${ID_tumor}_R2.fastq

# Workflow in tumoral sample
bwa mem -R "@RG\tID:$ID_tumor\tSM:$ID_tumor\tLB:library1\tPL:ILLUMINA" -t $threads "$reference_genome" ${ID_tumor}_R1.fastq ${ID_tumor}_R2.fastq > aligned.$ID_tumor.sam 2> aligned.$ID_tumor.log

# Filter and compress the SAM file
samtools view -F 256 -F 0x4 -b -h -@${threads} aligned.${ID_tumor}.sam > aligned.${ID_tumor}_f.bam 2> output_error_${ID_tumor}.txt
samtools sort -@${threads} aligned.${ID_tumor}_f.bam > aligned.${ID_tumor}_fs.bam
samtools index aligned.${ID_tumor}_fs.bam

# Variant calling: force-calling mode
gatk Mutect2 -R $reference_genome -I aligned.${ID_tumor}_fs.bam -alleles ${ID_control}_HC_hetero.vcf.gz -O ${ID_tumor}_FC.vcf -A ExcessHet -A CountNs -A BaseQuality \
-A VariantType -A RMSMappingQuality -A IndelClassify -A IndelLength -A AlleleFraction -A Coverage -A FisherStrand -A MappingQuality \
-A MappingQualityZero > output_FC_${ID_tumor}.log 2>&1

# Filtering
bcftools view -i 'FILTER=="PASS" && VARIANT_TYPE=="snp" && FORMAT/DP>=4' ${ID_tumor}_FC.vcf > ${ID_tumor}_FC_final.vcf
