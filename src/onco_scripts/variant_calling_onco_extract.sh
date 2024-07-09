#!/usr/bin/env sh

# Debugging
set -e
set -o pipefail

# BASH script for variant calling in control and tumoral sample, after filtering with trimmomatics
ID_control=$1
ID_tumor=$2
reference_genome=$3
read_control_1=$4
read_control_2=$5
read_tumor_1=$6
read_tumor_2=$7
known_sites=$8 # for BQSR, such as dbSNP151 database
threads=$9

# Workflow in control sample
bwa mem -R "@RG\tID:$ID_control\tSM:$ID_control\tLB:library1\tPL:ILLUMINA" -t $threads "$reference_genome" "$read_control_1" "$read_control_2" > aligned.$ID_control.sam 2> aligned.$ID_control.log

# Filter and compress the SAM file
samtools view -F 256 -F 0x4 -b -h -@$threads aligned.$ID_control.sam > aligned.${ID_control}_f.bam 2> output_error_$ID_control.txt
samtools sort -@$threads aligned.${ID_control}_f.bam > aligned.${ID_control}_fs.bam
samtools index aligned.${ID_control}_fs.bam
picard MarkDuplicates I=aligned.${ID_control}_fs.bam REMOVE_DUPLICATES=true O=aligned.${ID_control}_fs_dedup.bam M=metrics.${ID_control}_fs_dedup.txt 2> markduplicates_$ID_control.log

# BQSR
gatk BaseRecalibrator -I aligned.${ID_control}_fs_dedup.bam -R $reference_genome --known-sites $known_sites -O recal_${ID_control}.table
gatk ApplyBQSR -R $reference_genome -I aligned.${ID_control}_fs_dedup.bam --bqsr-recal-file recal_${ID_control}.table -O aligned.${ID_control}_fs_dedup_recal.bam 2> applybqsr_$ID_control.log

# Variant calling
gatk HaplotypeCaller -R $reference_genome -I aligned.${ID_control}_fs_dedup_recal.bam -O ${ID_control}_HC.vcf -A ExcessHet -A CountNs \
-A BaseQuality -A VariantType -A RMSMappingQuality -A IndelClassify -A IndelLength -A GcContent -A AlleleFraction -A Coverage \
-A FisherStrand -A MappingQuality -A MappingQualityZero > output_HC_${ID_control}.log 2>&1

# Filtering VCF files
bcftools view -i 'QUAL>=20 && VARIANT_TYPE=="snp" && FORMAT/DP>=10' ${ID_control}_HC.vcf > ${ID_control}_HC_filtered.vcf
bcftools view -i 'FORMAT/AF>=0.25 && FORMAT/AF<=0.75' -Oz -o ${ID_control}_HC_hetero.vcf.gz ${ID_control}_HC_filtered.vcf
bcftools index -t ${ID_control}_HC_hetero.vcf.gz

# Workflow in tumoral sample
bwa mem -R "@RG\tID:$ID_tumor\tSM:$ID_tumor\tLB:library1\tPL:ILLUMINA" -t $threads "$reference_genome" "$read_tumor_1" "$read_tumor_2" > aligned.$ID_tumor.sam 2> aligned.$ID_tumor.log

# Filter and compress the SAM file
samtools view -F 256 -F 0x4 -b -h -@$threads aligned.${ID_tumor}.sam > aligned.${ID_tumor}_f.bam 2> output_error_${ID_tumor}.txt
samtools sort -@$threads aligned.${ID_tumor}_f.bam > aligned.${ID_tumor}_fs.bam
samtools index aligned.${ID_tumor}_fs.bam
picard MarkDuplicates I=aligned.${ID_tumor}_fs.bam REMOVE_DUPLICATES=true O=aligned.${ID_tumor}_fs_dedup.bam M=metrics.${ID_tumor}_fs_dedup.txt 2> markduplicates_$ID_tumor.log

# BQSR
gatk BaseRecalibrator -I aligned.${ID_tumor}_fs_dedup.bam -R $reference_genome --known-sites $known_sites -O recal_${ID_tumor}.table
gatk ApplyBQSR -R $reference_genome -I aligned.${ID_tumor}_fs_dedup.bam --bqsr-recal-file recal_${ID_tumor}.table -O aligned.${ID_tumor}_fs_dedup_recal.bam 2> applybqsr_$ID_tumor.log

# Variant calling: force-calling mode
gatk Mutect2 -R $reference_genome -I aligned.${ID_tumor}_fs_dedup_recal.bam -alleles ${ID_control}_HC_hetero.vcf.gz -O ${ID_tumor}_FC.vcf -A ExcessHet -A CountNs -A BaseQuality \
-A VariantType -A RMSMappingQuality -A IndelClassify -A IndelLength -A AlleleFraction -A Coverage -A FisherStrand -A MappingQuality \
-A MappingQualityZero > output_FC_${ID_tumor}.log 2>&1

# Index and filter the VCF file
gatk IndexFeatureFile -I ${ID_tumor}_FC.vcf
gatk FilterMutectCalls --min-reads-per-strand 4 -R $reference_genome -V ${ID_tumor}_FC.vcf -O ${ID_tumor}_FC_filtered.vcf

# Further filtering
bcftools view -i 'FILTER=="PASS" && VARIANT_TYPE=="snp" && FORMAT/DP>=10' ${ID_tumor}_FC_filtered.vcf > ${ID_tumor}_FC_final.vcf
