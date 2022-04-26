#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48

WD="/gpfs/projects/bsc40/current/mschiavi/jloh"
cd ${WD}/scripts/wf-real_data

module load gcc/7.2.0 hisat2/2.1.0 samtools/1.13 bcftools/1.8 htslib/1.8 bedtools/2.29.2

nextflow \
-C run_with_real_data.config \
run \
run_with_real_data.nf \
-resume \
-work-dir ${WD}/scripts/wf-real_data/work \
-with-report ${WD}/scripts/wf-real_data/run_with_real_data.nf.report \
--adapters /gpfs/projects/bsc40/mschiavi/software/trimmomatic/0.39/adapters/TruSeq2_and_3.PE.fa \
--output_dir ${WD}/results-real_data \
--input_data ${WD}/raw_data/Bendixsen_2021/Z.metadata.tsv
