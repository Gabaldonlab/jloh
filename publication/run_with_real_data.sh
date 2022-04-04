#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48

cd /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/publication

module load gcc/7.2.0 hisat2/2.1.0 samtools/1.13 bcftools/1.8 htslib/1.8 bedtools/2.29.2

nextflow \
-C run_with_real_data.config \
run \
run_with_real_data.nf \
-resume \
-work-dir /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/work-real_data \
-with-report /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/publication/run_with_real_data.nf.report \
--output_dir /gpfs/projects/bsc40/current/mschiavi/jloh/results_real-data
