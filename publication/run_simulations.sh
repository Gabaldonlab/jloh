#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48

cd /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/wf-sim_data

module load gcc/7.2.0 hisat2/2.1.0 samtools/1.13 bcftools/1.8 htslib/1.8 bedtools/2.29.2

nextflow \
-C run_tests.config \
run \
run_tests.nf \
-resume \
-work-dir /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/wf-sim_data/work \
-with-report /gpfs/projects/bsc40/current/mschiavi/jloh/scripts/run_tests.nf.report \
--ref_genome /gpfs/projects/bsc40/current/mschiavi/jloh/raw_data/S_cerevisiae.fa \
--output_dir /gpfs/projects/bsc40/current/mschiavi/jloh/results-sim_data
