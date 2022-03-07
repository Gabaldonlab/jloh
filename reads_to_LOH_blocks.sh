#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -J JLOH
#SBATCH --time=48:00:00
#SBATCH --error=reads_to_LOH_blocks.run.stderr
#SBATCH --output=reads_to_LOH_blocks.run.stdout

module purge
module load intel mkl/2017.4 impi/2017.4
module load gcc/7.2.0 hisat2/2.1.0 samtools/1.13 bcftools/1.6 htslib/1.8 bedtools/2.25.0 graphviz/2.40.1

SAMPLE_ID="bendixsen_data"

nextflow \
-C reads_to_LOH_blocks.config \
run \
reads_to_LOH_blocks.nf \
-resume \
-work-dir ${SAMPLE_ID}.work \
-with-dag ${SAMPLE_ID}.svg \
-with-report ${SAMPLE_ID}.report \
-with-timeline ${SAMPLE_ID}.timeline \
--input_data /gpfs/projects/bsc40/current/mschiavi/hybridiv/raw_data/public/association_table.tsv \
--reads_dir /gpfs/projects/bsc40/current/mschiavi/hybridiv/raw_data/public \
--threads 48 \
--run_id ${SAMPLE_ID} \
--output_dir /gpfs/projects/bsc40/current/mschiavi/hybridiv/public_data 
