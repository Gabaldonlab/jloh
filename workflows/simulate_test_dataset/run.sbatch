#!/usr/bin/env bash 
#SBATCH -J sim_data 
#SBATCH -N 1 
#SBATCH -c 48 
#SBATCH --time 03:00:00
#SBATCH --error=/gpfs/projects/bsc40/current/mschiavi/jloh/scripts/wf-test_data_set/run.sbatch.stderr
#SBATCH --output=/gpfs/projects/bsc40/current/mschiavi/jloh/scripts/wf-test_data_set/run.sbatch.stdout

# import modules in your HPC cluster
# or simply have these executables in your $PATH 
# IMPORTANT: make sure that the executables declared in "nextflow.config" match the ones in your $PATH
module load gcc/7.2.0 hisat2/2.2.1 samtools/1.15 bcftools/1.15.1 htslib/1.15.1 bedtools/2.29.2

# change directory to where the workflow is 
WORKFLOW_DIR=""
cd $WORKFLOW_DIR

# run workflow 
nextflow \
-C nextflow.config \
run \
simulate_test_dataset.nf \
-resume \
-work-dir work 
