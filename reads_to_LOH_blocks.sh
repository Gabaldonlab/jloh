#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -J JLOH
#SBATCH --time=48:00:00
#SBATCH --error=reads_to_LOH_blocks.run.stderr
#SBATCH --output=reads_to_LOH_blocks.run.stdout

SAMPLE_ID="put_a_SAMPLE_ID_for_your_output_files"

nextflow \
-C reads_to_LOH_blocks.config \
run \
reads_to_LOH_blocks.nf \
-resume \
-work-dir ${SAMPLE_ID}.work \
-with-dag ${SAMPLE_ID}.svg \
-with-report ${SAMPLE_ID}.report \
-with-timeline ${SAMPLE_ID}.timeline \
--read_type PE \
--reads_for <path_to_reads_forward> \
--reads_rev <path_to_reads_reverse> \
--ref_genome <path_to_ref_genome> \
--threads 48 \
--run_id ${SAMPLE_ID} \
--output_dir ${SAMPLE_ID}
