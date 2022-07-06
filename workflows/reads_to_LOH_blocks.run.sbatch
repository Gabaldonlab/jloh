#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -c 48
#SBATCH -J JLOH

### INSTRUCTIONS ###

# 1. copy this sbatch file, the example input data tsv, and the config file to your directories
# 2. edit the config file using the parameters you want to use (otherwise leave default)
# 3. prepare the input data table following the example in /workflows
#    1) sample_id 2) reads_for 3) reads_rev 4) parent_A_id 5) parent_A_fasta 6) parent_B_id 7) parent_B_fasta
# 4. create the output directory where you want these results to be put
# 5. edit the fields in the next section of this file
# 6. submit this sbatch file with "sbatch reads_to_LOH_blocks.run.sbatch"
# 7. if it crashes: submit it again (will restart from where it crashed)
# 8. if it shows errors on the code side: open an issue on the github page and i will fix it!
# 9. once completed:
#    - remove ${OUTDIR}/nf/${SAMPLE_ID}.work directory (tmp files)
#    - notify me for pipeline directory cleanup

### EDIT FROM HERE ###

SAMPLE_ID="<put here a sample ID (no whitespaces) to identify output files>"
OUTDIR="<put here an existing output directory where everything will be written>"
INPUT_DATA="<path to the table with input data, see example in /workflows>"
ADAPTERS="<path to the FASTA file with the adapter sequences to trim, see /workflows for a default>"
CONFIG_FILE="<path to the config file, find template in /workflows>"
THREADS="<number of parallel threads, max 48>"

### STOP EDITING FROM HERE ###

module purge
module load intel mkl/2017.4 impi/2017.4
module load gcc/7.2.0 hisat2/2.2.1 samtools/1.15 bcftools/1.15.1 htslib/1.15.1 bedtools/2.29.2
module load graphviz/2.40.1 java/8u131 fastqc/0.11.9 trimmomatic/0.39 bioawk/1.0

if [ ! -d ${OUTDIR} ]; then mkdir ${OUTDIR}; fi
if [ ! -d ${OUTDIR}/nf ]; then mkdir ${OUTDIR}/nf; fi

nextflow \
-C ${CONFIG_FILE} \
run \
/gpfs/projects/bsc40/project/pipelines/JLOH/workflows/reads_to_LOH_blocks.nf \
-resume \
-work-dir ${OUTDIR}/nf/${SAMPLE_ID}.work \
-with-dag ${OUTDIR}/nf/${SAMPLE_ID}.svg \
-with-report ${OUTDIR}/nf/${SAMPLE_ID}.report \
-with-timeline ${OUTDIR}/nf/${SAMPLE_ID}.timeline \
--input_data ${INPUT_DATA} \
--adapters ${ADAPTERS} \
--threads ${THREADS} \
--run_id ${SAMPLE_ID} \
--output_dir ${OUTDIR}
