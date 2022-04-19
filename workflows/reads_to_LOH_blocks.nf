#!/usr/bin/env nextflow


if (params.help) {
println """

  // ------------------------------------------------------------------//
  // Workflow to extract LOH blocks from public paired-end hybrid data //
  // ------------------------------------------------------------------//

  [input]

--help              Print this help section                                     [off]
--output_dir        Name of output directory                                    [JLOH_run]
--threads           Number of parallel threads                                  [4]
--input_data        Table containing input data paths and info (see github)     [!]
--reads_dir         Directory containing the reads (for double-checking)        [!]

  [trimming]

--adapters          Path to the adapters file to pass to trimmomatic            [off]
--leading           Trimmomatic "LEADING" setting                               [20]
--trailing          Trimmomatic "TRAILING" setting                              [20]
--sliding_window    Trimmomatic "SLIDINGWINDOW" setting                         [4:25]
--avgqual           Trimmomatic "AVGQUAL" setting                               [20]
--minlen            Trimmomatic "MINLEN" setting                                [35]

  [mapping]

--scoring_fun       hisat2 scoring function                                     [L,0.0,-1.0]
--mm_penalties      hisat2 mismatch penalties                                   [6,2]
--gap_penalties     hisat2 gap penalties                                        [5,3]
--min_isize         hisat2 minimum insert size                                  [0]
--max_isize         hisat2 maximum insert size                                  [1000]

  [variant filtering]

--min_qual          Minimum variant quality (QUAL)                              [20]
--min_alt_frac      Minimum variant alternative allele fraction (AF)            [0.05]
--min_depth         Minimum variant coverage depth (DP)                         [10]
--mq0f              Maximum "mapping-quality-0" fraction of reads               [0.05]

  [genome divergence]

--min_mum_length    Minimum length to retain a genome alignment                 [1000]
--min_mum_id        Minimum sequence identity to retain a genome mapping        [95]

  [JLOH]

--min_snps          Min. number of proximal het SNPs to consider a het block    [2]
--snp_distance      Max. distance between SNPs to assign to same block          [100]
--min_loh_size      Min. size (bp) of the candidate LOH blocks                  [1000]
--min_af            Min. allele frequency to consider heterozygous              [0.2]
--max_af            Max. allele frequency to consider heterozygous              [0.8]
--min_frac_cov      Min. fraction of LOH block that has to be covered by reads  [0.5]
--hemizygous_cov    Frac. of the mean coverage under which LOH is hemizygous    [0.75]

"""
exit 0
}


// -----------------------------------------------------------------------------
// PAIRED END READS

// read input data

def lines = new File("${params.input_data}").collect {it}

Channel
  .fromList(lines)
  .map{ it -> it.tokenize("\t") }
  .map{ it -> [ it[0], file(it[1]), file(it[2]), it[3], file(it[4]), it[5], file(it[6]) ] }
  .set{ Input_data }

Channel
  .fromFilePairs("${params.reads_dir}/*_{1,2}.{fq,fastq}", flat: true)
  .set{ Reads }

Input_data
  .join(Reads)
  .map{ it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6] ] }
  .into{ Input_pretrim; Input_trim }


// -----------------------------------------------------------------------------
// reads trimming

process PE_quality_check_before_trim {

  executor = "local"
  cpus = 4
  maxForks = params.threads

  publishDir "${params.output_dir}/trimmed_reads/quality_check/before", mode: "copy"

  input:
    tuple val(sample_id), file(reads_for), file(reads_rev), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Input_pretrim

  output:
    path "${sample_id}"

  script:
    """
    if [ ! -d ${sample_id} ]; then mkdir ${sample_id}; fi &&
    ${FASTQC} \
    --threads ${params.threads} \
    --outdir ${sample_id} \
    ${reads_for} ${reads_rev}
    """
}

process PE_trim_reads {

  executor = "local"
  cpus = 8
  maxForks = params.threads

  publishDir  "${params.output_dir}/trimmed_reads",
              mode: "copy", pattern: "*.{P1,P2,U1,U2}.fastq"

  input:
    tuple val(sample_id), file(reads_for), file(reads_rev), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Input_trim

  output:
    tuple val(sample_id), \
    file("${sample_id}.P1.fastq"), \
    file("${sample_id}.P2.fastq"), \
    file("${sample_id}.U.fastq"), \
    val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
    into Reads_PE_trimmed

  script:
    """
    ${TRIMMOMATIC} \
    PE \
    -threads ${params.threads} \
    -summary ${sample_id}.summary.tsv \
    ${reads_for} ${reads_rev} \
    ${sample_id}.P1.fastq ${sample_id}.U1.fastq \
    ${sample_id}.P2.fastq ${sample_id}.U2.fastq \
    ILLUMINACLIP:${params.adapters}:2:30:10 \
    LEADING:${params.leading} \
    TRAILING:${params.trailing} \
    SLIDINGWINDOW:${params.sliding_window} \
    AVGQUAL:${params.avgqual} \
    MINLEN:${params.minlen} &&
    cat ${sample_id}.U1.fastq ${sample_id}.U2.fastq \
    > ${sample_id}.U.fastq
    """
}

Reads_PE_trimmed.into{ Reads_PE_trimmed_MAP; Reads_PE_trimmed_QC }

process PE_quality_check_after_trim {

  executor = "local"
  cpus = 4
  maxForks = params.threads

  publishDir "${params.output_dir}/trimmed_reads/quality_check/after", mode: "copy"

  input:
    tuple val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Reads_PE_trimmed_QC

  output:
    path "${sample_id}"

  script:
    """
    if [ ! -d ${sample_id} ]; then mkdir ${sample_id}; fi &&
    ${FASTQC} \
    --threads ${params.threads} \
    --outdir ${sample_id} \
    ${reads_for} ${reads_rev} ${reads_unpaired}
    """
}


process PE_map_trimmed_reads {

  executor = "local"
  cpus = params.threads
  maxForks = 1

  input:
    tuple val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Reads_PE_trimmed_MAP

  output:
    tuple val(sample_id), \
          file("${sample_id}.${genome_A_name}.sam"), file("${sample_id}.${genome_B_name}.sam"), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          into Hisat2_out

  script:
    """
    ${HISAT2_BUILD} -p ${params.threads} ${genome_A} ${genome_A_name}.idx &&
    ${HISAT2_BUILD} -p ${params.threads} ${genome_B} ${genome_B_name}.idx &&
    ${HISAT2} -p ${params.threads} \
    --no-spliced-alignment \
    --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
    --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
    -I ${params.min_isize} -X ${params.max_isize} \
    -x ${genome_A_name}.idx \
    -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
    -S ${sample_id}.${genome_A_name}.sam &&
    ${HISAT2} -p ${params.threads} \
    --no-spliced-alignment \
    --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
    --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
    -I ${params.min_isize} -X ${params.max_isize} \
    -x ${genome_B_name}.idx \
    -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
    -S ${sample_id}.${genome_B_name}.sam
    """
}

process filter_and_sort_bam {

  executor = "local"
  cpus = 4
  maxForks = params.threads

  publishDir  "${params.output_dir}/mapping",
              mode: "copy", pattern: "*.{fs.bam,fs.bam.bai}"

  input:
    tuple val(sample_id), \
          file(par_A_sam), file(par_B_sam), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Hisat2_out

  output:
    tuple val(sample_id), \
          file("${sample_id}.${genome_A_name}.fs.bam"), file("${sample_id}.${genome_A_name}.fs.bam.bai"), \
          file("${sample_id}.${genome_B_name}.fs.bam"), file("${sample_id}.${genome_B_name}.fs.bam.bai"), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          into Filt_Sort_Bam

  script:
    """
    ${SAMTOOLS} view --threads ${params.threads} -F 0x0100 -F 0x4 -h -b \
    -o ${sample_id}.${genome_A_name}.f.bam ${par_A_sam} &&
    ${SAMTOOLS} view --threads ${params.threads} -F 0x0100 -F 0x4 -h -b \
    -o ${sample_id}.${genome_B_name}.f.bam ${par_B_sam} &&
    ${SAMTOOLS} sort --threads ${params.threads} --output-fmt bam \
    -T ${sample_id}.${genome_A_name} \
    -o ${sample_id}.${genome_A_name}.fs.bam \
    ${sample_id}.${genome_A_name}.f.bam &&
    ${SAMTOOLS} sort --threads ${params.threads} --output-fmt bam \
    -T ${sample_id}.${genome_B_name} \
    -o ${sample_id}.${genome_B_name}.fs.bam \
    ${sample_id}.${genome_B_name}.f.bam &&
    ${SAMTOOLS} index -@ ${params.threads} -b \
    ${sample_id}.${genome_A_name}.fs.bam &&
    ${SAMTOOLS} index -@ ${params.threads} -b \
    ${sample_id}.${genome_B_name}.fs.bam
    """
}


process pileup_reads {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  input:
    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Filt_Sort_Bam

  output:
    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), file("${sample_id}.${genome_A_name}.mpileup.vcf"), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), file("${sample_id}.${genome_B_name}.mpileup.vcf"), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          into Pileups

  script:
    """
    ${SAMTOOLS} faidx ${genome_A} &&
    ${SAMTOOLS} faidx ${genome_B} &&
    ${BCFTOOLS} mpileup --fasta-ref ${genome_A} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${sample_id}.${genome_A_name}.mpileup.vcf --output-type v --skip-indels \
    ${par_A_fs_bam} &&
    ${BCFTOOLS} mpileup --fasta-ref ${genome_B} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${sample_id}.${genome_B_name}.mpileup.vcf --output-type v --skip-indels \
    ${par_B_fs_bam}
    """
}

process call_short_variants {

  executor = "local"
  cpus = 8
  maxForks = params.threads

  input:
    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), file(par_A_pileup), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), file(par_B_pileup), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Pileups

  output:
    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), file("${sample_id}.${genome_A_name}.raw.vcf"), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), file("${sample_id}.${genome_B_name}.raw.vcf"), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          into Raw_vcfs

  script:
    """
    ${BCFTOOLS} call \
    --threads ${params.threads} \
    --multiallelic-caller --variants-only \
    --output ${sample_id}.${genome_A_name}.raw.vcf --output-type v \
    ${par_A_pileup} &&
    ${BCFTOOLS} call \
    --threads ${params.threads} \
    --multiallelic-caller --variants-only \
    --output ${sample_id}.${genome_B_name}.raw.vcf --output-type v \
    ${par_B_pileup}
    """
}

process filter_short_variants {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.{ff.vcf,report.txt}"

  input:
    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), file(par_A_raw_vcf), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), file(par_B_raw_vcf), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          from Raw_vcfs

  output:
    tuple file("${sample_id}.${genome_A_name}.report.txt"), \
          file("${sample_id}.${genome_B_name}.report.txt") \
          into Reports

    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), file("${sample_id}.${genome_A_name}.ff.vcf"), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), file("${sample_id}.${genome_B_name}.ff.vcf"), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B) \
          into Filt_vcfs

  script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file ${par_A_raw_vcf} \
    --output-file ${sample_id}.${genome_A_name}.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_alt_frac} \
    --min-depth ${params.min_depth} \
    --map-qual-zero-frac ${params.mq0f} \
    &> ${sample_id}.${genome_A_name}.report.txt &&
    ${ALL2VCF} frequency \
    --in ${sample_id}.${genome_A_name}.f.vcf \
    --out ${sample_id}.${genome_A_name}.ff.vcf &&
    ${ALL2VCF} filter_vcf \
    --input-file ${par_B_raw_vcf} \
    --output-file ${sample_id}.${genome_B_name}.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_alt_frac} \
    --min-depth ${params.min_depth} \
    --map-qual-zero-frac ${params.mq0f} \
    &> ${sample_id}.${genome_B_name}.report.txt &&
    ${ALL2VCF} frequency \
    --in ${sample_id}.${genome_B_name}.f.vcf \
    --out ${sample_id}.${genome_B_name}.ff.vcf
    """
}


process call_LOH_blocks {

  executor = "local"
  cpus = 16
  maxForks = 3

  publishDir "${params.output_dir}/LOH", mode: "copy", pattern: "*/*.{tsv,vcf}"

  input:
    tuple val(sample_id), \
          file(par_A_fs_bam), file(par_A_fs_bam_index), file(par_A_vcf), \
          file(par_B_fs_bam), file(par_B_fs_bam_index), file(par_B_vcf), \
          val(genome_A_name), file(genome_A), val(genome_B_name), file(genome_B), \
          file(bed_mask) \
          from Filt_vcfs

  output:
    tuple val(sample_id), val(genome_A_name), val(genome_B_name), file(bed_mask), \
          file("${sample_id}/${sample_id}.exp_A.LOH_blocks.tsv"), \
          file("${sample_id}/${sample_id}.exp_B.LOH_blocks.tsv"), \
          file("${sample_id}/${sample_id}.exp_A.LOH_blocks.bed"), \
          file("${sample_id}/${sample_id}.exp_B.LOH_blocks.bed"), \
          file("${sample_id}/${sample_id}.exp_A.chrom_coverage.tsv"), \
          file("${sample_id}/${sample_id}.exp_B.chrom_coverage.tsv") \
          into Jloh_out

  script:
    """
    ${JLOH} extract \
    --threads ${params.threads} \
    --vcfs ${par_A_vcf} ${par_B_vcf} \
    --bams ${par_A_fs_bam} ${par_B_fs_bam} \
    --refs ${genome_A} ${genome_B} \
    --hybrid \
    --sample ${sample_id} \
    --output-dir ${sample_id} \
    --filter-mode all \
    --min-af ${params.min_af} \
    --max-af ${params.max_af} \
    --min-frac-cov ${params.min_frac_cov} \
    --min-snps ${params.min_snps} \
    --snp-distance ${params.snp_distance} \
    --min-length ${params.min_loh_size} \
    --hemi ${params.hemizygous_cov} \
    --min-uncov ${params.min_uncovered}
    """
}
