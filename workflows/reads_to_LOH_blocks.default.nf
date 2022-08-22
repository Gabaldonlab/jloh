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

--skip_trimming     Skip this step                                              [off]
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
// input reads

// sample_id, read_1, read_2, genome_seq

Channel
  .fromPath("${params.input_data}")
  .splitCsv(sep: "\t", header: false)
  .map{ it -> [it[0], file(it[1]), file(it[2]), file(it[3])] }
  .into{ Input_pretrim; Input_trim; Input_notrim }

// -----------------------------------------------------------------------------
// reads trimming

if (!params.skip_trimming) {

  process PE_quality_check_before_trim {

    executor = "local"
    cpus = 4
    maxForks = params.threads

    publishDir "${params.output_dir}/trimmed_reads/quality_check/before", mode: "copy"

    input:
    tuple \
    val(sample_id), file(reads_for), file(reads_rev), file(ref) \
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
    tuple \
    val(sample_id), file(reads_for), file(reads_rev), file(ref) \
    from Input_trim

    output:
    tuple \
    val(sample_id), file("${sample_id}.P1.fastq"), file("${sample_id}.P2.fastq"), \
    file("${sample_id}.U.fastq"), file(ref) \
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
    tuple \
    val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired), file(ref) \
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
    tuple \
    val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired), file(ref) \
    from Reads_PE_trimmed_MAP

    output:
    tuple \
    val(sample_id), file("${sample_id}.sam"), file(ref) \
    into Hisat2_out

    script:
    """
    ${HISAT2_BUILD} -p ${params.threads} ${ref} ref.idx &&
    ${HISAT2} -p ${params.threads} \
    --no-spliced-alignment \
    --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
    --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
    -I ${params.min_isize} -X ${params.max_isize} \
    -x ref.idx \
    -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
    -S ${sample_id}.sam
    """
  }
} else {

  process PE_map_reads {

    executor = "local"
    cpus = params.threads
    maxForks = 1

    input:
    tuple \
    val(sample_id), file(reads_for), file(reads_rev), file(ref) \
    from Input_notrim

    output:
    tuple \
    val(sample_id), file("${sample_id}.sam"), file(ref) \
    into Hisat2_out

    script:
    """
    ${HISAT2_BUILD} -p ${params.threads} ${ref} ref.idx &&
    ${HISAT2} -p ${params.threads} \
    --no-spliced-alignment \
    --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
    --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
    -I ${params.min_isize} -X ${params.max_isize} \
    -x ref.idx \
    -1 ${reads_for} -2 ${reads_rev} \
    -S ${sample_id}.sam
    """
  }
}


process filter_and_sort_bam {

  executor = "local"
  cpus = 8
  maxForks = 6

  publishDir  "${params.output_dir}/mapping",
              mode: "copy", pattern: "*.{fs.bam,fs.bam.bai}"

  input:
    tuple \
    val(sample_id), file(sam), file(ref) \
    from Hisat2_out

  output:
    tuple \
    val(sample_id), file("${sample_id}.fs.bam"), file("${sample_id}.fs.bam.bai"), file(ref) \
    into Filt_Sort_Bam

  script:
    """
    ${SAMTOOLS} view --threads ${params.threads} -F 0x0100 -F 0x4 -h -b \
    -o ${sample_id}.f.bam ${sam} &&
    ${SAMTOOLS} sort --threads 8 --output-fmt bam \
    -T ${sample_id} \
    -o ${sample_id}.fs.bam \
    ${sample_id}.f.bam &&
    ${SAMTOOLS} index -@ ${params.threads} -b \
    ${sample_id}.fs.bam
    """
}


process pileup_reads {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  input:
    tuple \
    val(sample_id), file(bam), file(bai), file(ref) \
    from Filt_Sort_Bam

  output:
    tuple \
    val(sample_id), file(bam), file(bai), file("${sample_id}.mpileup.vcf"), file(ref) \
    into Pileups

  script:
    """
    ${SAMTOOLS} faidx ${ref} &&
    ${BCFTOOLS} mpileup --fasta-ref ${ref} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${sample_id}.mpileup.vcf --output-type v --skip-indels \
    ${bam}
    """
}

process call_short_variants {

  executor = "local"
  cpus = 8
  maxForks = params.threads

  input:
    tuple \
    val(sample_id), file(bam), file(bai), file(pileup), file(ref) \
    from Pileups

  output:
    tuple \
    val(sample_id), file(bam), file(bai), file("${sample_id}.raw.vcf"), file(ref) \
    into Raw_vcfs

  script:
    """
    ${BCFTOOLS} call \
    --threads ${params.threads} \
    --multiallelic-caller --variants-only \
    --output ${sample_id}.raw.vcf --output-type v \
    ${pileup}
    """
}

process filter_short_variants {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.{ff.vcf,report.txt}"

  input:
    tuple \
    val(sample_id), file(bam), file(bai), file(raw_vcf), file(ref) \
    from Raw_vcfs

  output:
    file "${sample_id}.report.txt" into Reports

    tuple \
    val(sample_id), file(bam), file(bai), file("${sample_id}.ff.vcf"), file(ref) \
    into Filt_vcfs

  script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file ${raw_vcf} \
    --output-file ${sample_id}.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_alt_frac} \
    --min-depth ${params.min_depth} \
    --map-qual-zero-frac ${params.mq0f} \
    &> ${sample_id}.report.txt &&
    ${ALL2VCF} frequency \
    --in ${sample_id}.f.vcf \
    --out ${sample_id}.ff.vcf
    """
}


process call_LOH_blocks {

  executor = "local"
  cpus = 16
  maxForks = 3

  publishDir "${params.output_dir}/LOH", mode: "copy", pattern: "*/*.{tsv,bed}"

  input:
    tuple \
    val(sample_id), file(bam), file(bai), file(vcf), file(ref) \
    from Filt_vcfs

  output:
    tuple \
    val(sample_id), \
    file("${sample_id}/${sample_id}.LOH_blocks.tsv"), \
    file("${sample_id}/${sample_id}.LOH_blocks.bed"), \
    file("${sample_id}/${sample_id}.exp.het_blocks.bed"), \
    file("${sample_id}/${sample_id}.exp.chrom_coverage.tsv") \
    into Jloh_out

  script:
    """
    ${JLOH} extract \
    --threads ${params.threads} \
    --vcf ${vcf} \
    --bam ${bam} \
    --ref ${ref} \
    --sample ${sample_id} \
    --output-dir ${sample_id} \
    --filter-mode ${params.var_filter} \
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
