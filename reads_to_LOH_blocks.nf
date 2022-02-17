#!/usr/bin/env nextflow


if (params.help) {
println """

  // ----------------------------------------------------------------------//
  // Workflow to extract LOH blocks leveraging JLOH and other common tools //
  // ----------------------------------------------------------------------//

  // refer to:
  // https://github.com/Gabaldonlab/jloh
  // Author: Matteo Schiavinato
  // DOI: (yet to be published)

  [input]

--help              Print this help section                                     [off]
--run_id            ID of the processed sample (or the run, when >1 sample)     [smpl]
--output_dir        Name of output directory                                    [JLOH_run]
--threads           Number of parallel threads                                  [4]
--read_type         Select paired-end (PE) or single-end (SE)                   [PE]
--reads_for         If paired-end reads, forward reads                          [!]
--reads_rev         If paired-end reads, reverse reads                          [!]
--reads             If single-end reads, the only reads file                    [!]
--reads_dir         Directory containing *_1.fastq and *_2.fastq read files
                    (if paired) or all read files (if single)                   [off]
--ref_genome        Path to the reference genome in FASTA format                [!]
--ref_genome_dir    Directory containing all ref genomes to map the reads on
                    (must end in *.fa or in *.fasta)                            [off]

  [trimming]

--skip_trimming     Specify this option if you want to skip the trimming        [off]
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

  [JLOH]

--min_het_snps      Min. number of het SNPs to consider a het bloock            [2]
--snp_distance      Size (bp) of the window to use for LOH detection            [100]
--min_loh_size      Min. size (bp) of the candidate LOH blocks                  [100]
--block_distance    Combine LOH blocks into one if closer than this distance    [100]
--min_het_af        Min. allele frequency to consider heterozygous              [0.30]
--max_het_af        Max. allele frequency to consider heterozygous              [0.70]
--min_frac_cov      Min. fraction of LOH block that has to be covered by reads  [0.5]
--hemizygous_cov    Frac. of the mean coverage under which LOH is hemizygous    [0.75]
--alpha             Two-sided max. diff. from chromosomal cluster mean cov      [0.05]

"""
exit 0
}


// -----------------------------------------------------------------------------
// PAIRED END READS

// reading genomes

if ((!params.ref_genome_dir) && (params.ref_genome)) {

  Channel
    .fromPath("${params.ref_genome}")
    .map{ it ->  [it.baseName, it ]}
    .set{ Ref_genomes }

} else if ((params.ref_genome_dir) && (!params.ref_genome)) {

  Channel
    .fromPath("${params.ref_genome_dir}/*.{fa,fasta}")
    .map{ it ->  [it.baseName, it ]}
    .set{ Ref_genomes }

} else {

  println "ERROR: use either --ref_genome or --ref_genome_dir, but not both"
  exit 0

}

if (params.read_type == "PE") {

  // read input read files

  if ((!params.reads_dir) && (params.reads_for && params.reads_rev)) {

    Channel
      .fromPath("${params.reads_for}")
      .map{ it -> [ it.baseName, it ] }
      .set{ Reads_1 }

    Channel
      .fromPath("${params.reads_rev}")
      .map{ it -> [ it.baseName, it ] }
      .set{ Reads_2 }

    Reads_1
      .join(Reads_2)
      .into{ Reads_PE_QC; Reads_PE_TRIM; Reads_PE_RAW }

  } else if ((params.reads_dir) && !(params.reads_for && params.reads_rev)) {

    Channel
      .fromFilePairs("${params.reads_dir}/*_{1,2}.{fq,fastq}", flat: true)
      .into{ Reads_PE_QC; Reads_PE_TRIM; Reads_PE_RAW }

  } else {

    println "ERROR: use either --reads_for/--reads_rev or --reads_dir, but not both"
    exit 0

  }

  // trim reads if needed
  // then map them

  if (!params.skip_trimming) {

    process PE_quality_check_before_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check/before", mode: "copy"

      input:
        tuple val(sample_id), file(reads_for), file(reads_rev) from Reads_PE_QC

      output:
        path "${params.run_id}.${sample_id}"

      script:
        """
        if [ ! -d ${params.run_id}.${sample_id} ]; then mkdir ${params.run_id}.${sample_id}; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir ${params.run_id}.${sample_id} \
        ${reads_for} ${reads_rev}
        """
    }

    process PE_trim_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir  "${params.output_dir}/trimmed_reads",
                  mode: "copy", pattern: "*.{P1,P2,U1,U2}.fastq"

      input:
        tuple val(sample_id), file(reads_for), file(reads_rev) from Reads_PE_TRIM

      output:
        tuple val(sample_id), \
        file("${params.run_id}.${sample_id}.P1.fastq"), \
        file("${params.run_id}.${sample_id}.P2.fastq"), \
        file("${params.run_id}.${sample_id}.U.fastq") into Reads_PE_trimmed

      script:
        """
        ${TRIMMOMATIC} \
        PE \
        -threads ${params.threads} \
        -summary ${sample_id}.summary.tsv \
        ${reads_for} ${reads_rev} \
        ${params.run_id}.${sample_id}.P1.fastq ${params.run_id}.${sample_id}.U1.fastq \
        ${params.run_id}.${sample_id}.P2.fastq ${params.run_id}.${sample_id}.U2.fastq \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:${params.leading} \
        TRAILING:${params.trailing} \
        SLIDINGWINDOW:${params.sliding_window} \
        AVGQUAL:${params.avgqual} \
        MINLEN:${params.minlen} &&
        cat ${params.run_id}.${sample_id}.U1.fastq ${params.run_id}.${sample_id}.U2.fastq \
        > ${params.run_id}.${sample_id}.U.fastq
        """
    }

    Reads_PE_trimmed.into{ Reads_PE_trimmed_MAP; Reads_PE_trimmed_QC }

    process PE_quality_check_after_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check/after", mode: "copy"

      input:
        tuple val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired) from Reads_PE_trimmed_QC

      output:
        path "${params.run_id}.${sample_id}"

      script:
        """
        if [ ! -d ${params.run_id}.${sample_id} ]; then mkdir ${params.run_id}.${sample_id}; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir ${params.run_id}.${sample_id} \
        ${reads_for} ${reads_rev} ${reads_unpaired}
        """
    }

    // cartesian product of the reads x genomes channels
    Reads_PE_trimmed_MAP
      .combine(Ref_genomes)
      .set{ Hisat2_in }

    // output:
    // tuple val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired), \
    //       val(genome_id), file(ref_genome)

    process PE_map_trimmed_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      // publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

      input:
        tuple val(sample_id), file(reads_for), file(reads_rev), file(reads_unpaired), \
              val(genome_id), file(ref_genome) from Hisat2_in

      output:
        tuple val(sample_id), val(genome_id), file(ref_genome), \
              file("${params.run_id}.${sample_id}.${genome_id}.sam") into Hisat2_out

      script:
        """
        ${HISAT2_BUILD} -p ${params.threads} ${ref_genome} ${genome_id}.idx &&
        ${HISAT2} -p ${params.threads} \
        --no-spliced-alignment \
        --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
        --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
        -I ${params.min_isize} -X ${params.max_isize} \
        -x ${genome_id}.idx \
        -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
        -S ${params.run_id}.${sample_id}.${genome_id}.sam
        """
    }
  } else {

    // cartesian product of the reads x genomes channels
    Reads_PE_RAW
      .combine(Ref_genomes)
      .set{ Hisat2_in }

    // output:
    // tuple val(sample_id), file(reads_for), file(reads_rev), \
    //       val(genome_id), file(ref_genome)

    process PE_map_untrimmed_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      // publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

      input:
        tuple val(sample_id), file(reads_for), file(reads_rev), \
              val(genome_id), file(ref_genome) from Hisat2_in

      output:
        tuple val(sample_id), val(genome_id), file(ref_genome), \
              file("${params.run_id}.${sample_id}.${genome_id}.sam") into Hisat2_out

      script:
        """
        ${HISAT2} -p ${params.threads} \
        --no-spliced-alignment \
        --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
        --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
        -I ${params.min_isize} -X ${params.max_isize} \
        -x ${genome_id}.idx \
        -1 ${reads_for} -2 ${reads_rev} \
        -S ${params.run_id}.${sample_id}.${genome_id}.sam
        """
    }
  }


// -----------------------------------------------------------------------------
// SINGLE END READS

} else if (params.read_type == "SE") {

  if ((!params.reads_dir) && (params.reads)) {

    Channel
      .fromPath("${params.reads}")
      .map{ it -> [it.baseName, it]}
      .into{ Reads_SE_QC; Reads_SE_TRIM; Reads_SE_RAW }

  } else if ((params.reads_dir) && (!params.reads)) {

    Channel
      .fromPath("${params.reads}*.{fq,fastq}")
      .map{ it -> [it.baseName, it]}
      .into{ Reads_SE_QC; Reads_SE_TRIM; Reads_SE_RAW }

  } else {

    println "ERROR: use either --reads or --reads_dir, but not both"
    exit 0

  }


  if (!params.skip_trimming) {

    process SE_quality_check_before_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check/before", mode: "copy"

      input:
        tuple val(sample_id), file(reads) from Reads_SE_QC

      output:
        path "${params.run_id}.${sample_id}"

      script:
        """
        if [ ! -d ${params.run_id}.${sample_id} ]; then mkdir ${params.run_id}.${sample_id}; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir ${params.run_id}.${sample_id} \
        ${reads}
        """
    }

    process SE_trim_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir  "${params.output_dir}/trimmed_reads",
                  mode: "copy", pattern: "*.trimmed.fastq"

      input:
        tuple val(sample_id), file(reads) from Reads_SE_TRIM

      output:
        tuple val(sample_id), file("${params.run_id}.${sample_id}.trimmed.fastq") into Reads_SE_trimmed

      script:
        """
        ${TRIMMOMATIC} \
        SE \
        -threads ${params.threads} \
        -summary ${sample_id}.summary.tsv \
        ${reads} ${params.run_id}.${sample_id}.trimmed.fastq \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:${params.leading} \
        TRAILING:${params.trailing} \
        SLIDINGWINDOW:${params.sliding_window} \
        AVGQUAL:${params.avgqual} \
        MINLEN:${params.minlen}
        """
    }

    Reads_SE_trimmed.into{ Reads_SE_trimmed_MAP; Reads_SE_trimmed_QC }

    process SE_quality_check_after_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check/after", mode: "copy"

      input:
        tuple val(sample_id), file(reads) from Reads_SE_trimmed_QC

      output:
        path "${params.run_id}.${sample_id}"

      script:
        """
        if [ ! -d ${params.run_id}.${sample_id} ]; then mkdir ${params.run_id}.${sample_id}; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir ${params.run_id}.${sample_id} \
        ${reads}
        """
    }

    // cartesian product of the reads x genomes channels
    Reads_SE_trimmed_MAP
      .combine(Ref_genomes)
      .set{ Hisat2_in }

    // output:
    // tuple val(sample_id), file(reads), val(genome_id), file(ref_genome)

    process SE_map_trimmed_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      // publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

      input:
        tuple val(sample_id), file(reads), val(genome_id), file(ref_genome) from Hisat2_in

      output:
        tuple val(sample_id), val(genome_id), file(ref_genome), \
              file("${params.run_id}.${sample_id}.${genome_id}.sam") into Hisat2_out

      script:
        """
        ${HISAT2_BUILD} -p ${params.threads} ${ref_genome} ${genome_id}.idx &&
        ${HISAT2} -p ${params.threads} \
        --no-spliced-alignment \
        --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
        --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
        -I ${params.min_isize} -X ${params.max_isize} \
        -x ${genome_id}.idx \
        -U ${reads} \
        -S ${params.run_id}.${sample_id}.${genome_id}.sam
        """
      }

  } else {

    // cartesian product of the reads x genomes channels
    Reads_SE_RAW
      .combine(Ref_genomes)
      .set{ Hisat2_in }

    // output:
    // tuple val(sample_id), file(reads), val(genome_id), file(ref_genome)

    process SE_map_untrimmed_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      // publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

      input:
        tuple val(sample_id), file(reads), val(genome_id), file(ref_genome) from Hisat2_in

      output:
        tuple val(sample_id), val(genome_id), file(ref_genome), \
              file("${params.run_id}.${sample_id}.${genome_id}.sam") into Hisat2_out

      script:
        """
        ${HISAT2_BUILD} -p ${params.threads} ${ref_genome} ${genome_id}.idx &&
        ${HISAT2} -p ${params.threads} \
        --no-spliced-alignment \
        --score-min ${params.scoring_fun} --mp ${params.mm_penalties} \
        --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
        -I ${params.min_isize} -X ${params.max_isize} \
        -x ${genome_id}.idx \
        -U ${reads} \
        -S ${params.run_id}.${sample_id}.${genome_id}.sam
        """
    }
  }
}

// ---------------------------------------------------------------------------
// COMMON STEPS

process filter_and_sort_bam {

  executor = "local"
  cpus = params.threads
  maxForks = params.threads

  publishDir  "${params.output_dir}/mapping",
              mode: "copy", pattern: "*.{fs.bam,fs.bam.bai}"

  input:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(sam) from Hisat2_out

  output:
    tuple val(sample_id), val(genome_id), file(ref_genome), \
          file("${params.run_id}.${sample_id}.${genome_id}.fs.bam"), \
          file("${params.run_id}.${sample_id}.${genome_id}.fs.bam.bai") \
          into Filt_Sort_Bam

  script:
    """
    ${SAMTOOLS} view --threads ${params.threads} -F 0x0100 -F 0x4 -h -b \
    -o ${params.run_id}.${sample_id}.${genome_id}.f.bam ${sam} &&
    ${SAMTOOLS} sort --threads ${params.threads} --output-fmt bam \
    -T ${params.run_id}.${sample_id}.${genome_id} \
    -o ${params.run_id}.${sample_id}.${genome_id}.fs.bam \
    ${params.run_id}.${sample_id}.${genome_id}.f.bam &&
    ${SAMTOOLS} index -@ ${params.threads} -b \
    ${params.run_id}.${sample_id}.${genome_id}.fs.bam
    """
}


process pileup_reads {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  // publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.mpileup.vcf"

  input:
    tuple val(sample_id), val(genome_id), file(ref_genome), \
          file(bam_fs), file(bam_fs_bai) from Filt_Sort_Bam

  output:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), file(bam_fs_bai), \
          file("${params.run_id}.${sample_id}.${genome_id}.mpileup.vcf") into Pileups

  script:
    """
    ${SAMTOOLS} faidx ${ref_genome} &&
    ${BCFTOOLS} mpileup --fasta-ref ${ref_genome} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${params.run_id}.${sample_id}.${genome_id}.mpileup.vcf --output-type v --skip-indels \
    ${bam_fs}
    """
}

process call_short_variants {

  executor = "local"
  cpus = params.threads
  maxForks = params.threads

  // publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.raw.vcf"

  input:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), \
          file(bam_fs_bai), file(pileup) from Pileups

  output:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), file(bam_fs_bai), \
          file("${params.run_id}.${sample_id}.${genome_id}.raw.vcf") into Raw_vcfs

  script:
    """
    ${BCFTOOLS} call \
    --threads ${params.threads} \
    --multiallelic-caller --variants-only \
    --output ${params.run_id}.${sample_id}.${genome_id}.raw.vcf --output-type v \
    ${pileup}
    """
}

process filter_short_variants {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.{ff.vcf,report.txt}"

  input:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), file(bam_fs_bai), \
          file(raw_vcf) from Raw_vcfs

  output:
    file "${params.run_id}.${sample_id}.${genome_id}.report.txt" into Reports

    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), file(bam_fs_bai), \
          file("${params.run_id}.${sample_id}.${genome_id}.ff.vcf") into Filt_vcfs

  script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file ${raw_vcf} \
    --output-file ${params.run_id}.${sample_id}.${genome_id}.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_alt_frac} \
    --min-depth ${params.min_depth} \
    --map-qual-zero-frac ${params.mq0f} \
    --threads 1 \
    --report ${params.run_id}.${sample_id}.${genome_id}.report.txt &&
    ${ALL2VCF} frequency \
    --in ${params.run_id}.${sample_id}.${genome_id}.f.vcf \
    --out ${params.run_id}.${sample_id}.${genome_id}.ff.vcf

    """
}

process call_LOH_blocks {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/LOH", mode: "copy"

  input:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), file(bam_fs_bai), \
          file(filt_vcf) from Filt_vcfs

  output:
    tuple val(sample_id), val(genome_id), file(ref_genome), file(bam_fs), file(bam_fs_bai), \
          path("${params.run_id}.${sample_id}.${genome_id}") into Jloh_out

  script:
    """
    ${BIOAWK} -c fastx '{print \$name\"\\t\"length(\$seq)}' ${ref_genome} \
    > ${params.run_id}.genome_file.tsv &&
    ${JLOH} \
    --threads ${params.threads} \
    --vcf ${filt_vcf} \
    --bam ${bam_fs} \
    --genome-file ${params.run_id}.genome_file.tsv \
    --sample ${params.run_id}.${sample_id}.${genome_id} \
    --output-dir ${params.run_id}.${sample_id}.${genome_id} \
    --filter-mode all \
    --min-af ${params.min_het_af} \
    --max-af ${params.max_het_af} \
    --min-frac-cov ${params.min_frac_cov} \
    --min-het-snps ${params.min_het_snps} \
    --snp-distance ${params.snp_distance} \
    --block-dist ${params.block_distance} \
    --min-size ${params.min_loh_size} \
    --hemi ${params.hemizygous_cov} \
    --alpha ${params.alpha}

    """
}
