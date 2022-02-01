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
--read_type         Select paired-end (PE) or single-end (SE)                   [PE]
--reads_for         If paired-end reads, forward reads                          [!]
--reads_rev         If paired-end reads, reverse reads                          [!]
--reads             If single-end reads, the only reads file                    [!]
--ref_genome        Path to the reference genome in FASTA format                [!]
--threads           Number of parallel threads                                  [4]
--output_dir        Name of output directory                                    [JLOH_run]
--sample_id         ID of the processed sample                                  [smpl]

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

  --min_het_snps    Min. number of het SNPs to consider a het bloock            [2]
  --snp_distance    Size (bp) of the window to use for LOH detection            [100]
  --min_loh_size    Min. size (bp) of the candidate LOH blocks                  [100]
  --block_distance  Combine LOH blocks into one if closer than this distance    [100]
  --min_het_af      Min. allele frequency to consider heterozygous              [0.30]
  --max_het_af      Max. allele frequency to consider heterozygous              [0.70]
  --min_frac_cov    Min. fraction of LOH block that has to be covered by reads  [0.5]
  --hemizygous_cov  Frac. of the mean coverage under which LOH is hemizygous    [0.75]
  --alpha           Two-sided max. diff. from chromosomal cluster mean cov      [0.05]

"""
exit 0
}


// -----------------------------------------------------------------------------
// PAIRED END READS

process index_ref_genome {

  executor = "local"
  cpus = params.threads
  maxForks = params.threads

  output:
    file "${params.sample_id}.genome_index" into genome_index

  script:
    """
    ${HISAT2_BUILD} -p ${params.threads} ${params.ref_genome} ${params.sample_id}.genome_index
    """
}

if (params.read_type == "PE") {

  if (!params.skip_trimming) {

    process PE_quality_check_before_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check", mode: "copy"

      output:
        path "before_trimming"

      script:
        """
        if [ ! -d before_trimming ]; then mkdir before_trimming; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir before_trimming \
        ${params.reads_for} ${params.reads_rev}
        """
    }

    process PE_trim_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir  "${params.output_dir}/trimmed_reads",
                  mode: "copy", pattern: "*.{P1,P2,U1,U2}.fastq"

      output:
        tuple file("${params.sample_id}.P1.fastq"), file("${params.sample_id}.P2.fastq"), \
        file("${params.sample_id}.U.fastq") into Trimmed_reads

      script:
        """
        ${TRIMMOMATIC} \
        PE \
        -threads ${params.threads} \
        -summary trim.summary.tsv \
        ${params.reads_for} ${params.reads_rev} \
        ${params.sample_id}.P1.fastq ${params.sample_id}.U1.fastq \
        ${params.sample_id}.P2.fastq ${params.sample_id}.U2.fastq \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:${params.leading} \
        TRAILING:${params.trailing} \
        SLIDINGWINDOW:${params.sliding_window} \
        AVGQUAL:${params.avgqual} \
        MINLEN:${params.minlen} &&
        cat ${params.sample_id}.U1.fastq ${params.sample_id}.U2.fastq \
        > ${params.sample_id}.U.fastq
        """
    }

    Trimmed_reads.into{ Trimmed_reads_MAP; Trimmed_reads_QC }

    process PE_quality_check_after_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check", mode: "copy"

      input:
        tuple file(reads_for), file(reads_rev), file(reads_unpaired) from Trimmed_reads_QC

      output:
        path "after_trimming"

      script:
        """
        if [ ! -d after_trimming ]; then mkdir after_trimming; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir after_trimming \
        ${reads_for} ${reads_rev} ${reads_unpaired}
        """
    }

    process PE_map_trimmed_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

      input:
        tuple file(reads_for), file(reads_rev), file(reads_unpaired) from Trimmed_reads_MAP
        file genome_index

      output:
        file "${params.sample_id}.sam" into sam

      script:
        """
        ${HISAT2} -p ${params.threads} \
        --no-spliced-alignment \
        --score-min ${params.relaxed_score} --mp ${params.mm_penalties} \
        --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
        -I ${params.min_isize} -X ${params.max_isize} \
        -x ${genome_index} \
        -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
        -S ${params.sample_id}.sam
        """
    }
  } else {

    process PE_map_untrimmed_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

      input:
        file genome_index

      output:
        file "${params.sample_id}.sam" into sam

      script:
        """
        ${HISAT2} -p ${params.threads} \
        --no-spliced-alignment \
        --score-min ${params.relaxed_score} --mp ${params.mm_penalties} \
        --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
        -I ${params.min_isize} -X ${params.max_isize} \
        -x ${genome_index} \
        -1 ${params.reads_for} -2 ${params.reads_rev} \
        -S ${params.sample_id}.sam
        """
    }
  }


// -----------------------------------------------------------------------------
// SINGLE END READS

} else if (params.read_type == "SE") {

  if (!params.skip_trimming) {

    process SE_quality_check_before_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check", mode: "copy"

      output:
        path "before_trimming"

      script:
        """
        if [ ! -d before_trimming ]; then mkdir before_trimming; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir before_trimming \
        ${params.reads}
        """
    }

    process SE_trim_reads {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir  "${params.output_dir}/trimmed_reads",
                  mode: "copy", pattern: "*.{P1,P2,U1,U2}.fastq"

      output:
        file "${params.sample_id}.trimmed.fastq" into reads
        file "${params.sample_id}.trimmed.fastq" into reads_QC

      script:
        """
        ${TRIMMOMATIC} \
        PE \
        -threads ${params.threads} \
        -summary trim.summary.tsv \
        ${params.reads} ${params.sample_id}.trimmed.fastq \
        ILLUMINACLIP:${params.adapters}:2:30:10 \
        LEADING:${params.leading} \
        TRAILING:${params.trailing} \
        SLIDINGWINDOW:${params.sliding_window} \
        AVGQUAL:${params.avgqual} \
        MINLEN:${params.minlen}
        """
    }

    process SE_quality_check_after_trim {

      executor = "local"
      cpus = params.threads
      maxForks = params.threads

      publishDir "${params.output_dir}/trimmed_reads/quality_check", mode: "copy"

      input:
        file reads_QC

      output:
        path "after_trimming"

      script:
        """
        if [ ! -d after_trimming ]; then mkdir after_trimming; fi &&
        ${FASTQC} \
        --threads ${params.threads} \
        --outdir after_trimming \
        ${reads_QC}
        """
    }
  } else {

    reads = params.reads
  }

  process map_reads {

    executor = "local"
    cpus = params.threads
    maxForks = params.threads

    publishDir  "${params.output_dir}/mapping", mode: "copy", pattern: "*.sam"

    input:
      file reads
      file genome_index

    output:
      file "${params.sample_id}.sam" into sam

    script:
      """
      ${HISAT2} -p ${params.threads} \
      --no-spliced-alignment \
      --score-min ${params.relaxed_score} --mp ${params.mm_penalties} \
      --rdg ${params.gap_penalties} --rfg ${params.gap_penalties} \
      -I ${params.min_isize} -X ${params.max_isize} \
      -x ${genome_index} \
      -U ${reads} \
      -S ${params.sample_id}.sam
      """
    }
  }

process filter_and_sort_bam {

  executor = "local"
  cpus = params.threads
  maxForks = params.threads

  publishDir  "${params.output_dir}/mapping",
              mode: "copy", pattern: "*.{fs.bam,fs.bam.bai}"

  input:
    file sam

  output:
    file "${params.sample_id}.f.bam" into bam_f
    file "${params.sample_id}.fs.bam" into bam_fs
    file "${params.sample_id}.fs.bam" into bam_fs_for_jloh
    file "${params.sample_id}.fs.bam.bai" into bam_fs_index
    file "${params.sample_id}.fs.bam.bai" into bam_fs_index_for_jloh

  script:
    """
    ${SAMTOOLS} view --threads ${params.threads} -F 0x0100 -F 0x4 -h -b \
    -o ${params.sample_id}.f.bam ${sam} &&
    ${SAMTOOLS} sort --threads ${params.threads} --output-fmt bam -T ${params.sample_id} \
    -o ${params.sample_id}.fs.bam \
    ${params.sample_id}.f.bam &&
    ${SAMTOOLS} index -@ ${params.threads} -b \
    ${params.sample_id}.fs.bam
    """
}

// ---------------------------------------------------------------------------
// COMMON STEPS

process pileup_reads {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.mpileup.vcf"

  input:
    file bam_fs
    file bam_fs_index

  output:
    file "${params.sample_id}.mpileup.vcf" into pileup

  script:
    """
    ${SAMTOOLS} faidx ${params.ref_genome} &&
    ${BCFTOOLS} mpileup --fasta-ref ${params.ref_genome} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${params.sample_id}.mpileup.vcf --output-type v --skip-indels \
    ${bam_fs}
    """
}

process call_short_variants {

  executor = "local"
  cpus = params.threads
  maxForks = params.threads

  publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.raw.vcf"

  input:
    file pileup

  output:
    file "${params.sample_id}.raw.vcf" into raw_vcf

  script:
    """
    ${BCFTOOLS} call --threads ${params.threads} --multiallelic-caller --variants-only \
    --output ${params.sample_id}.raw.vcf --output-type v \
    ${pileup}
    """
}

process filter_short_variants {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir  "${params.output_dir}/LOH/process", mode: "copy", pattern: "*.{ff.vcf,report.txt}"

  input:
    file raw_vcf

  output:
    file "${sample_id}.report.txt" into filtering_report
    file "${sample_id}.ff.vcf" into filt_vcf

  script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file ${raw_vcf} \
    --output-file ${params.sample_id}.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_allelefreq} \
    --min-depth ${params.min_depth} \
    --map-qual-zero-frac ${params.mq0f} \
    --threads 1 \
    --report ${params.sample_id}.report.txt &&
    ${ALL2VCF} frequency \
    --in ${params.sample_id}.f.vcf \
    --out ${params.sample_id}.ff.vcf

    """
}

process call_LOH_blocks {

  executor = "local"
  cpus = 1
  maxForks = params.threads

  publishDir "${params.output_dir}/LOH", mode: "copy", pattern: "*.{vcf,tsv,bed}"

  input:
    file filt_vcf
    file bam_fs_for_jloh
    file bam_fs_index_for_jloh

  output:
    file "${params.sample_id}_jloh/*"

  script:
    """
    ${BIOAWK} -c fastx '{print \$name\"\\t\"length(\$seq)}' ${params.ref_genome} \
    > ${params.sample_id}.genome_file.tsv &&
    ${JLOH} \
    --vcf ${filt_vcf} \
    --bam ${bam_fs_for_jloh} \
    --genome-file ${params.sample_id}.genome_file.tsv \
    --sample ${params.sample_id} \
    --output-dir ${params.sample_id}_jloh \
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
