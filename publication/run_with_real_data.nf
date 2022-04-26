#!/usr/bin/env nextflow

Channel
  .fromPath(params.input_data)
  .splitCsv(header: true, sep: "\t")
  .map{ it -> [it.SAMPLE_ID, it.ACCESSION, it.REF_A, it.REF_B, it.READS_FOR, it.READS_REV] }
  .map{ it -> [it[0], it[1], file(it[2]), file(it[3]), file(it[2]).baseName, file(it[3]).baseName, file(it[4]), file(it[5])] }
  .set { Samples }

// index genomes and map reads

process trim_reads {

  executor="slurm"
  maxForks=4
  cpus=12
  maxRetries=3
  clusterOptions="--ntasks-per-node=4"

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(reads_for), file(reads_rev) \
    from Samples

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file("${sample_id}.P1.fastq"), file("${sample_id}.P2.fastq"), file("${sample_id}.U.fastq") \
    into Trimmed_reads

  script:
    """
    ${TRIMMOMATIC} \
    PE \
    -threads 12 \
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

process map_reads {

  executor="slurm"
  maxForks=1
  cpus=48
  maxRetries=3
  clusterOptions="--ntasks-per-node=1"

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(reads_for), file(reads_rev), file(reads_unpaired) \
    from Trimmed_reads

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file("${sample_id}.${ref_A_name}.sam"), file("${sample_id}.${ref_B_name}.sam") \
    into Hisat2_out

  script:
    """
    ${HISAT2_BUILD} -p 48 ${ref_A} ${ref_A_name} &&
    ${HISAT2_BUILD} -p 48 ${ref_B} ${ref_B_name} &&
    ${HISAT2} -p 48 \
    --no-spliced-alignment \
    --score-min L,0.0,-1.0 \
    -I 0 -X 1000 \
    -x ${ref_A_name} -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
    -S ${sample_id}.${ref_A_name}.sam &&
    ${HISAT2} -p 48 \
    --no-spliced-alignment \
    -I 0 -X 1000 \
    --score-min L,0.0,-1.0 \
    -x ${ref_B_name} -1 ${reads_for} -2 ${reads_rev} -U ${reads_unpaired} \
    -S ${sample_id}.${ref_B_name}.sam
    """
}


// filter and sort read mapping results
// removing secondary alignments (-F 0x0100)
// and records from unmapped reads (-F 0x4)

process filter_and_sort_bams {

  executor="local"
  maxForks=12
  cpus=4
  maxRetries=3

  publishDir "${params.output_dir}/${sample_id}/mapping", \
  mode: "copy", \
  pattern: "${sample_id}.{${ref_A_name},${ref_B_name}}.fs.{bam,bam.bai}"

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(sam_A), file(sam_B) \
    from Hisat2_out

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file("${sample_id}.${ref_A_name}.fs.bam"), file("${sample_id}.${ref_B_name}.fs.bam") \
    into Filt_sort_bams

  script:
    """
    ${SAMTOOLS} view -@ 4 -h -b -F 0x0100 -F 0x4 --output ${sample_id}.${ref_A_name}.f.bam ${sam_A} &&
    ${SAMTOOLS} view -@ 4 -h -b -F 0x0100 -F 0x4 --output ${sample_id}.${ref_B_name}.f.bam ${sam_B} &&
    ${SAMTOOLS} sort -@ 4 -T Sc -O bam -o ${sample_id}.${ref_A_name}.fs.bam ${sample_id}.${ref_A_name}.f.bam &&
    ${SAMTOOLS} sort -@ 4 -T Sc -O bam -o ${sample_id}.${ref_B_name}.fs.bam ${sample_id}.${ref_B_name}.f.bam &&
    ${SAMTOOLS} index ${sample_id}.${ref_A_name}.fs.bam &&
    ${SAMTOOLS} index ${sample_id}.${ref_B_name}.fs.bam
    """
}

// perform pileup of the reads
// this is the first step in the variant calling process
// it also is the longest one
// and cannot be multi-threaded
// so this is where we lost the most time

process perform_pileup {

  executor="local"
  maxForks=48
  cpus=1
  maxRetries=3

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B) \
    from Filt_sort_bams

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B), \
    file("${sample_id}.${ref_A_name}.mpileup.vcf"), file("${sample_id}.${ref_B_name}.mpileup.vcf") \
    into Pileups

  script:
    """
    ${BCFTOOLS} mpileup \
    --fasta-ref ${ref_A} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${sample_id}.${ref_A_name}.mpileup.vcf \
    --output-type v --skip-indels \
    ${fs_bam_A} &&
    ${BCFTOOLS} mpileup \
    --fasta-ref ${ref_B} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output ${sample_id}.${ref_B_name}.mpileup.vcf \
    --output-type v --skip-indels \
    ${fs_bam_B}
    """
}

// call variants
// here only SNPs are kept, as indels were not part of the simulation

process perform_snp_calling {

  executor="local"
  maxForks=48
  cpus=1
  maxRetries=3

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B), \
    file(pileup_A), file(pileup_B) \
    from Pileups

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B), \
    file("${sample_id}.${ref_A_name}.raw.vcf"), file("${sample_id}.${ref_B_name}.raw.vcf") \
    into Var_calls

  script:
    """
    ${BCFTOOLS} call \
    --multiallelic-caller --variants-only \
    --output ${sample_id}.${ref_A_name}.raw.vcf --output-type v ${pileup_A} &&
    ${BCFTOOLS} call \
    --multiallelic-caller --variants-only \
    --output ${sample_id}.${ref_B_name}.raw.vcf --output-type v ${pileup_B}
    """
}

// filter SNPs
// removing those with:
// > 0.05 MQ0F
// < 20 QUAL
// < 4 DP
// < 0.05 AF

process filter_variants {

  executor="local"
  maxForks=48
  cpus=1
  maxRetries=3

  publishDir "${params.output_dir}/${sample_id}/variants", \
  mode: "copy", \
  pattern: "${sample_id}.{${ref_A_name},${ref_B_name}}.ff.vcf"

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B), \
    file(raw_vcf_A), file(raw_vcf_B) \
    from Var_calls

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B), \
    file("${sample_id}.${ref_A_name}.ff.vcf"), file("${sample_id}.${ref_B_name}.ff.vcf") \
    into Filt_var_calls

  script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file ${raw_vcf_A} --map-qual-zero-frac ${params.max_mq0f} \
    --output-file ${sample_id}.${ref_A_name}.f.vcf --quality ${params.min_qual} --alt-frac ${params.min_af} --min-depth ${params.min_depth} &&
    ${ALL2VCF} frequency \
    --in ${sample_id}.${ref_A_name}.f.vcf --out ${sample_id}.${ref_A_name}.ff.vcf &&
    ${ALL2VCF} filter_vcf \
    --input-file ${raw_vcf_B} --map-qual-zero-frac ${params.max_mq0f} \
    --output-file ${sample_id}.${ref_B_name}.f.vcf --quality ${params.min_qual} --alt-frac ${params.min_af} --min-depth ${params.min_depth} &&
    ${ALL2VCF} frequency \
    --in ${sample_id}.${ref_B_name}.f.vcf --out ${sample_id}.${ref_B_name}.ff.vcf
    """
}

// extract LOH blocks

process run_jloh_extract {

  // executor="slurm"
  // maxForks=2
  // cpus=24
  // maxRetries=3
  // clusterOptions="--ntasks-per-node=2"

  executor="local"
  maxForks=1
  cpus=48

  publishDir "${params.output_dir}/${sample_id}/LOH", \
  mode: "copy", \
  pattern: "*.{tsv,bed}"

  input:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    val(ref_A_name), val(ref_B_name), \
    file(fs_bam_A), file(fs_bam_B), \
    file(ff_vcf_A), file(ff_vcf_B) \
    from Filt_var_calls

  output:
    tuple \
    val(sample_id), val(accession), file(ref_A), file(ref_B), \
    file("${sample_id}.LOH_blocks.tsv"), file("${sample_id}.LOH_blocks.bed"), \
    file("${sample_id}.LOH_candidates.tsv"), file("${sample_id}.LOH_candidates.bed") \
    into Jloh_extract_out

  script:
    """
    ${JLOH} extract \
    --hybrid --threads 48 \
    --vcfs ${ff_vcf_A} ${ff_vcf_B} \
    --bams ${fs_bam_A} ${fs_bam_B} \
    --refs ${ref_A} ${ref_B} \
    --sample ${sample_id} --output-dir . \
    --min-length ${params.min_len} --min-snps ${params.min_snps} &&
    { cat ${sample_id}.LOH_blocks.A.tsv; tail -n+2 ${sample_id}.LOH_blocks.B.tsv; } \
    > ${sample_id}.LOH_blocks.tsv &&
    { cat ${sample_id}.LOH_candidates.A.tsv; tail -n+2 ${sample_id}.LOH_candidates.B.tsv; } \
    > ${sample_id}.LOH_candidates.tsv &&
    cat ${sample_id}.LOH_blocks.A.bed ${sample_id}.LOH_blocks.B.bed \
    > ${sample_id}.LOH_blocks.bed &&
    cat ${sample_id}.LOH_candidates.A.bed ${sample_id}.LOH_candidates.B.bed > ${sample_id}.LOH_candidates.bed
    """
}
