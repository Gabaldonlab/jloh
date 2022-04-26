#!/usr/bin/env nextflow

Channel.fromPath("${params.ref_genome}").set{ Ref_genome }
Channel.from(0.01, 0.03, 0.05, 0.10, 0.15, 0.20).set{ Divergence }
Channel.from(0.1, 0.2, 0.3, 0.4).set{ Loh }
Channel.from(100, 1000).set{ Minlens }
Channel.from(100, 500, 1000).set{ Min_dist }
Channel.from(["30X", 3000000], ["10X", 1000000]).set{ Covs }
Channel.from(2, 4, 6).set{ Snps }

Ref_genome
  .combine(Divergence)
  .combine(Loh)
  .combine(Minlens)
  .combine(Covs)
  .combine(Snps)
  .combine(Min_dist)
  .map{ it -> [ it[0], "snps_${it[6]}.mindist_${it[7]}.minlen_${it[3]}.cov_${it[4]}.run.div_${it[1]}.loh_${it[2]}",
                it[1], it[2], it[3], it[4], it[5], it[6], it[7]] }
  .set{ Test_conditions }


// mutate genome
// introducing LOH blocks and SNPs
// saving coordinates of blocks for true positive evaluation later
// saving SNP coordinates
// saving non-divergent genome regions

process mutate_genome {

  executor="local"
  maxForks=12
  cpus=4
  maxRetries=3

  publishDir "${params.output_dir}/${sample_out_dir}/sim_genome", \
  mode: "copy"

  input:
    tuple file(ref_genome), val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Test_conditions

  output:
    tuple file(ref_genome), \
    file("Sd.fa"), file("Sd.fa.lohs"), file("Sd.fa.non_divergent"), file("Sd.fa.snps"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Mut_ref

  script:
    """
    ${JLOH} sim --fasta ${params.ref_genome} \
    --threads 4 \
    --out Sd.fa \
    --hybrid \
    --divergence ${div} \
    --loh ${loh} \
    --loh-mean-length ${params.loh_mean_length} \
    --chr-name-mod S_cere,S_dive
    """
}

// concatenate genomes
// for IGV viewing later on

process concatenate_genomes {

  executor="local"
  maxForks=48
  cpus=1
  maxRetries=3

  publishDir "${params.output_dir}/${sample_out_dir}/sim_genome", \
  mode: "copy", \
  pattern: "Sc_Sd.*"

  input:
    tuple \
    file(ref_genome), file(mut_genome), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Mut_ref

  output:
    tuple \
    file(ref_genome), file(mut_genome), file("Sc_Sd.fa"), file("Sc_Sd.fa.fai"), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Conc_ref

  script:
    """
    cat ${ref_genome} ${mut_genome} > Sc_Sd.fa
    samtools faidx Sc_Sd.fa
    """
}

// simulate reads
// these reads are the groundwork of the study
// they will be used to attempt recovery of LOH blocks from seq data
// and to calculate precision and recall at the end

process simulate_reads {

  executor="local"
  maxForks=48
  cpus=1
  maxRetries=3

  input:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Conc_ref

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file("Sc_Sd.R1.fastq"), file("Sc_Sd.R2.fastq"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Sim_reads

  script:
    """
    ${WGSIM} \
    -e 0.002 \
    -d 1000 -s 200 \
    -N ${num_reads} \
    -1 125 -2 125 \
    -r 0.0 -R 0.0 \
    -X 0.0 \
    ${conc_genome} \
    Sc_Sd.R1.fastq Sc_Sd.R2.fastq
    """
}

// index genomes and map reads

process map_reads {

  executor="slurm"
  maxForks=16
  cpus=48
  maxRetries=3

  input:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(reads_for), file(reads_rev), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Sim_reads

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(reads_for), file(reads_rev), \
    file("Sc.sam"), file("Sd.sam"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Hisat2_out

  script:
    """
    ${HISAT2_BUILD} \
    -p 48 \
    ${ref_genome} \
    Sc.fa &&
    ${HISAT2_BUILD} \
    -p 48 \
    ${mut_genome} \
    Sd.fa &&
    ${HISAT2} \
    -p 48 \
    --ignore-quals --no-spliced-alignment \
    --score-min L,0.0,-1.0 \
    -I 800 -X 1200 \
    -x Sc.fa -1 ${reads_for} -2 ${reads_rev} \
    -S Sc.sam &&
    ${HISAT2} \
    -p 48 \
    --ignore-quals --no-spliced-alignment \
    --score-min L,0.0,-1.0 \
    -I 800 -X 1200 \
    -x Sd.fa -1 ${reads_for} -2 ${reads_rev} \
    -S Sd.sam
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

  publishDir "${params.output_dir}/${sample_out_dir}/mapping", \
  mode: "copy", \
  pattern: "{Sc,Sd}.fs.{bam,bam.bai}"

  input:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(reads_for), file(reads_rev), \
    file(sc_sam), file(sd_sam), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Hisat2_out

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file("Sc.fs.bam"), file("Sd.fs.bam"), \
    file("Sc.fs.bam.bai"), file("Sd.fs.bam.bai"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Filt_sort_bams

  script:
    """
    ${SAMTOOLS} view -@ 4 -h -b -F 0x0100 -F 0x4 --output Sc.f.bam Sc.sam &&
    ${SAMTOOLS} view -@ 4 -h -b -F 0x0100 -F 0x4 --output Sd.f.bam Sd.sam &&
    ${SAMTOOLS} sort -@ 4 -T Sc -O bam -o Sc.fs.bam Sc.f.bam &&
    ${SAMTOOLS} sort -@ 4 -T Sd -O bam -o Sd.fs.bam Sd.f.bam &&
    ${SAMTOOLS} index Sc.fs.bam &&
    ${SAMTOOLS} index Sd.fs.bam
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
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Filt_sort_bams

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file("Sc.mpileup.vcf"), file("Sd.mpileup.vcf"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Pileups

  script:
    """
    ${BCFTOOLS} mpileup \
    --fasta-ref ${ref_genome} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output Sc.mpileup.vcf \
    --output-type v --skip-indels \
    ${sc_fs_bam} &&
    ${BCFTOOLS} mpileup \
    --fasta-ref ${mut_genome} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output Sd.mpileup.vcf \
    --output-type v --skip-indels \
    ${sd_fs_bam}
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
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file(sc_pileup), file(sd_pileup), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Pileups

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file("Sc.raw.vcf"), file("Sd.raw.vcf"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Var_calls

  script:
    """
    ${BCFTOOLS} call \
    --multiallelic-caller --variants-only \
    --output Sc.raw.vcf --output-type v ${sc_pileup} &&
    ${BCFTOOLS} call \
    --multiallelic-caller --variants-only \
    --output Sd.raw.vcf --output-type v ${sd_pileup}
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

  publishDir "${params.output_dir}/${sample_out_dir}/variants", \
  mode: "copy", \
  pattern: "{Sc,Sd}.ff.vcf"

  input:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file(sc_raw_vcf), file(sd_raw_vcf), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Var_calls

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file("Sc.ff.vcf"), file("Sd.ff.vcf"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Filt_var_calls

  script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file Sc.raw.vcf --map-qual-zero-frac ${params.max_mq0f} \
    --output-file Sc.f.vcf --quality ${params.min_qual} --alt-frac ${params.min_af} --min-depth ${params.min_depth} &&
    ${ALL2VCF} frequency \
    --in Sc.f.vcf --out Sc.ff.vcf &&
    ${ALL2VCF} filter_vcf \
    --input-file Sd.raw.vcf --map-qual-zero-frac ${params.max_mq0f} \
    --output-file Sd.f.vcf --quality ${params.min_qual} --alt-frac ${params.min_af} --min-depth ${params.min_depth} &&
    ${ALL2VCF} frequency \
    --in Sd.f.vcf --out Sd.ff.vcf
    """
}


// extract LOH blocks

process run_jloh_extract {

  executor="local"
  maxForks=1
  cpus=48
  maxRetries=3

  publishDir "${params.output_dir}/${sample_out_dir}/LOH", \
  mode: "copy", \
  pattern: "*.{tsv,bed}"

  input:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file(sc_ff_vcf), file(sd_ff_vcf), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Filt_var_calls

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file(sc_ff_vcf), file(sd_ff_vcf), \
    file("Sc_Sd.LOH_blocks.tsv"), file("Sc_Sd.LOH_blocks.bed"), file("Sc_Sd.LOH_candidates.bed"), \
    file("Sc_Sd.exp_A.het_snps.vcf"), file("Sc_Sd.exp_B.het_snps.vcf"), \
    file("Sc_Sd.exp_A.homo_snps.vcf"), file("Sc_Sd.exp_B.homo_snps.vcf"), \
    file("Sc_Sd.exp_A.chrom_coverage.tsv"), file("Sc_Sd.exp_B.chrom_coverage.tsv"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Jloh_extract_out

  script:
    """
    ${JLOH} extract \
    --hybrid --threads 48 \
    --vcfs ${sc_ff_vcf} ${sd_ff_vcf} \
    --bams ${sc_fs_bam} ${sd_fs_bam} \
    --refs ${ref_genome} ${mut_genome} \
    --sample Sc_Sd --output-dir . \
    --min-length ${min_len} --min-snps ${min_snps} --snp-distance ${min_dist} &&
    { cat Sc_Sd.LOH_blocks.A.tsv; tail -n+2 Sc_Sd.LOH_blocks.B.tsv; } \
    > Sc_Sd.LOH_blocks.tsv &&
    cat Sc_Sd.LOH_blocks.A.bed Sc_Sd.LOH_blocks.B.bed \
    > Sc_Sd.LOH_blocks.bed &&
    cat Sc_Sd.LOH_candidates.A.bed Sc_Sd.LOH_candidates.B.bed > Sc_Sd.LOH_candidates.bed
    """
}


// make statistics

process get_stats_from_run {

  executor="local"
  maxForks=48
  cpus=1
  maxRetries=3

  publishDir "${params.output_dir}/${sample_out_dir}/TP_TN", \
  mode: "copy", \
  pattern: "{generated,predicted,stats}.{bed,tsv}"

  input:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file(sc_ff_vcf), file(sd_ff_vcf), \
    file(blocks_tsv), file(blocks_bed), file(candidates_bed), \
    file(het_snps_vcf_A), file(het_snps_vcf_B), \
    file(homo_snps_vcf_A), file(homo_snps_vcf_B), \
    file(chrom_cov_A), file(chrom_cov_B), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    from Jloh_extract_out

  output:
    tuple \
    file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
    file(true_lohs), file(true_non_divergent), file(true_snps), \
    file(sc_fs_bam), file(sd_fs_bam), \
    file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
    file(sc_ff_vcf), file(sd_ff_vcf), \
    file(blocks_tsv), file(blocks_bed), file(candidates_bed), \
    file(het_snps_vcf_A), file(het_snps_vcf_B), \
    file(homo_snps_vcf_A), file(homo_snps_vcf_B), \
    file(chrom_cov_A), file(chrom_cov_B), \
    file("generated.bed"), file("predicted.bed"), file("stats.tsv"), file("genome_file.tsv"), \
    val(sample_out_dir), val(div), val(loh), \
    val(min_len), val(read_cov), val(num_reads), val(min_snps), val(min_dist) \
    into Stats_out

  script:
    """
    ${BIOAWK} -c fastx '{print \$name"\t"length(\$seq)}' ${conc_genome} > genome_file.tsv &&
    cat ${candidates_bed} | awk '\$1 ~ /S_dive/' > predicted.bed &&
    { cat ${true_non_divergent} | awk -v x=${min_len} '\$3-\$2 >= x'; \
    cat ${true_lohs} | awk -v x=${min_len} '\$3-\$2 >= x'; } | \
    sort -k1,1 -k2n,2 -k3nr,3 > generated.bed &&
    python3 ${TP_TN_SCRIPT} \
    --generated generated.bed \
    --predicted predicted.bed \
    --genome-file genome_file.tsv \
    --out stats.tsv
    """
}
