#!/usr/bin/env nextflow

Channel
    .fromPath(params.input_data)
    .splitCsv(header:true)
    .map{ row-> tuple(row.SAMPLE_ID, row.ACCESSION, file(row.REF_A), file(row.REF_B), file(row.READS_FOR), file(row.READS_REV)) }
    .set { Samples }

// index genomes and map reads

// process map_reads {
//
//   executor="local"
//   maxForks=1
//   cpus=48
//
//   input:
//     tuple \
//     val(sample_id), file(ref_A), file(ref_B), file(reads_for), file(reads_rev) \
//     from Samples
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(reads_for), file(reads_rev), \
//     file("Sc.sam"), file("Sd.sam"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Hisat2_out
//
//   script:
//     """
//     ${HISAT2_BUILD} \
//     -p 48 \
//     ${ref_genome} \
//     Sc.fa &&
//     ${HISAT2_BUILD} \
//     -p 48 \
//     ${mut_genome} \
//     Sd.fa &&
//     ${HISAT2} \
//     -p 48 \
//     --ignore-quals --no-spliced-alignment \
//     --score-min L,0.0,-1.0 \
//     -I 800 -X 1200 \
//     -x Sc.fa -1 ${reads_for} -2 ${reads_rev} \
//     -S Sc.sam &&
//     ${HISAT2} \
//     -p 48 \
//     --ignore-quals --no-spliced-alignment \
//     --score-min L,0.0,-1.0 \
//     -I 800 -X 1200 \
//     -x Sd.fa -1 ${reads_for} -2 ${reads_rev} \
//     -S Sd.sam
//     """
// }
//
//
// // filter and sort read mapping results
// // removing secondary alignments (-F 0x0100)
// // and records from unmapped reads (-F 0x4)
//
// process filter_and_sort_bams {
//
//   executor="local"
//   maxForks=12
//   cpus=4
//
//   publishDir "${params.output_dir}/${sample_out_dir}/mapping", \
//   mode: "copy", \
//   pattern: "{Sc,Sd}.fs.{bam,bam.bai}"
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(reads_for), file(reads_rev), \
//     file(sc_sam), file(sd_sam), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Hisat2_out
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file("Sc.fs.bam"), file("Sd.fs.bam"), \
//     file("Sc.fs.bam.bai"), file("Sd.fs.bam.bai"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Filt_sort_bams
//
//   script:
//     """
//     ${SAMTOOLS} view -@ 4 -h -b -F 0x0100 -F 0x4 --output Sc.f.bam Sc.sam &&
//     ${SAMTOOLS} view -@ 4 -h -b -F 0x0100 -F 0x4 --output Sd.f.bam Sd.sam &&
//     ${SAMTOOLS} sort -@ 4 -T Sc -O bam -o Sc.fs.bam Sc.f.bam &&
//     ${SAMTOOLS} sort -@ 4 -T Sd -O bam -o Sd.fs.bam Sd.f.bam &&
//     ${SAMTOOLS} index Sc.fs.bam &&
//     ${SAMTOOLS} index Sd.fs.bam
//     """
// }
//
// // perform pileup of the reads
// // this is the first step in the variant calling process
// // it also is the longest one
// // and cannot be multi-threaded
// // so this is where we lost the most time
//
// process perform_pileup {
//
//   executor="local"
//   maxForks=48
//   cpus=1
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Filt_sort_bams
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file("Sc.mpileup.vcf"), file("Sd.mpileup.vcf"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Pileups
//
//   script:
//     """
//     ${BCFTOOLS} mpileup \
//     --fasta-ref ${ref_genome} \
//     --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
//     --output Sc.mpileup.vcf \
//     --output-type v --skip-indels \
//     ${sc_fs_bam} &&
//     ${BCFTOOLS} mpileup \
//     --fasta-ref ${mut_genome} \
//     --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
//     --output Sd.mpileup.vcf \
//     --output-type v --skip-indels \
//     ${sd_fs_bam}
//     """
// }
//
// // call variants
// // here only SNPs are kept, as indels were not part of the simulation
//
// process perform_snp_calling {
//
//   executor="local"
//   maxForks=48
//   cpus=1
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_pileup), file(sd_pileup), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Pileups
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file("Sc.raw.vcf"), file("Sd.raw.vcf"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Var_calls
//
//   script:
//     """
//     ${BCFTOOLS} call \
//     --multiallelic-caller --variants-only \
//     --output Sc.raw.vcf --output-type v ${sc_pileup} &&
//     ${BCFTOOLS} call \
//     --multiallelic-caller --variants-only \
//     --output Sd.raw.vcf --output-type v ${sd_pileup}
//     """
// }
//
// // filter SNPs
// // removing those with:
// // > 0.05 MQ0F
// // < 20 QUAL
// // < 4 DP
// // < 0.05 AF
//
// process filter_variants {
//
//   executor="local"
//   maxForks=12
//   cpus=4
//
//   publishDir "${params.output_dir}/${sample_out_dir}/variants", \
//   mode: "copy", \
//   pattern: "{Sc,Sd}.ff.vcf"
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_raw_vcf), file(sd_raw_vcf), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Var_calls
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file("Sc.ff.vcf"), file("Sd.ff.vcf"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Filt_var_calls
//
//   script:
//     """
//     ${ALL2VCF} filter_vcf \
//     --input-file Sc.raw.vcf --map-qual-zero-frac ${params.max_mq0f} \
//     --output-file Sc.f.vcf --quality ${params.min_qual} --alt-frac ${params.min_af} --min-depth ${params.min_depth} &&
//     ${ALL2VCF} frequency \
//     --in Sc.f.vcf --out Sc.ff.vcf &&
//     ${ALL2VCF} filter_vcf \
//     --input-file Sd.raw.vcf --map-qual-zero-frac ${params.max_mq0f} \
//     --output-file Sd.f.vcf --quality ${params.min_qual} --alt-frac ${params.min_af} --min-depth ${params.min_depth} &&
//     ${ALL2VCF} frequency \
//     --in Sd.f.vcf --out Sd.ff.vcf
//     """
// }
//
// // perform genome to genome alignment
// // to exclude regions where the two reference genomes
// // are identical to begin with
// // which are regions where we can expect no heterozygosity to happen
// // this doesn't matter with simulated data
// // since we do a copy of the S. cerevisiae genome to get the Sd genome
// // every part that looks exactly the same is bound to be an LOH block
// // this will matter with real data where we cannot expect regions to be
// // exactly identical
//
// process perform_g2g_mapping {
//
//   executor="local"
//   maxForks=48
//   cpus=1
//
//   publishDir "${params.output_dir}/${sample_out_dir}/LOH", \
//   mode: "copy", \
//   pattern: "Sc_Sd.regions.bed"
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_ff_vcf), file(sd_ff_vcf), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Filt_var_calls
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_ff_vcf), file(sd_ff_vcf), file("Sc_Sd.regions.bed"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Regions_out
//
//   script:
//     """
//     ${JLOH} g2g \
//     --sensitive \
//     --ref-A ${ref_genome} --ref-B ${mut_genome} --est-divergence ${div} \
//     --out . \
//     --min-length ${min_len} &&
//     cat A.bed B.bed > Sc_Sd.regions.bed
//     """
// }
//
//
// //
//
// process run_jloh_extract {
//
//   executor="local"
//   maxForks=1
//   cpus=48
//
//   publishDir "${params.output_dir}/${sample_out_dir}/LOH", \
//   mode: "copy", \
//   pattern: "*.{tsv,bed}"
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_ff_vcf), file(sd_ff_vcf), file(regions_file), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Regions_out
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_ff_vcf), file(sd_ff_vcf), file(regions_file), \
//     file("Sc_Sd.LOH_blocks.tsv"), file("Sc_Sd.LOH_blocks.bed"), file("Sc_Sd.LOH_candidates.bed"), \
//     file("Sc_Sd.exp_A.het_snps.vcf"), file("Sc_Sd.exp_B.het_snps.vcf"), \
//     file("Sc_Sd.exp_A.homo_snps.vcf"), file("Sc_Sd.exp_B.homo_snps.vcf"), \
//     file("Sc_Sd.exp_A.chrom_coverage.tsv"), file("Sc_Sd.exp_B.chrom_coverage.tsv"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Jloh_extract_out
//
//   script:
//     """
//     ${JLOH} extract \
//     --hybrid --threads 48 \
//     --vcfs ${sc_ff_vcf} ${sd_ff_vcf} \
//     --bams ${sc_fs_bam} ${sd_fs_bam} \
//     --refs ${ref_genome} ${mut_genome} \
//     --sample Sc_Sd --output-dir . --regions ${regions_file} \
//     --min-length ${min_len} --min-snps ${min_snps} &&
//     { cat Sc_Sd.LOH_blocks.A.tsv; tail -n+2 Sc_Sd.LOH_blocks.B.tsv; } \
//     > Sc_Sd.LOH_blocks.tsv &&
//     cat Sc_Sd.LOH_blocks.A.bed Sc_Sd.LOH_blocks.B.bed \
//     > Sc_Sd.LOH_blocks.bed &&
//     cat Sc_Sd.LOH_candidates.A.bed Sc_Sd.LOH_candidates.B.bed > Sc_Sd.LOH_candidates.bed
//     """
// }
//
//
// // make statistics
//
// process get_stats_from_run {
//
//   executor="local"
//   maxForks=48
//   cpus=1
//
//   publishDir "${params.output_dir}/${sample_out_dir}/TP_FP", \
//   mode: "copy", \
//   pattern: "{generated,predicted,stats}.{bed,tsv}"
//
//   input:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_ff_vcf), file(sd_ff_vcf), file(regions_file), \
//     file(blocks_tsv), file(blocks_bed), file(candidates_bed), \
//     file(het_snps_vcf_A), file(het_snps_vcf_B), \
//     file(homo_snps_vcf_A), file(homo_snps_vcf_B), \
//     file(chrom_cov_A), file(chrom_cov_B), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     from Jloh_extract_out
//
//   output:
//     tuple \
//     file(ref_genome), file(mut_genome), file(conc_genome), file(conc_genome_idx), \
//     file(true_lohs), file(true_non_divergent), file(true_snps), \
//     file(sc_fs_bam), file(sd_fs_bam), \
//     file(sc_fs_bam_idx), file(sd_fs_bam_idx), \
//     file(sc_ff_vcf), file(sd_ff_vcf), file(regions_file), \
//     file(blocks_tsv), file(blocks_bed), file(candidates_bed), \
//     file(het_snps_vcf_A), file(het_snps_vcf_B), \
//     file(homo_snps_vcf_A), file(homo_snps_vcf_B), \
//     file(chrom_cov_A), file(chrom_cov_B), \
//     file("generated.bed"), file("predicted.bed"), file("stats.tsv"), file("genome_file.tsv"), \
//     val(sample_out_dir), val(div), val(loh), \
//     val(min_len), val(read_cov), val(num_reads), val(min_snps) \
//     into Stats_out
//
//   script:
//     """
//     ${BIOAWK} -c fastx '{print \$name"\t"length(\$seq)}' ${conc_genome} > genome_file.tsv &&
//     cat ${candidates_bed} | awk '\$1 ~ /S_dive/' > predicted.bed &&
//     { cat ${true_non_divergent} | awk -v x=${min_len} '\$3-\$2 >= x'; \
//     cat ${true_lohs} | awk -v x=${min_len} '\$3-\$2 >= x'; } | \
//     sort -k1,1 -k2n,2 -k3nr,3 > generated.bed &&
//     python3 ${TP_FP_SCRIPT} \
//     --generated generated.bed \
//     --predicted predicted.bed \
//     --genome-file genome_file.tsv \
//     --out stats.tsv
//     """
// }
