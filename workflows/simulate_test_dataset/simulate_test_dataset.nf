#!/usr/bin/env nextflow 

// simulate genome at 10% divergence 
// ---------------------------------
// only divergence 
// no LOH 
// use jloh sim 

process simulate_div_genome {

    executor "local"
    cpus 48
    maxForks 1
    maxRetries 3

    output:
    file "sim.fa" into sim_genome 

    script:
    """
    ${JLOH} sim \
    --threads 48 \
    --fasta ${params.fasta} \
    --out-fasta sim.fa \
    --out-haplotypes sim.haplotypes.tsv \
    --mean-haplotype-size 25000 \
    --min-haplotype-length 5000 \
    --divergence ${params.divergence} \
    --loh 0.0 \
    --chrom-name-replace S_cere S_dive
    """
}


// introduce some LOH between the two generated genomes 
// ----------------------------------------------------
// this is not done with jloh sim 
// because it requires just swapping a few portions
// make sure to save true positives 

process introduce_LOH {

    executor "local"
    cpus 48
    maxForks 1
    maxRetries 3

    publishDir "${params.output_dir}", 
    mode: "copy", 
    pattern: "*.{fa,fa.true_loh_blocks.bed}"

    input:
    file sim_genome

    output:
    tuple \
    file("A.fa"), file("B.fa") \
    into Genomes_w_loh

    tuple \
    file("A.fa.true_loh_blocks.bed"), file("B.fa.true_loh_blocks.bed") \
    into True_loh_blocks 

    script:
    """
    ${projectDir}/src/swap-blocks.py \
    --genome-A ${params.fasta} \
    --genome-B ${sim_genome} \
    --loh-fraction ${params.loh_fraction} \
    --mean-haplotype-size 25000 \
    --min-haplotype-length 5000 \
    --out-A A.fa \
    --out-B B.fa
    """
}


// simulate reads 
// --------------
// using wgsim
// with some extra mutation rate (0.5%)
// to add some noise 

process simulate_reads {

    executor "local"
    maxForks 48
    cpus 1
    maxRetries 3

    publishDir "${params.output_dir}",
    mode: "copy", 
    pattern: "hybrid.{R1,R2}.fastq"

    input:
    tuple \
    file(genome_A), file(genome_B) \
    from Genomes_w_loh

    output:
    tuple \
    file(genome_A), file(genome_B), \
    file("hybrid.R1.fastq"), file("hybrid.R2.fastq") \
    into Read_sim_out

    script:
    """
    cat ${genome_A} ${genome_B} > hybrid.fa &&
    ${WGSIM} \
    -d 1000 -s 200 \
    -N ${params.num_reads} \
    -1 ${params.read_length} -2 ${params.read_length} \
    -r ${params.mutation_rate} \
    hybrid.fa \
    hybrid.R1.fastq hybrid.R2.fastq
    """
}


// index genomes and map reads
// ---------------------------
// using hisat2 
// to make it fast 
// for the testing 

process map_reads {

    executor "local"
    maxForks 1
    cpus 48
    maxRetries 3

    input:
    tuple \
    file(genome_A), file(genome_B), \
    file(R1), file(R2) \
    from Read_sim_out

    output:
    tuple \
    file(genome_A), file(genome_B), \
    file("hybrid.onto_A.sam"), file("hybrid.onto_B.sam") \
    into Mapping_out

    script:
    """
    ${HISAT2_BUILD} \
    -p 48 \
    ${genome_A} \
    A_index &&
    ${HISAT2_BUILD} \
    -p 48 \
    ${genome_B} \
    B_index &&
    ${HISAT2} \
    -p 48 \
    --ignore-quals --no-spliced-alignment \
    --score-min L,0.0,-2.0 \
    -I 800 -X 1200 \
    -x A_index -1 ${R1} -2 ${R2} \
    -S hybrid.onto_A.sam &&
    ${HISAT2} \
    -p 48 \
    --ignore-quals --no-spliced-alignment \
    --score-min L,0.0,-2.0 \
    -I 800 -X 1200 \
    -x B_index -1 ${R1} -2 ${R2} \
    -S hybrid.onto_B.sam
    """
}


// filter and sort read mapping results
// ------------------------------------
// removing secondary alignments (-F 0x0100)
// and records from unmapped reads (-F 0x4)

process filter_and_sort_bams {

    executor "local"
    maxForks 1
    cpus 48
    maxRetries 3

    publishDir "${params.output_dir}",
    mode: "copy",
    pattern: "*.fs.{bam,bam.bai}"

    input:
    tuple \
    file(genome_A), file(genome_B), \
    file(genome_A_sam), file(genome_B_sam) \
    from Mapping_out

    output:
    tuple \
    file(genome_A), file(genome_B), \
    file("hybrid.onto_A.fs.bam"), file("hybrid.onto_A.fs.bam.bai"), \
    file("hybrid.onto_B.fs.bam"), file("hybrid.onto_B.fs.bam.bai") \
    into Filt_sort_bams

    script:
    """
    ${SAMTOOLS} view -@ 48 -h -b -F 0x0100 -F 0x4 \
    --output hybrid.onto_A.f.bam ${genome_A_sam} &&
    ${SAMTOOLS} view -@ 48 -h -b -F 0x0100 -F 0x4 \
    --output hybrid.onto_B.f.bam ${genome_B_sam} &&
    ${SAMTOOLS} sort -@ 48 -T A -O bam -o hybrid.onto_A.fs.bam \
    hybrid.onto_A.f.bam &&
    ${SAMTOOLS} sort -@ 48 -T B -O bam -o hybrid.onto_B.fs.bam \
    hybrid.onto_B.f.bam &&
    ${SAMTOOLS} index hybrid.onto_A.fs.bam &&
    ${SAMTOOLS} index hybrid.onto_B.fs.bam
    """
}


// perform pileup of the reads
// ---------------------------
// annotating the relevant fields that are required from JLOH 

process perform_pileup {

    executor "local"
    maxForks 1
    cpus 48
    maxRetries 3

    input:
    tuple \
    file(genome_A), file(genome_B), \
    file(A_bam), file(A_bai), \
    file(B_bam), file(B_bai) \
    from Filt_sort_bams

    output:
    tuple \
    file(genome_A), file(genome_B), \
    file(A_bam), file(A_bai), \
    file(B_bam), file(B_bai), \
    file("hybrid.onto_A.mpileup.vcf"), file("hybrid.onto_B.mpileup.vcf") \
    into Pileups

    script:
    """
    ${BCFTOOLS} mpileup \
    --fasta-ref ${genome_A} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output hybrid.onto_A.mpileup.vcf \
    --output-type v --skip-indels \
    ${A_bam} &&
    ${BCFTOOLS} mpileup \
    --fasta-ref ${genome_B} \
    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
    --output hybrid.onto_B.mpileup.vcf \
    --output-type v --skip-indels \
    ${B_bam}
    """
}


// call variants
// -------------

process perform_snp_calling {

    executor "local"
    maxForks 1
    cpus 48
    maxRetries 3

    input:
    tuple \
    file(genome_A), file(genome_B), \
    file(A_bam), file(A_bai), \
    file(B_bam), file(B_bai), \
    file(A_pileup), file(B_pileup) \
    from Pileups

    output:
    tuple \
    file(genome_A), file(genome_B), \
    file(A_bam), file(A_bai), \
    file(B_bam), file(B_bai), \
    file("hybrid.onto_A.raw.vcf"), file("hybrid.onto_B.raw.vcf") \
    into Var_calls

    script:
    """
    ${BCFTOOLS} call \
    --multiallelic-caller --variants-only \
    --output hybrid.onto_A.raw.vcf --output-type v ${A_pileup} &&
    ${BCFTOOLS} call \
    --multiallelic-caller --variants-only \
    --output hybrid.onto_B.raw.vcf --output-type v ${B_pileup} 
    """
}


// filter SNPs
// -----------

process filter_variants {

    executor "local"
    maxForks 1
    cpus 48
    maxRetries 3

    publishDir "${params.output_dir}",
    mode: "copy",
    pattern: "*.ff.vcf"

    input:
    tuple \
    file(genome_A), file(genome_B), \
    file(A_bam), file(A_bai), \
    file(B_bam), file(B_bai), \
    file(A_vcf), file(B_vcf) \
    from Var_calls

    output:
    tuple \
    file(genome_A), file(genome_B), \
    file(A_bam), file(A_bai), \
    file(B_bam), file(B_bai), \
    file("hybrid.onto_A.ff.vcf"), file("hybrid.onto_B.ff.vcf") \
    into Filt_var_calls

    script:
    """
    ${ALL2VCF} filter_vcf \
    --input-file ${A_vcf} \
    --map-qual-zero-frac ${params.max_mq0f} \
    --output-file hybrid.onto_A.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_af} \
    --min-depth ${params.min_depth} &&
    ${ALL2VCF} frequency \
    --in hybrid.onto_A.f.vcf --out hybrid.onto_A.ff.vcf &&
    ${ALL2VCF} filter_vcf \
    --input-file ${B_vcf} \
    --map-qual-zero-frac ${params.max_mq0f} \
    --output-file hybrid.onto_B.f.vcf \
    --quality ${params.min_qual} \
    --alt-frac ${params.min_af} \
    --min-depth ${params.min_depth} &&
    ${ALL2VCF} frequency \
    --in hybrid.onto_B.f.vcf --out hybrid.onto_B.ff.vcf
    """
}
