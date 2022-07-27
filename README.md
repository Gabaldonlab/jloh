# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract, filter, and manage blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome. Optionally works with a paired set of VCFs, BAMs and FASTAs produced from a hybrid genome (`--hybrid` mode).

![JLOH workflow](images/j_loh.png)

## Table of Contents

- [J LOH](#j-loh)
  * [Table of Contents](#table-of-contents)
  * [Install](#install)
  * [Run](#run)
    + [Nextflow workflow](#nextflow-workflow)
- [Modules](#modules)
  * [JLOH sim](#jloh-sim)
  * [JLOH g2g](#jloh-g2g)
  * [JLOH extract](#jloh-extract)
    + [Sorting of SNPs by zygosity](#sorting-of-snps-by-zygosity)
    + [Extraction of heterozygous regions](#extraction-of-heterozygous-regions)
    + [Extraction of candidate blocks](#extraction-of-candidate-blocks)
    + [Coverage trimming](#coverage-trimming)
    + [Determination of block zygosity](#determination-of-block-zygosity)
    + [Output](#output)
      - [Interpreting output](#interpreting-output)
  * [JLOH filter](#jloh-filter)
  * [JLOH density](#jloh-density)
- [Hybrid mode - step by step guide](#hybrid-mode---step-by-step-guide)


## Install and run

Detailed information [can be found here](docs/INSTALL.md).

## JLOH extract

This is the most important module of JLOH. Its functions are well described in the figure on the top. It is used to extract LOH blocks starting from VCF, BAM, and FASTA files. Detailed information on this algorithm functioning [can be found here](docs/EXTRACT.md).

### output

`jloh extract` functions in two modes (hybrid, and default). In default mode, it produces these output files:

| Output file                 | Description |
|-----------------------------|-|
| jloh.exp.het_snps.vcf       | VCF file containing the SNPs that were labelled as "heterozygous" by the algorithm |
| jloh.exp.homo_snps.vcf      | VCF file containing the SNPs that were labelled as "homozygous" by the algorithm |
| jloh.exp.het_blocks.bed     | BED file with heterozygous blocks, i.e. intervals in the genome carrying heterozygous SNPs |
| jloh.exp.genome_file.tsv    | TSV file with chromosomes and their length |
| jloh.exp.chrom_coverage.tsv | TSV file with chromosomes and their average coverage per position |
| jloh.LOH_candidates.bed     | BED file with candidate LOH blocks before coverage and length screening |
| jloh.LOH_candidates.tsv     | TSV file with candidate LOH blocks before coverage and length screening |
| jloh.LOH_blocks.bed         | BED file with selected LOH blocks |
| jloh.LOH_blocks.tsv         | TSV file with selected LOH blocks. *This is the main output of the program* |

The output can be assessed in a genome viewer together with the input BAM files, producing a profile like this one:

![Example](images/example.png)

## JLOH filter

This tool filters the output produced by `JLOH extract` according to several criteria. The user can select to filter the LOH blocks based on their coverage, on their SNP count, SNP density, length, or extract individual regions just like in samtools.

## JLOH intersect

This tool perform ensemble operations with two JLOH output files, namely intersection, complement, and unique elements extraction.

## JLOH chimeric

This tool is a specific module to extract genes that overlap LOH blocks from two different origins (i.e. chimeric genes). It does not assume anything about the lists of LOH blocks that are passed. If a gene has two LOH blocks in its sequence (one from each list) it will be considered a candidate chimeric gene.

The usage involves two sets of LOH blocks produced by `jloh extract`, plus the `*het_blocks.bed` file also produced by `jloh extract`. In case the "extract" module was run in `--hybrid` mode, the heterozygous blocks files are two (A and B). In that case, they must be concatenated into a single file using `cat` and provided as a single file to `jloh chimeric`.

### output

`jloh chimeric` produces two files:

| Output file               | Description |
|---------------------------|-|
| out.chimeric.features.gff | GFF file containing rows that identify genes marked as chimeric. This file can be directly loaded into IGV or GBrowse for manual inspection, making your job easier. |
| out.chimeric.IDs.txt      | TXT file containing the gene IDs of all the genes found as chimeric. These can be of three types: 1) A-B chimeras, from a mixture of the two LOH block sources; 2) A-H chimeras, from a mixture of the A source and the heterozygous blocks, 3) B-H chimeras, from a mixture of the B source and the heterozygous blocks. |

## JLOH g2g

This program finds diverging regions between two genomes in FASTA format. The input are the two sequences and an estimated divergence value, and the output is a bed file representing the regions that contain SNPs between the two genomes.

`JLOH g2g` runs more than one tool from the MUMmer arsenal to map the two genomes, filter the results, extract the SNPs. Then, it uses `all2vcf mummer` to convert the MUMmer output to VCF format, and `bedtools merge` to generate BED intervals from SNPs. Intervals are expanded as long as there are overlaps.

These regions are a good `--regions` file to pass to `JLOH extract` in `--hybrid` mode, excluding false positives that may arise by mapping.

## JLOH density

This tool computes the densities of all SNPs, heterozygous SNPs, and homozygous SNPs over the genome sequence.

## JLOH sim

This module generates a copy of a reference sequence in FASTA format, introducing a series of mutations selected randomly over the sequence of each scaffold/chromosome. Optionally, the module can include a series of LOH blocks defined by percentage of the whole genome (e.g. 20%).

The output is:
- the mutated reference sequence in FASTA format
- a tab-separated file with all the introduced SNPs, with positions annotated in 1-based format including reference and alternative allele
- a BED file with all the introduced LOH blocks, hence in 0-based half-open format
- another BED file with all the regions that are deemed as non-divergent based on the declared `--snp-distance`.

# Hybrid mode - step by step guide for data preparation

*The following steps (apart from the installation) are automatically performed by the nextflow workflow described in the previous section. Maybe have a look at that before doing everything by hand!* ;)

The first step is making sure that everything has been installed properly. Hence, follow the instructions in the [Install](#install) section prior to trying this guide.

The hybrid mode is used to extract LOH blocks in hybrid genomes. Hybrid genomes have two subgenomes that carry substantial differences between each other, while still being somewhat close in terms of sequence identity. Single-nucleotide positional differences between two subgenomes can be identified as **heterozygous SNPs** if reads are mapped against one subgenome at a time, allowing reads from both subgenomes to map onto it. Regions where little heterozygous SNPs are found are good LOH block candidates.

We start by mapping a set of quality-trimmed paired-end reads sequenced from a hybrid species against both of its parental genome sequences. We do it with **HISAT2** but you can use the tool you prefer. For the whole paragraph we will use four threads, represented here as `-p 4`.

```
hisat2-build -p 4 parent_A.fasta idx_A
hisat2-build -p 4 parent_B.fasta idx_B

hisat2-align-s -p 4 --score-min L,0.0,-1.0 -x idx_A -1 hybrid_reads.R1.fastq -2 hybrid_reads.R2.fastq -S A.sam
hisat2-align-s -p 4 --score-min L,0.0,-1.0 -x idx_B -1 hybrid_reads.R1.fastq -2 hybrid_reads.R2.fastq -S B.sam
```

Note the particularly relaxed mapping criteria (`--score-min L,0.0,-1.0`). These will allow for most reads to map on both genomes, regardless of which subgenome they come from. The result is that we can simulate heterozygosity in the derived SNPs. When mapping the reads, make sure that you're mapping most of them (85-90%). A fraction may still not map due to the target region being missing from the genome sequence.

After mapping, we filter, sort and index the output **sam** files using **samtools** with 4 threads.

```
samtools view -@ 4 -h -b -F 0x0100 -F 0x4 A.sam | samtools sort -@ 4 -T tmp_A > A.bam
samtools view -@ 4 -h -b -F 0x0100 -F 0x4 B.sam | samtools sort -@ 4 -T tmp_B > B.bam

samtools index A.bam
samtools index B.bam
```

We then index the reference genomes using **samtools** and we use the filtered mapping records to perform a pileup of the reads using **bcftools**.

```
samtools faidx parent_A.fasta
samtools faidx parent_B.fasta

bcftools mpileup --fasta-ref parent_A.fasta \
--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
--output-type v --skip-indels --output A.mpileup.vcf A.bam

bcftools mpileup --fasta-ref parent_B.fasta \
--annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
--output-type v --skip-indels --output B.mpileup.vcf B.bam
```

Note that we don't care about indels in this workflow, using the `--skip-indels` option, and that we care instead to annotate a lot of INFO/ and FORMAT/ fields of the generated VCF file. We then call the SNPs using **bcftools** again.

```
bcftools call --threads 4 --multiallelic-caller --variants-only \
--output A.raw.vcf --output-type v A.mpileup.vcf

bcftools call --threads 4 --multiallelic-caller --variants-only \
--output B.raw.vcf --output-type v B.mpileup.vcf
```

Once we have called the SNPs, we filter them using **all2vcf** and we add the allele frequency subfield in the INFO field also using **all2vcf**.

```
all2vcf filter_vcf --input-file A.raw.vcf --output-file A.f.vcf --quality 20 --alt-frac 0.05 \
--min-depth 4 --map-qual-zero-frac 0.05
all2vcf filter_vcf --input-file B.raw.vcf --output-file B.f.vcf --quality 20 --alt-frac 0.05 \
--min-depth 4 --map-qual-zero-frac 0.05

all2vcf frequency --in A.f.vcf --out A.ff.vcf
all2vcf frequency --in B.f.vcf --out B.ff.vcf
```

Now, we extract the regions containing heterozygous SNPs between the two parental progenitors using **JLOH g2g**. These regions are going to be used to limit the output of **JLOH extract**. The `--est-divergence` parameter is the expected divergence between the two subgenomes (range: 0.0-1.0). The cross-mapping between the two parental genomes will be allowed up to to twice the indicated divergence.

```
jloh g2g --ref-A parent_A.fasta --ref-B parent_B.fasta --est-divergence 0.05 --min-length 1000 > A_and_B.regions.bed
```

Finally, we call LOH blocks using **JLOH extract**, limiting the output to regions showing divergence in genome-to-genome mapping.

```
jloh extract --threads 4 --vcfs A.ff.vcf B.ff.vcf --bams A.bam B.bam \
--refs parent_A.fasta parent_B.fasta --regions A_and_B.regions.bed \
--sample <STRING> --output-dir <PATH>
```

The `<sample>.LOH_blocks.{A,B}.tsv` files contained in `--output-dir` will contain all *bona fide* blocks found with this approach in either A or B genome. The `<sample>.LOH_candidates.{A,B}.tsv` files, instead, contain all the blocks that were called from the program, regardless of the regions specified in `--regions`.
