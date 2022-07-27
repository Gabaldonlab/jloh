# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract, filter, and manage blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome. Optionally works with a paired set of VCFs, BAMs and FASTAs produced from a hybrid genome (`--hybrid` mode).

![JLOH workflow](images/j_loh.png)

## Table of Contents

- [J LOH](#j-loh)
  * [Table of Contents](#table-of-contents)
  * [Install and run](#install-and-run)
  * [Basic usage](#basic-usage)
- [Modules detailed description](#modules-detailed-description)
  * [JLOH extract](#jloh-extract)
  * [JLOH filter](#jloh-filter)
  * [JLOH intersect](#jloh-intersect)
  * [JLOH chimeric](#jloh-chimeric)
  * [JLOH g2g](#jloh-g2g)
  * [JLOH density](#jloh-density)
  * [JLOH sim](#jloh-sim)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>


## Install and run

Detailed information [can be found here](docs/INSTALL.md).

## Basic usage

**jloh** has several modules that can be run independently. Many of them depend on the output of another module, making the analysis intercompatible. This is the list of modules, each of which is described more in detail below.

```
 -- Extraction
    extract             Extract LOH blocks from VCF, BAM and FASTA files
    filter              Filter extracted LOH blocks
    intersect           Perform intersection/removal operations with output files
    chimeric            Extract genes featuring LOH blocks from different haplotypes
    g2g                 Align two genomes to find regions that should carry SNPs

 -- Calculations
    density             Estimate heterozygous and homozygous SNP densities

 -- Simulations
    sim                 Simulate a divergent copy of a genome/protein sequence(s)
```

# Modules detailed description

## JLOH extract

This is the most important module of JLOH. Its functions are well described in the figure on the top. It is used to extract LOH blocks starting from VCF, BAM, and FASTA files. Detailed information on this algorithm functioning [can be found here](docs/EXTRACT.md).

**Output**

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

In `--hybrid` mode, it produces the same output files but repeated twice: once per each parent. More information on the `--hybrid` mode, and how to use it properly, [can be found here](docs/HYBRID.md).

## JLOH filter

This tool filters the output produced by `JLOH extract` according to several criteria. The user can select to filter the LOH blocks based on their coverage, on their SNP count, SNP density, length, or extract individual regions just like in samtools.

## JLOH intersect

This tool perform ensemble operations with two JLOH output files, namely intersection, complement, and unique elements extraction.

## JLOH chimeric

This tool is a specific module to extract genes that overlap LOH blocks from two different origins (i.e. chimeric genes). It does not assume anything about the lists of LOH blocks that are passed. If a gene has two LOH blocks in its sequence (one from each list) it will be considered a candidate chimeric gene.

The usage involves two sets of LOH blocks produced by `jloh extract`, plus the `*het_blocks.bed` file also produced by `jloh extract`. In case the "extract" module was run in `--hybrid` mode, the heterozygous blocks files are two (A and B). In that case, they must be concatenated into a single file using `cat` and provided as a single file to `jloh chimeric`.

**Output**

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
