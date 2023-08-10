# J LOH
![Latest Version](https://img.shields.io/github/v/tag/gabaldonlab/jloh?label=Latest%20Version)
[![DOI](https://zenodo.org/badge/425015409.svg)](https://zenodo.org/badge/latestdoi/425015409)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/repository/docker/cgenomics/jloh)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

### What is JLOH?

JLOH is a tool to extract, filter, and analyse loss of heterozygosity (LOH) blocks based on single-nucleotide polymorphisms (SNPs), read mapping information, and a reference genome sequence.

### What input does it need?

JLOH only needs three file types as input: 
- **VCF** file with called heterozygous and homozygous SNPs 
- **BAM/SAM** file with read mapping results, from which the variants were called
- **FASTA** file with a reference genome sequence where reads were mapped to get the BAM and VCF files

### I have a hybrid, does it work with it?

Yes, it does. The extraction module `jloh extract` can work with reads from a hybrid species mapped onto the hybrid's reference genome. It will produce LOH blocks although you won't know to which subgenome they belong to. If you have the parental genomes of the hybrid, however, `jloh extract` has an `--assign-blocks` option which allows to use them to assign LOH blocks to subgenomes. In that case, follow the instructions in [this guide](docs/ASSIGN_BLOCKS.md) to provide two reference genomes, two BAM/SAM files containing mapping results, and two VCF files containing SNPs. Each pair corresponds to the two parentals. See details [here](#jloh-extract).

![JLOH workflow](images/j_loh.png)

Schematic of the "jloh extract" module workflow.

## Table of Contents

- [J LOH](#j-loh)
  * [Install and run](#install-and-run)
  * [Basic usage](#basic-usage)
- [Modules detailed description](#modules-detailed-description)
  * [Extraction tools](#extraction-tools)
    + [JLOH extract](#jloh-extract)
    + [JLOH filter](#jloh-filter)
    + [JLOH intersect](#jloh-intersect)
    + [JLOH chimeric](#jloh-chimeric)
    + [JLOH g2g](#jloh-g2g)
  * [Calculation tools](#calculation-tools)
    + [JLOH stats](#jloh-stats)
    + [JLOH junctions](#jloh-junctions)
  * [Simulation tools](#simulation-tools)
    + [JLOH sim](#jloh-sim)
  * [Visualization tools](#visualization-tools)
    + [JLOH plot](#jloh-plot)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Install and run

Detailed information [can be found here](docs/INSTALL.md).

## Basic usage

**jloh** has several modules that can be run independently. Many of them depend on the output of another module, making the analysis intercompatible. This is the list of modules, each of which is described more in detail below.

```
 -- Extraction
    g2g                 Align two genomes to find regions that should carry SNPs
    extract             Extract LOH blocks from VCF, BAM and FASTA files
    filter              Filter extracted LOH blocks
    intersect           Perform intersection/removal operations with output files
    chimeric            Extract genes featuring LOH blocks from different haplotypes

 -- Calculations
    stats               Estimate heterozygous and homozygous SNP statistics
    junctions           Calculate number of block-to-block junctions over the genome

 -- Simulation
    sim                 Simulate a divergent copy of a genome/protein sequence(s)

 -- Visualization
    plot                Make an LOH propensity plot from "extract" output file(s)

```

# Modules detailed description

## Extraction tools 

### JLOH extract

This is the most important module of JLOH. Its functions are well described in the figure on the top. It is used to extract LOH blocks starting from VCF, BAM, and FASTA file(s). Detailed information on how this algorithm works [can be found here](docs/EXTRACT.md).

This module has three ways of operating: 

| Parameter setting    | Input reference(s)     |
|----------------------|------------------------|
| default (no setting) | same organism as reads | 
| `--assign-blocks`    | parentals of a hybrid  |
| `--one-parent`       | one parental only      |

**Output**

`jloh extract` functions in three modes (default, `--assign-blocks` and `--one-parent`). 

In default mode, it produces these output files:

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

In `--assign-blocks` mode, it produces the same output files but repeated twice: once per each parent. More information on the `--assign-blocks` mode, and how to use it properly, [can be found here](docs/ASSIGN_BLOCKS.md).

In `--one-parent` mode, it produces the same output of the default mode. However, it will imply that the provided reference genome is from one of the two parental progenitors, and use it to annotate blocks towards the REF or the ALT allele as in the `--assign-blocks` mode. The difference from the latter is that only one parental genome is provided. The difference from the default mode is that the provided genome is not from the hybrid itself but from one of the parents. 

### JLOH filter

This tool filters the output produced by `JLOH extract` according to several criteria. The user can select to filter the LOH blocks based on their coverage, on their SNP count, SNP density, length, or extract individual regions just like in samtools.

### JLOH intersect

This tool perform ensemble operations with two JLOH output files, namely intersection, complement, and unique elements extraction.

### JLOH chimeric

This tool is a specific module to extract genes that overlap LOH blocks from two different origins (i.e. chimeric genes). It does not assume anything about the lists of LOH blocks that are passed. If a gene has two LOH blocks in its sequence (one from each list) it will be considered a candidate chimeric gene.

The usage involves two sets of LOH blocks produced by `jloh extract`, plus the `*het_blocks.bed` file also produced by `jloh extract`. In case the "extract" module was run in `--assign-blocks` mode, the heterozygous blocks files are two (A and B). In that case, they must be concatenated into a single file using `cat` and provided as a single file to `jloh chimeric`.

**Output**

`jloh chimeric` produces two files:

| Output file               | Description |
|---------------------------|-|
| out.chimeric.features.gff | GFF file containing rows that identify genes marked as chimeric. This file can be directly loaded into IGV or GBrowse for manual inspection, making your job easier. |
| out.chimeric.IDs.txt      | TXT file containing the gene IDs of all the genes found as chimeric. These can be of three types: 1) A-B chimeras, from a mixture of the two LOH block sources; 2) A-H chimeras, from a mixture of the A source and the heterozygous blocks, 3) B-H chimeras, from a mixture of the B source and the heterozygous blocks. |

### JLOH g2g

This program finds diverging regions between two genomes in FASTA format. The input are the two sequences and an estimated divergence value, and the output is a bed file representing the regions that contain SNPs between the two genomes. This module is useful only when using the `--assign-blocks` mode in `jloh extract`, as the output BED file of `g2g` can be directly passed to the `--regions` parameter.

`JLOH g2g` runs more than one tool from the MUMmer arsenal to map the two genomes, filter the results, extract the SNPs. Then, it uses `all2vcf mummer` to convert the MUMmer output to VCF format, and `bedtools merge` to generate BED intervals from SNPs. Intervals are expanded as long as there are overlaps.

## Calculation tools 

### JLOH stats

This tool computes the densities of all SNPs, heterozygous SNPs, and homozygous SNPs over the genome sequence.

Besides the mean snp density, it calculates a distribution and extracts the most relevant quantiles.

These quantiles are useful to choose what value to set in the `--min-snps-kbp` parameter of `jloh extract`. 

A detailed description on how to interpret these quantiles [can be found here](docs/QUANTILES.md).

### JLOH junctions 

This tool is used to extract intervals within a genome sequence where blocks from two different sources are found. These are the candidate junction regions between different haplotypes.

This module is similar to `jloh chimeric`, with the difference that it's focus is not genes, but rather genome regions. 

## Simulation tools 

### JLOH sim

This module generates a copy of a reference sequence in FASTA format, introducing a series of mutations selected randomly over the sequence of each scaffold/chromosome. Optionally, the module can include a series of LOH blocks defined by percentage of the whole genome (e.g. 20%).

The output is:
- the mutated reference sequence in FASTA format
- a tab-separated file with all the introduced SNPs, with positions annotated in 1-based format including reference and alternative allele
- a BED file with all the introduced LOH blocks, hence in 0-based half-open format
- another BED file with all the regions that are deemed as non-divergent based on the declared `--snp-distance`.

## Visualization tools 

### JLOH plot 

This module is used to represent graphically the results obtained from `jloh extract`. We suggest to use it with multiple input data, as it will produce an individual plot for each chromosome with all the samples included within it. 

The resulting plot is a heatmap which is easy to interpret and highlights the LOH tendency of each region of each chromome, be it the reference allele (REF) or the alternative allele (ALT). 
