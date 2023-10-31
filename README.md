# J LOH
![Latest Version](https://img.shields.io/github/v/tag/gabaldonlab/jloh?label=Latest%20Version)
[![DOI](https://zenodo.org/badge/425015409.svg)](https://zenodo.org/badge/latestdoi/425015409)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.1-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://hub.docker.com/repository/docker/cgenomics/jloh)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

### What is JLOH?

JLOH is a tool to extract, filter, and analyse loss of heterozygosity (LOH) blocks based on single-nucleotide polymorphisms (SNPs), read mapping information, and a reference genome sequence.

*WARNING*: JLOH is made to assess LOH in genomes with at least 1% heterozygosity between their homolog chromosomes (or subgenomes, if hybrids). As of now, it isn't suitable for cancer data setups. 

As of September 2023, this is just a landing page. Read the full documentation on [jloh.readthedocs.io](http://jloh.readthedocs.io).

### What input does it need?

JLOH only needs three file types as input: 
- **VCF** file with called heterozygous and homozygous SNPs 
- **BAM/SAM** file with read mapping results, from which the variants were called
- **FASTA** file with a reference genome sequence where reads were mapped to get the BAM and VCF files

### I have a hybrid, does it work with it?

Yes, it does. The extraction module `jloh extract` can work with reads from a hybrid species mapped onto the hybrid's reference genome. See more about it on [jloh.readthedocs.io](http://jloh.readthedocs.io). 

