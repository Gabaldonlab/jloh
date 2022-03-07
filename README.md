  # J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract, filter, and manage blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome.

![JLOH workflow](images/j_loh.png)

## Table of Contents

- [Install](#Install)
- [Run](#Run)
- [Implementation](#Implementation)
- [Output](#Output)

## Install

As simple as: `git clone https://github.com/Gabaldonlab/jloh.git`
And it's ready to go! But there are a few dependencies:

| Program     | Type        | Version | Links      |
|-------------|-------------|---------|------------|
| Biopython   | Module      | 1.79    | [source](https://biopython.org/), [cite](https://doi.org/10.1093/bioinformatics/btp163) |
| pybedtools  | Module      | 0.8.2   | [source](https://daler.github.io/pybedtools/main.html), [cite](https://doi.org/10.1093/bioinformatics/btr539) |
| pysam       | Module      | 0.1.7   | [source](https://pypi.org/project/pysam/), [cite](https://github.com/pysam-developers/pysam) |
| Python      | Interpreter | 3.6.1   | [source](https://www.python.org/downloads/release/python-397/), [cite](http://citebay.com/how-to-cite/python/) |

Note that **pybedtools** will look for [bedtools](https://bedtools.readthedocs.io/en/latest/) in the `$PATH`, while **pysam** will look for [samtools](http://www.htslib.org/). Hence, you may have to install these two programs as well.  

## Run

JLOH has different tools which are used to extract and perform different operations on LOH blocks. The main tool is `JLOH extract`: its basic usage is as simple as:

```
./jloh extract --threads <num_threads> --vcf <VCF> --ref <FASTA> --bam <BAM> [options]
```

To see all available tools simply run:

```
./jloh -h
```

### Nextflow workflow

Together with the **JLOH** tool, we provide also a [Nextflow](http://nextflow.io/) workflow that you can use to run your samples directly from raw reads to LOH blocks. All you have to do is to edit the configuration file of the workflow (\*config) and the running script (\*sh). Edit them according to your own computer / server, and then run:

`bash reads_to_LOH_blocks.sh`

In case you're working on a cluster with a slurm queuing system, you can edit the `#SBATCH` lines at the beginning and then run:

`sbatch reads_to_LOH_blocks.sh`


## Implementation

This section describes in detail what you can achieve with **JLOH**. Later in this guide, a nextflow workflow to go from reads to LOH blocks is also described.

### Sorting of SNPs by zygosity

The variants passed with `--vcf` are scanned, subdividing heterozygous and homozygous SNPs into two separate files: `<sample>.het_snps.vcf` and `<sample>.homo_snps.vcf`. Indels and other types of variation are discarded. The heterozygous SNPs are used to extract regions containing heterozygosity, while the homozygous SNPs are used to assign homozygous regions to either the alternative (ALT) or the reference (REF) allele. The selection of heterozygous SNPs is conducted based on their `FORMAT` (e.g. `GT 0/1` or `1/2` for heterozygous SNPs). Selected SNPs should also have an allele frequency (`AF`) annotation, and are retained if their `AF` is larger than `--min-af` and lower than `--max-af`. The default parameters (`--min-af 0.2 --max-af 0.8`) are probably ok for most users.

Missing allele frequency? Try using [all2vcf](https://github.com/MatteoSchiavinato/all2vcf).

### Extraction of heterozygous regions

Heterozygous regions are extracted based on clusters of heterozygous SNPs. The maximum SNP distance (`--snp-distance`) defines how far can the SNPs be while still being considered part of the same cluster, while the `--min-snps` parameter defines how many SNPs must be found within a region to consider it. This produces a list of heterozygous regions that will then be ignored in LOH block detection.

### Extraction of candidate blocks

Everything that did not include sufficient heterozygous SNPs is then screened as a potential LOH block. Clusters of **homozygous** SNPs are extracted the same way as for clusters of heterozygous SNPs. Blocks with sufficient SNPs in terms of `--min-snps` and `--snp-distance` will be considered as alternative allele blocks (i.e. `ALT`). Every region that is not a heterozygous SNP cluster or a homozygous SNP cluster is considered a reference allele homozygous block (`REF`).

At the end of this block of operations JLOH has a list of potential LOH blocks in terms of chromosome, start, end, number of SNPs, and length. Blocks shorter than `--min-length` are filtered out at this point. The remaining ones are screened against the initial heterozygous regions, trimming any overlapping region. The default parameters (`--min-snps 2 --snp-distance 100 --min-length 1000`) are probably ok for most users.

### Coverage trimming

Candidate blocks are screened against the coverage profiles of each chromosome using the BAM file passed with `--bam`. Regions of candidate blocks that don't have any coverage are trimmed, reducing the candidate block coordinates only to the covered portion of it. In case this makes it too short, i.e. shorter than `--min-length`, the block is consequently discarded.

Each block must also pass a *covered fraction* filter (the fraction of positions actually covered by reads). If the fraction is lower than `--min-frac-cov`, the block is discarded. Uncovered regions are used for this trimming step only if longer than `--min-uncov`. That is, if less than `--min-uncov` bases are uncovered within a block, the block stays untouched. The default parameters (`--min-frac-cov 0.5 --min-uncov 10`) are probably ok for most users.

### Determination of block zygosity

The coverage information is then used to infer the zygosity of a block. The upstream and downstream regions of each block are extracted. The extent of the up/downstream region is regulated by the `--overhang` and the `--min-overhang` parameters. The first defines how many bp have to be considered for this coverage test, while the second one defines the minimum fraction of them to assign zygosity. For example, if the user sets `--overhang 5000 --min-overhang 0.5`, a tolerance of 50% (2500 bp) is applied to the overhang length both up- and downstream of the block. If less than 2500 bp are available either up- or downstream, however, no zygosity is inferred. This is to avoid false positives at the beginning or end of a chromosome.

There are two types of zygosity predicted by JLOH: `homo` and `hemi`. These are defined by the `--hemi` parameter, which is a threshold from 0 to 1 applied to the coverage ratios calculated between the block and its flanking regions. The comparison with the upstream and the downstream regions returns two ratios (range: 0-1). If both are below `--hemi`, the block is considered hemizygous. Otherwise, it is considered homozygous. If they show contrasting results, the block is labelled as "NA". This combines together with the presence or absence of homozygous SNPs in four possible scenarios (see image in the **Interpreting Output** chapter below). The default parameters (`--overhang 5000 --min-overhang 0.9 --hemi 0.75`) are probably ok for most users.

## Output

These are the important output files of JLOH:  

- `<sample>.LOH_blocks.tsv`: the main output of the file. The first three columns represent the genomic ranges where LOH blocks have been found. These can be easily cut to produce a BED file. It also contains: the relative coverage with respect to the genome mean coverage (`#4`), the length of the block (`#5`), the allele to which it has been assigned (`#6`), the number of homozygous SNPs found within it (`#7`) and the number of heterozygous SNPs found within it (`#8`).

- `<sample>.exp.het_snps.vcf`: VCF file containing the heterozygous SNPs isolated from the input VCF.

- `<sample>.exp.homo_snps.vcf`: VCF file containing the homozygous SNPs isolated from the input VCF.

### Interpreting output

JLOH's main output is a table containing all candidate LOH blocks. These blocks are annotated with various information on them, but the interpretation is strongly case-specific, hence is left to the user. We strongly advice the user to load the input VCF, BAM, and the output TSV to a genome browser such as [IGV](https://software.broadinstitute.org/software/igv/). The output TSV's first three columns are already a BED file (0-based, half-open intervals) and can be used as they are.

Ideally, LOH blocks assigned as `REF` should match regions that are covered by reads but depleted of *homozygous* SNPs, while those assigned as `ALT` should match regions dense in *homozygous* SNPs. Detected LOH blocks (regardless of the annotation) should not overlap regions that are dense in *heterozygous SNPs*.

In terms of zygosity, regions annotated as *homo* should have a read coverage that is higher or equal to the value set with the `--hemi` parameter (default: 0.75, i.e. 75%). Regions annotated as *hemi* should have a read coverage that is below this value.

An example as seen in IGV is provided below.  

![Example](images/example.png)
