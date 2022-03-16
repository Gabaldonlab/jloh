# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract, filter, and manage blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome.

![JLOH workflow](images/j_loh.png)

## Table of Contents

- [J LOH](#j-loh)
  * [Table of Contents](#table-of-contents)
  * [Install](#install)
  * [Run](#run)
    + [Nextflow workflow](#nextflow-workflow)
  * [Example run - step by step guide](#example-run---step-by-step-guide)
- [Implementation](#implementation)
  * [JLOH extract](#jloh-extract)
    + [Sorting of SNPs by zygosity](#sorting-of-snps-by-zygosity)
    + [Extraction of heterozygous regions](#extraction-of-heterozygous-regions)
    + [Extraction of candidate blocks](#extraction-of-candidate-blocks)
    + [Coverage trimming](#coverage-trimming)
    + [Determination of block zygosity](#determination-of-block-zygosity)
    + [Output](#output)
      - [Interpreting output](#interpreting-output)
- [Modules](#modules)
  * [JLOH sim](#jloh-sim)
  * [JLOH g2g](#jloh-g2g)
  * [JLOH extract](#jloh-extract-1)
  * [JLOH filter](#jloh-filter)
  * [JLOH density](#jloh-density)


## Install

As simple as: `git clone https://github.com/Gabaldonlab/jloh.git`
And it's ready to go! But there are a few dependencies:

| Program     | Type        | Version | Links      |
|-------------|-------------|---------|------------|
| all2vcf     | Program     | 0.7.3   | [source](https://github.com/MatteoSchiavinato/all2vcf), [cite](https://github.com/MatteoSchiavinato/all2vcf) |
| bedtools    | Program     | 2.30    | [source](https://bedtools.readthedocs.io/en/latest/), [cite](https://doi.org/10.1002/0471250953.bi1112s47) |
| Biopython   | Module      | 1.79    | [source](https://biopython.org/), [cite](https://doi.org/10.1093/bioinformatics/btp163) |
| MUMmer      | Program     | 3.1     | [source](https://anaconda.org/bioconda/mummer), [cite](https://doi.org/10.1186%2Fgb-2004-5-2-r12) |
| numpy       | Module      | 1.21.4  | [source](https://numpy.org/), [cite](https://doi.org/10.1038/s41586-020-2649-2) |
| pandas      | Module      | 1.3.5   | [source](https://pandas.pydata.org/), [cite](https://doi.org/10.5281/zenodo.3509134) |
| pybedtools  | Module      | 0.8.2   | [source](https://daler.github.io/pybedtools/main.html), [cite](https://doi.org/10.1093/bioinformatics/btr539) |
| pysam       | Module      | 0.1.7   | [source](https://pypi.org/project/pysam/), [cite](https://github.com/pysam-developers/pysam) |
| Python      | Interpreter | 3.6.1   | [source](https://www.python.org/downloads/release/python-397/), [cite](http://citebay.com/how-to-cite/python/) |
| samtools    |  Program    | 1.13    | [source](http://www.htslib.org/), [cite](https://doi.org/10.1093/gigascience/giab008) |

Note that **pybedtools** will look for [bedtools](https://bedtools.readthedocs.io/en/latest/) in the `$PATH`, while **pysam** will look for [samtools](http://www.htslib.org/).

The installation of **MUMmer** can be easily done via conda (`conda install -c bioconda mummer`). This will place in the `$PATH` all the toolkit from the MUMmer arsenal, in particular the tools needed by JLOH to run: nucmer, delta-filter, show-snps.

The installation of **all2vcf** is very straightforward. First clone the repository:

```
git clone https://github.com/MatteoSchiavinato/all2vcf
```

Then make a symbolic link of the main `all2vcf` executable (not of the whole folder!) into your `/bin`:

```
ln -s /path/to/all2vcf /path/to/bin
```

## Run

JLOH has different tools which are used to extract and perform different operations on LOH blocks. To see all available tools simply run:

```
./jloh -h
```

### Nextflow workflow

Together with the **JLOH** tool, we provide also a [Nextflow](http://nextflow.io/) workflow that you can use to run your samples directly from raw reads to LOH blocks. All you have to do is to edit the configuration file of the workflow (\*config) and the running script (\*sh). Edit them according to your own computer / server, and then run:

`bash reads_to_LOH_blocks.sh`

In case you're working on a cluster with a slurm queuing system, you can edit the `#SBATCH` lines at the beginning and then run:

`sbatch reads_to_LOH_blocks.sh`


## Example run - step by step guide

The following steps (apart from the installation) are performed by the nextflow workflow described in the previous section. Maybe have a look at that before doing everything by hand! ;)

The first step is of course making sure that everything has been installed properly. Hence, follow the instructions in the [Install](#install) section prior to trying this guide.

Once everything is installed properly, we start by mapping a set of quality-trimmed paired-end reads sequenced from a hybrid species against both of its parental genome sequences. We do it with **HISAT2** but you can use the tool you prefer. For the whole paragraph we will use four threads, represented here as `-p 4`.

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

Note that we don't care about indels in this workflow, using the `--skip-indels` option, and that we care instead to annotate a lot of INFO/ and FORMAT/ fields of the generated VCF file.

We then call the SNPs using **bcftools** again.

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

Now, we extract the regions of high sequence identity between the two hybrid genomes used in the analysis with **JLOH g2g**. These regions are going to be masked from the produced LOH blocks later on due to the fact that we can't be sure they are true LOH blocks and not just regions that never were heterozygous in the first place. The `--min-identity` parameters should reflect the expected divergence between the two subgenomes, e.g. with 5% divergence you can set it to 95.

```
jloh g2g --ref-A parent_A.fasta --ref-B parent_B.fasta --min-identity 95 --min-length 1000 > A_and_B.mask.bed
```

Finally, we call LOH blocks using **JLOH extract**.

```
jloh extract --threads 4 --vcfs A.ff.vcf B.ff.vcf --bams A.bam B.bam \
--refs parent_A.fasta parent_B.fasta --mask A_and_B.mask.bed \
--sample <STRING> --output-dir <PATH>
```

The `<sample>.LOH_blocks.tsv` file contained in `--output-dir` will contain all *bona fide* blocks found with this approach, and is the output of the workflow. **CAREFUL**: these blocks are annotated in 0-based coordinates (e.g. positions from 1 to 20 are annotated as 0:20). A bed version is provided together with it.

# Implementation

This section describes in detail what you can achieve with **JLOH**. Later in this guide, a nextflow workflow to go from reads to LOH blocks is also described.

## JLOH extract

This tool extracts LOH blocks from a VCF file, a BAM file, and their reference genome in FASTA format.

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

### Output

These are the important output files of JLOH:  

- `<sample>.LOH_blocks.tsv`: the main output of the file. The first three columns represent the genomic ranges where LOH blocks have been found. These can be easily cut to produce a BED file. It also contains: the relative coverage with respect to the genome mean coverage (`#4`), the length of the block (`#5`), the allele to which it has been assigned (`#6`), the number of homozygous SNPs found within it (`#7`) and the number of heterozygous SNPs found within it (`#8`).

- `<sample>.exp.het_snps.vcf`: VCF file containing the heterozygous SNPs isolated from the input VCF.

- `<sample>.exp.homo_snps.vcf`: VCF file containing the homozygous SNPs isolated from the input VCF.

#### Interpreting output

JLOH's main output is a table containing all candidate LOH blocks. These blocks are annotated with various information on them, but the interpretation is strongly case-specific, hence is left to the user. We strongly advice the user to load the input VCF, BAM, and the output TSV to a genome browser such as [IGV](https://software.broadinstitute.org/software/igv/). The output TSV's first three columns are already a BED file (0-based, half-open intervals) and can be used as they are.

Ideally, LOH blocks assigned as `REF` should match regions that are covered by reads but depleted of *homozygous* SNPs, while those assigned as `ALT` should match regions dense in *homozygous* SNPs. Detected LOH blocks (regardless of the annotation) should not overlap regions that are dense in *heterozygous SNPs*.

In terms of zygosity, regions annotated as *homo* should have a read coverage that is higher or equal to the value set with the `--hemi` parameter (default: 0.75, i.e. 75%). Regions annotated as *hemi* should have a read coverage that is below this value.

An example as seen in IGV is provided below.  

![Example](images/example.png)

# Modules

## JLOH sim

This module generates a copy of a reference sequence in FASTA format, introducing a series of mutations selected randomly over the sequence of each scaffold/chromosome. Optionally, the module can include a series of LOH blocks defined by percentage of the whole genome (e.g. 20%).

The output is:
- the mutated reference sequence in FASTA format
- a tab-separated file with all the introduced SNPs, with positions annotated in 1-based format including reference and alternative allele
- a BED file with all the introduced LOH blocks, hence in 0-based half-open format

## JLOH g2g

This program is made for extracting regions of high sequence identity between two genome sequences in FASTA format. The input are the two sequences, and the output is a bed file representing the regions that are depleted of SNPs between the two genomes.

`JLOH g2g` runs more than one tool from the MUMmer arsenal to map the two genomes, filter the results, extract the SNPs. Then, it uses `all2vcf mummer` to convert the MUMmer output to VCF format, and `bedtools merge` to generate BED intervals from SNPs. Intervals are expanded as long as there are overlaps, and at the end are reversed, to find regions without SNPs.

These regions are a good `--mask` to pass to `JLOH extract` in `--hybrid` mode.

## JLOH extract

This is the most important module of JLOH. Its functions are well described above. It is used to extract LOH blocks starting from VCF, BAM, and FASTA files.

## JLOH filter

This tool filters the output produced by `JLOH extract` according to several criteria. The user can select to filter the LOH blocks based on their coverage, on their SNP count, SNP density, length, or extract individual regions just like in samtools.

## JLOH density

This tool computes the densities of all SNPs, heterozygous SNPs, and homozygous SNPs over the genome sequence. It is meant to be an informative tool for the user.
