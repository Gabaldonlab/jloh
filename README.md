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

This will execute the command contained in the `*.sh` script, which in turn is a `nextflow run` command that runs the `*.nf` script with the `*.config` configuration file. In case you're working on a cluster with a slurm queuing system, you can edit the `#SBATCH` lines at the beginning and then run:

`sbatch reads_to_LOH_blocks.sh`

# Modules

## JLOH extract

This is the most important module of JLOH. Its functions are well described in the figure on the top. It is used to extract LOH blocks starting from VCF, BAM, and FASTA files.

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


### Sorting of SNPs by zygosity

The variants passed with `--vcf` are scanned, subdividing heterozygous and homozygous SNPs into two separate files: `<sample>.het_snps.vcf` and `<sample>.homo_snps.vcf`. Indels and other types of variation are discarded. The heterozygous SNPs are used to extract regions containing heterozygosity. Selected SNPs should also have an allele frequency (`AF`) annotation, and are retained if their `AF` is larger than `--min-af` and lower than `--max-af`. The default parameters (`--min-af 0.2 --max-af 0.8`) are probably ok for most users.

Missing allele frequency? Try using [all2vcf](https://github.com/MatteoSchiavinato/all2vcf).

If running in `--hybrid` mode, the homozygous SNPs are used to assign homozygous regions to either the alternative (ALT) or the reference (REF) allele. The selection of heterozygous SNPs is conducted based on their `FORMAT` (e.g. `GT 0/1` or `1/2` for heterozygous SNPs).

### Extraction of heterozygous regions

First, heterozygous SNPs are used to find heterozygous regions which are then masked as they cannot be LOH blocks. These are extracted based on clusters of heterozygous SNPs. The maximum SNP distance (`--snp-distance`) defines how far can the SNPs be while still being considered part of the same cluster, while the `--min-snps` parameter defines how many SNPs form a basic cluster.

The produced list of heterozygous regions that will then be ignored in the downstream LOH block detection.

### Extraction of candidate blocks

Everything that did not include sufficient heterozygous SNPs is then screened as a potential LOH block. Blocks shorter than `--min-length` are filtered out at this point. The remaining ones are screened against the initial heterozygous regions, trimming any overlapping region. The default parameters (`--min-snps 2 --snp-distance 100 --min-length 1000`) are probably ok for most users.

If running in `--hybrid` mode, clusters of **homozygous** SNPs are extracted the same way as for clusters of heterozygous SNPs. Blocks with sufficient SNPs in terms of `--min-snps` and `--snp-distance` will be considered as alternative allele blocks (i.e. `ALT`). Every region that is not a heterozygous SNP cluster or a homozygous SNP cluster is considered a reference allele homozygous block (`REF`).

### Coverage trimming

Candidate blocks are then screened in terms of coverage using the information in the `--bam` file. Blocks must be covered at least to a certain degree to be retained, which is controlled by `--min-frac-cov`.

Blocks that pass this filter are trimmed in their length according to coverage, removing their uncovered portions. Uncovered regions are used for this trimming step only if longer than `--min-uncov`. That is, if less than `--min-uncov` bases are uncovered within a block, the block stays untouched. In case this makes the block shorter than `--min-length`, the block is discarded. The default parameters (`--min-frac-cov 0.5 --min-uncov 10`) are probably ok for most users.

### Determination of block zygosity

The coverage information is then used to infer the zygosity of a block. The upstream and downstream regions of each block are extracted. The length of the up/downstream region to be considered is defined with the `--overhang` parameter. Some tolerance in the overhang length is accomodated by the `--min-overhang` parameter, defining a minimum fraction of `--overhang` to be considered sufficient. For example, if the user sets `--overhang 5000 --min-overhang 0.5`, a tolerance of 50% (2500 bp) is applied to the overhang length both up- and downstream of the block. If less than 2500 bp are available either up- or downstream, however, no zygosity is inferred. This is to avoid false positives at the beginning or end of a chromosome.

There are two types of zygosity predicted by JLOH: `homo` and `hemi`. These are defined by the `--hemi` parameter, which is a threshold from 0 to 1 applied to the coverage ratios calculated between the block and its flanking regions. The comparison with the upstream and the downstream regions returns two ratios. If both are below `--hemi`, the block is considered hemizygous. Otherwise, it is considered homozygous. If they show contrasting results, the block is labelled as "NA". This combines together with the presence or absence of homozygous SNPs in four possible scenarios (see image in the **Interpreting Output** chapter below). For example, if the user uses the default setting `--hemi 0.75`, if the block has less than 75% of the coverage of both its up- and downstream regions, it is labelled as hemizygous.

The default parameters (`--overhang 5000 --min-overhang 0.9 --hemi 0.75`) are probably ok for most users.

### Output

These are the five most important output files of JLOH:

- `<sample>.LOH_blocks.tsv`: the genomic ranges (1-based, closed) where LOH blocks have been found, together with other information such as coverage and number of SNPs. A BED file (0-based, half-open) is produced together with it.

- `<sample>.LOH_candidates.tsv`: the genomic ranges (1-based, closed) of **all** LOH blocks found, including those outside of `--regions`. A BED file (0-based, half-open) is produced together with it.

- `<sample>.het_blocks.bed`: the detected heterozygous regions which have been excluded from the LOH block extraction.

- `<sample>.exp.het_snps.vcf`: VCF file containing the heterozygous SNPs isolated from the input VCF.

- `<sample>.exp.homo_snps.vcf`: VCF file containing the homozygous SNPs isolated from the input VCF.


#### Interpreting output

JLOH's main output is a table containing all candidate LOH blocks. These blocks are annotated with various information on them, but the interpretation is strongly case-specific, hence is left to the user. We strongly advice the user to load the input VCF, BAM, and the output TSV to a genome browser such as [IGV](https://software.broadinstitute.org/software/igv/). The output TSV's first three columns are already a BED file (0-based, half-open intervals) and can be used as they are.

Ideally, LOH blocks assigned as `REF` should match regions that are covered by reads but depleted of *homozygous* SNPs, while those assigned as `ALT` should match regions dense in *homozygous* SNPs. Detected LOH blocks (regardless of the annotation) should not overlap regions that are dense in *heterozygous SNPs*.

In terms of zygosity, regions annotated as *homo* should have a read coverage that is higher or equal to the value set with the `--hemi` parameter (default: 0.75, i.e. 75%). Regions annotated as *hemi* should have a read coverage that is below this value.

An example as seen in IGV is provided below.  

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
