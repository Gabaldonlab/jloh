# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome.

![JLOH workflow](images/j_loh.png)

## Install

As simple as: `git clone https://github.com/Gabaldonlab/jloh.git` or `https://github.com/MatteoSchiavinato/jloh.git`
And it's ready to go! But there are a few dependencies:

| Program     | Type        | Version | Links      |
|-------------|-------------|---------|------------|
| pybedtools  | Module      | 0.8.2   | [source](https://daler.github.io/pybedtools/main.html), [cite](https://doi.org/10.1093/bioinformatics/btr539) |
| pysam       | Module      | 0.1.7   | [source](https://pypi.org/project/pysam/), [cite](https://github.com/pysam-developers/pysam) |
| Python      | Interpreter | 3.6.1   | [source](https://www.python.org/downloads/release/python-397/), [cite](http://citebay.com/how-to-cite/python/) |

Note that **pybedtools** will look for [bedtools](https://bedtools.readthedocs.io/en/latest/) in the `$PATH`, while **pysam** will look for [samtools](http://www.htslib.org/). Hence, you may have to install these two programs as well.  

## Run

The basic usage of the program is as simple as:

```
./jloh --vcf <VCF> --genome-file <GENOME_FILE> --bam <BAM> [options]
```

To produce a genome file, simply calculate the length of each sequence in your reference FASTA file and produce a file containing their name + length, structured in a tab-separated format that looks like this:

```
chr1  89214414
chr2  1231455
chr3  90804782
...
```

## The workflow, described

This section describes in detail what you can achieve with **JLOH**.

### Command-line arguments

JLOH has many command-line arguments that can be set to fine-tune the analysis:

```
[I/O]
--vcf               Input VCF file                                              [!]
--genome-file       File with chromosome lengths (chromosome <TAB> size)        [!]
--bam               BAM file used to call the --vcf variants                    [!]
--sample            Sample name / Strain name for output files                  [loh_blocks]
--output-dir        Output directory                                            [loh_blocks_out]
--debug             Activate generation of several intermediate files           [off]
--print-info        Show authors and edits with dates                           [off]

[modes]
--keep-secondary    Ignores secondary alignments from --bam file                [off]
--no-alleles        Don't use homozygous SNPs to assign LOH blocks to REF/ALT   [off]

[parameters]
--filter-mode       "pass" to keep only PASS variants, "all" to keep everything [all]
--min-het-snps      Min. num. heterozygous SNPs in heterozygous region          [2]
--snp-distance      Max. distance (bp between SNPs for blocks definition        [100]
--block-dist        Combine LOH blocks into one if closer than this distance    [100]
--min-size          Min. LOH block size                                         [100]
--min-af            Min. allele frequency to consider a variant heterozygous    [0.3]
--max-af            Max. allele frequency to consider a variant heterozygous    [0.7]
--min-frac-cov      Min. fraction of LOH block that has to be covered by reads  [0.5]
--hemi              Frac. of the mean coverage under which LOH is hemizygous    [0.75]

[pre-existing variation]
--t0-vcf            VCF with variants to ignore from --vcf                      [off]
--t0-bam            BAM file used to call the --t0-vcf variants                 [off]
--t0-filter-type    What to do with t0 LOH events? "keep" or "remove"           [remove]
```

### selection of SNPs

The variants passed with `--vcf` are filtered, retaining only heterozygous SNPs which are placed in an output file called `<sample>.het_snps.vcf`. Homozygous SNPs are placed in another output file called `<sample>.homo_snps.vcf`. Indels and other types of variation are discarded.

The heterozygous SNPs are used to extract regions containing heterozygosity, while the homozygous SNPs are used to assign homozygous regions to either the alternative (ALT) or the reference (REF) allele; this is done by default, unless the user passes the `--no-alleles` option (see later).

The selection of heterozygous SNPs is conducted based on their `FORMAT` (field number 9 and 10 of a VCF file). The first column of this field carries a series of annotations separated by colons (e.g. GT:AF) the values of which are annotated the same way on the second column (e.g.` 0/1:0.60`). If a SNP is annotated as heterozygous, it will carry a genotype (`GT`) such as `0/1` or `1/2`. It should also have an allele frequency (`AF`) annotation. SNPs are considered heterozygous if their `AF` annotation falls between the values specified with `--min-af` and `--max-af`.

Missing allele frequency? Try using [all2vcf](https://github.com/MatteoSchiavinato/all2vcf).

### extraction of heterozygosity regions

The first block of operations is aimed at the extraction of genomic regions rich in *heterozygous SNPs*. These regions cannot be LOH blocks and are identified in order to be masked. The heterozygous regions are extracted based on clusters of heterozygous SNPs. The minimum number of SNPs eliminates regions with too little SNPs to be considered true positives, while the maximum SNP distance defines how the SNPs are clustered. These parameters are controlled with the `--min-het-snps` and the `--snp-distance` parameters. This produces a list of genomic intervals defining the heterozygous clusters, including the number of SNPs contained in them. These intervals will be masked for later LOH detection.

### Extraction of candidate blocks

Everything that did not include sufficient heterozygous SNPs is then screened as a potential LOH block. First, the complementary intervals of the heterozygous ones are extracted. Then, the program makes a choice depending on the `--no-alleles` settings:

- If not set, the same procedure used for the heterozygous SNPs is repeated to extract clusters of **homozygous** SNPs. These indicate that the candidate LOH block has retained an alternative allele (i.e. `ALT`), different from the one in the reference genome. Regions without homozygous SNPs confirm, instead, the reference allele (i.e. `REF`). The identified intervals are then used to subdivide candidate blocks into **ALT** and **REF**.

- If set, JLOH does not perform the ALT/REF analysis and assigns an arbitrary `NA` allele annotation to all candidate blocks. Depending on the species you work with (hybrid, poorly assembled, etc.) this may be a safe choice.

Regardless of the choice, at the end of this block of operations JLOH has a list of potential LOH blocks i terms of chromosome, start, end, number of SNPs, and length. Blocks shorter than `--min-size` are at this point filtered out. The remaining ones are screened against the initial heterozygous regions: overlapping regions are cropped from the candidate LOH block coordinates.

### Determination of block zygosity

The third block of operations involves the determination of whether each candidate LOH block is homozygous or hemizygous. This is done through read coverage: homozygous blocks will have a read coverage that is comparable with that of the global mean coverage per position; hemizygous blocks will have lost one of the two copies (if diploid) and will have only half of the global mean coverage. This is annotated in the output file and can be controlled with the `--hemi` parameter.

#### How coverage is used

Each region that is considered as a candidate LOH region is screened by coverage using the BAM file passed with `--bam`. First, the global mean coverage per position is computed for the whole BAM file. To do that, JLOH checks if the BAM file is indexed, and if not, it indexes it using the **pysam** module. Then, the reads mapping inside each region are extracted using the **pysam** module, and compared against the global mean coverage.

Each block must also pass a *covered fraction* filter (the fraction of positions actually covered by reads). If the fraction is lower than `--min-frac-cov`, the block is discarded.

#### The `--keep-secondary` option

In the default scenario, the user should use only primary alignments for their analysis. This is strongly recommended and is commonly done, to avoid coverage overestimations determined by multi-mapping reads with many secondary alignments. If the user sets the `--keep-secondary` option, all provided reads are used. This may be useful to assess impact of secondary alignments on your analysis, but is not recommended unless you know what you're doing.

### deal with known pre-existing LOH blocks

The last block of operations is aimed at comparing the detected LOH blocks with any pre-existing ones. This is an *optional* step and is not performed unless the user passes another VCF file with the `-t0-vcf` option. This file will be considered as variation that pre-dates the one listed in `--vcf`.

This "t0" variation is used to extract LOH blocks the same way as described for the input VCF file. The user can then choose what to do with it with the `--t0-filter-type` option. The default ("remove") will remove any overlapping LOH block found between the "t0" VCF and the input VCF. This reduces the output LOH blocks only to those unique to the input VCF. Some users however may want to keep only blocks that are found in "t0" too. To do that, one just has to specity `--t0-filter-type keep`.

By default, JLOH won't expect a "t0" dataset and will go through this step without doing anything.

## Output

**JLOH** produces the following output files:

- `<sample>.LOH_blocks.tsv`: the main output of the file. The first three columns represent the genomic ranges where LOH blocks have been found. These can be easily cut to produce a BED file. It also contains: the relative coverage with respect to the genome mean coverage (`#4`), the length of the block (`#5`), the allele to which it has been assigned (`#6`), the number of homozygous SNPs found within it (`#7`) and the number of heterozygous SNPs found within it (`#8`).

- `<sample>.exp.LOH_candidates.tsv`: this file contains candidate LOH blocks that have been discarded due to any of the reasons why the script discards a block (coverage, length, number of SNPs contained). These blocks are placed in this file for debugging reasons and for the user to be aware of what is not kept.

- `<sample>.exp.het_snps.vcf`: VCF file containing the heterozygous SNPs isolated from the input VCF.

- `<sample>.exp.homo_snps.vcf`: VCF file containing the homozygous SNPs isolated from the input VCF.

Note: if a `--t0-vcf` and a `-t0-bam` file are passed, there will be extra files with the same naming system corresponding to these input files. These will carry the label "t0" instead of "exp".

## Interpreting output

JLOH's main output is a table containing all candidate LOH blocks. These blocks are annotated with various information on them, but the interpretation is strongly case-specific, hence is left to the user. We strongly advice the user to load the input VCF, BAM, and the output TSV to a genome browser such as [IGV](https://software.broadinstitute.org/software/igv/). The output TSV's first three columns are already a BED file (0-based, half-open intervals) and can be used as they are.

Ideally, LOH blocks assigned as `REF` should match regions that are depleted of *homozygous* SNPs, while those assigned as `ALT` should match regions dense in *homozygous* SNPs. Detected LOH blocks (regardless of the annotation) should not overlap regions that are dense in *heterozygous SNPs*.

In terms of zygosity, regions annotated as *homo* should have a read coverage that is higher or equal to the value set with the `--hemi` parameter (default: 0.75, i.e. 75%). Regions annotated as *hemi* should have a read coverage that is below this value.

An example as seen in IGV is provided below.  

![Example](images/example.png)
