# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs), read mapping, and a reference genome. The main assumption behind this tool is that a true LOH event will be found where two conditions are verified in a specific DNA region:

- No heterozygous SNPs are called
- The region mean coverage increases (or decreases) by at least a certain amount when compared to the global mean coverage

Additionally, JLOH allows you to pass a list of known variants that will be used to filter the output depending on the option you set with `--t0-filter-type` ("keep" or "remove").

![JLOH workflow](images/workflow.png)

## Install

As simple as: `git clone https://github.com/Gabaldonlab/jloh.git` or `https://github.com/MatteoSchiavinato/jloh.git`
And it's ready to go! But there are a few dependencies:

| Program     | Type        | Version | Links      |
|-------------|-------------|---------|------------|
| Bedtools    | Program     | 2.25    | [source](https://bedtools.readthedocs.io/en/latest/), [cite](https://doi.org/10.1093/bioinformatics/btq033) |
| pybedtools  | Module      | 0.8.2   | [source](https://daler.github.io/pybedtools/main.html), [cite](https://doi.org/10.1093/bioinformatics/btr539) |
| pysam       | Module      | 0.1.7   | [source](https://pypi.org/project/pysam/), [cite](https://github.com/pysam-developers/pysam) |
| Python      | Interpreter | 3.6.1   | [source](https://www.python.org/downloads/release/python-397/), [cite](http://citebay.com/how-to-cite/python/) |

Note that **pybedtools** will look for **bedtools** in the `$PATH`.

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

The program has many other options:

```
Usage:
./jloh --vcf <VCF> --genome-file <GENOME_FILE> --bam <BAM> [options]

[mandatory]
--vcf               Input VCF file containing all types of variants             [!]
--genome-file       File with chromosome lengths (chromosome <TAB> size)        [!]
--bam               BAM file used to call the --vcf variants                    [!]

[optional]
--sample            Sample name / Strain name for output files                  [loh_blocks]
--output-dir        Output directory                                            [loh_blocks_out]
--t0-vcf            VCF with variants to ignore from --vcf                      [off]
--t0-bam            BAM file used to call the --t0-vcf variants                 [off]
--t0-filter-type    What to do with t0 LOH events? "keep" or "remove"           [remove]
--min-het-snps      Min. num. of heterozygous SNPs in heterozygous region       [2]
--snp-distance      Min. distance (bp between SNPs for blocks definition        [100]
--min-size          Min. LOH block size                                         [100]
--block-lengths     Comma-sep. list of desired distances between LOH blocks     [100,1000,5000]
--min-af            Min. allele frequency to consider a variant heterozygous    [0.3]
--max-af            Max. allele frequency to consider a variant heterozygous    [0.7]
--min-frac-cov      Min. mean coverage fraction for LOH blocks                  [0.75]
                    (used only if --bam specified)
--max-frac-cov      Max. mean coverage fraction for LOH blocks                  [1.25]
                    (used only if --bam specified)
--print-info        Show authors and edits with dates                           [off]
```

## Output

### In Short

**JLOH** produces two output files and places them in the directory specified with `--output-dir`. Everything else is placed inside the `process` directory which, in turn, is inside `--output-dir`. The two output files are:

- `*.het_snps.vcf`: VCF file containing only the selected heterozygous SNPs that are used to extract LOH blocks. If the user specified a `--t0-vcf`, these will be two files: one deriving from `--vcf` (containing the `exp` keyword in the filename) and one derived from `-t0-vcf` (containing the `t0` keyword).

- `*.LOH_blocks.filt.bed`: BED file containing the identified candidate LOH blocks. The BED file has 5 columns: chromosome, start position, end position, coverage change, length. The "coverage change" indicates the LOH block's coverage relative to the global mean coverage (i.e. if it says 2.1, it means that the coverage of this LOH block was 210% of the global mean coverage).

Note that a proper BED file would not have coverage and length on the 4th and 5th column. If using this output file as input to another program, make sure to remove these two columns (e.g. `cut -f 1,2,3`).

### All output explained

#### selection of heterozygous SNPs

The variants passed with `--vcf` are filtered, retaining only heterozygous SNPs. Indels and homozygous SNPs are filtered out as they aren't normally used to extract LOH blocks. The selection of heterozygous SNPs is conducted based on their FORMAT field (field number 9 and 10 of a VCF file). The first column of this field carries a series of annotations separated by colons (e.g. GT:AF) the values of which are annotated the same way on the second column (e.g.` 0/1:0.60`). If a SNP is annotated as heterozygous, it will carry a genotype (`GT`) such as `0/1` or `1/2`. It should also have an allele frequency (`AF`) annotation. SNPs are considered heterozygous if their `AF` annotation falls between the values specified with `--min-af` and `--max-af`. The heterozygous SNPs are then saved in a VCF file called `<sample>.het_snps.vcf`, which is one of the two output files of the program. This file is placed in `--output-dir`.

#### extraction of heterozygosity regions

The putative heterozygosity regions are extracted based on the number of heterozygous SNPs they contain, and the minimum distance between them. The minimum number of SNPs eliminates regions with too little SNPs to be considered true positives, while the minimum SNP distance rules out read mapping artifacts. These parameters are controlled with the `--min-het-snps` and the `--snp-distance` parameters. This step produces three temporary files in the `<output_dir>/process` folder:

1) `<sample>.d<snp_distance>bp_provisory.bed`: BED intervals surrounding heterozygous SNPs, including the SNP count inside the interval.
2) `<sample>.d<snp_distance>bp_provisory.w_len.bed`: the same file but containing the length of each interval.
3) `<sample>.d<snp_distance>bp.bed`: BED intervals defining regions with sufficient number of SNPs (`--min-het-snps`) and sufficiently far apart (`--snp-distance`).

Note that the SNP distance has nothing to do with the length of an interval.

#### extract complementary homozygous regions

The objective of this tool is to extract LOH blocks, which are defined by the *loss* of  heterozygosity. Hence, the tool at this point selects the *complementary* intervals of the heterozygous ones. This step produces a temporary file in the `<output_dir>/process` folder called `<sample>.d<snp_distance>bp.complement.bed`. These regions are then filtered by their length, retaining only those with a length larger or equal to `--min-size`. This produces a temporary file in the `<output_dir>/process` folder called `<sample>.homo.d<snp_distance>bp.bed`.

#### filter by coverage

Each region that is considered as a candidate LOH region is screened by coverage using the BAM file passed with `--bam`. First, a mean coverage is computed for the whole BAM file. To do that, JLOH checks if the BAM file is indexed, and if not, it indexes it using the **pysam** module. Then, reads mapping inside each region are extracted using the **pysam** module, and compared against the general mean coverage. Candidate LOH blocks are retained in the output only if they have a coverage *below* or *above* what is considered "normal" coverage. The "normal" coverage range is defined via two parameters: `--min-frac-cov` and `--max-frac-cov`. The global mean coverage computed from the BAM file is multiplied by these two values to obtained the two boundaries of "normal" coverage. Every candidate LOH block falling above or below these two thresholds will be considered a candidate LOH block. LOH blocks are placed in a temporary file in the `<output_dir>/process` folder  called `<sample>.LOH_blocks.bed`.

#### deal with known pre-existing LOH blocks

If the user passes another VCF file with the `-t0-vcf` option, this file will be considered as variation that pre-dates the one listed in `--vcf`. This "t0" variation is used to extract LOH blocks the same way as described for the input VCF file. The user can then choose what to do with it with the `--t0-filter-type` option. The default ("remove") will remove any overlapping LOH block found between the "t0" VCF and the input VCF. This reduces the output LOH blocks only to those unique to the input VCF. Some users however may want to keep only blocks that are found in "t0" too. To do that, one just has to specity `--t0-filter-type keep`. In both cases, an output file is produced, called `<sample>.LOH_blocks.filt.bed`. This is the real output of the program, and is placed in `--output-dir`.

If the user did not specify any `--t0-vcf`, this file will be identical to `<sample>.LOH_blocks.bed`.
