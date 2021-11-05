# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs).

## Install

As simple as: `git clone https://github.com/MatteoSchiavinato/jloh`
And it's ready to go!

### Dependencies

| Bedtools    | Program     | 2.25    | [source](https://bedtools.readthedocs.io/en/latest/), [cite](https://doi.org/10.1093/bioinformatics/btq033) |
| pandas      | Module      | 1.3.4   | [source](https://pandas.pydata.org/), [cite](https://doi.org/10.5281/zenodo.5574486) |
| pysam       | Module      | 0.1.7   | [source](https://pypi.org/project/pysam/), [cite](https://github.com/pysam-developers/pysam) |
| Python      | Interpreter | 3.9.7   | [source](https://www.python.org/downloads/release/python-397/), [cite](http://citebay.com/how-to-cite/python/) |

## Run

The basic usage of the program is as simple as:

```
./jloh --vcf <VCF> --genome-file <GENOME_FILE> [options]
```

The program has many other options:

```
J LOH
Still the one from the block
-
Extact LOH blocks from SNPs and a reference genome
--------------------------------------------------------------------------------
by Matteo Schiavinato
based on the work of Leszek Pryszcz and Veronica Mixao
DOIs:
Pryszcz et al., 2014	https://doi.org/10.1093/gbe/evu082
Mixao et al., 2019		https://doi.org/10.3389/fgene.2019.00383
--------------------------------------------------------------------------------

Usage:
./jloh --vcf <VCF> --genome-file <GENOME_FILE> [options]

[mandatory]
--vcf               Input VCF file containing all types of variants             [!]
--genome-file       File with chromosome lengths (chromosome \t size)           [!]

[optional]
--sample            Sample name / Strain name for output files                  [loh_blocks]
--output-dir        Output directory                                            [loh_blocks_out]
--window            Window for heterozygous blocks definition                   [100]
--min-size          Min. LOH block size                                         [100]
--min-af            Min. allele frequency to consider a variant heterozygous    [off]
--max-af            Max. allele frequency to consider a variant heterozygous    [off]
--bam               BAM file (only required when filtering by coverage)         [off]
--min-frac-cov      Min. mean coverage fraction for LOH blocks                  [off]
--max-frac-cov      Max. mean coverage fraction for LOH blocks                  [off]
--bedtools          Path to the bedtools executable                             [bedtools]
--print-info        Show authors and edits with dates                           [off]
```

## Output

...
