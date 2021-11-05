# J LOH

*[Still the one from the block](https://www.youtube.com/watch?v=dly6p4Fu5TE)*

A tool to extract blocks of loss of heterozygosity (LOH) based on single-nucleotide polymorphisms (SNPs).

## Install

## Run

The program has many options, but only a few are mandatory:

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
--vcf             Input VCF file containing all types of variants               [!]
--genome-file     File with chromosome lengths (chromosome \t size)             [!]

[optional]
--sample          Sample name / Strain name for output files 					          [loh_blocks]
--output-dir      Output directory 											                        [loh_blocks_out]
--windowq         Window for heterozygous blocks definition                     [100]
--min-size        Min. LOH block size 										                      [100]
--min-af          Min. allele frequency to consider a variant heterozygous      [off]
--max-af          Max. allele frequency to consider a variant heterozygous      [off]
--bam             BAM file (only required when filtering by coverage)           [off]
--min-frac-cov		Min. mean coverage fraction for LOH blocks                    [off]
--max-frac-cov		Max. mean coverage fraction for LOH blocks                    [off]
--bedtools        Path to the bedtools executable                               [bedtools]
--print-info      Show authors and edits with dates                             [off]
```

## Output
