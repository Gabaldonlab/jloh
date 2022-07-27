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

![Example](../images/example.png)
