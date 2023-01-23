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

Then, we run the **JLOH stats** module to calculate the distribution of SNP densities across the genome in windows of adjustable size. Default is 10 kbp. We do this for *both parental genomes*. 

```
jloh stats --vcf A.ff.vcf 
jloh stats --vcf B.ff.vcf 
```

The suggested parameter setting from `jloh stats` is the one correspoonding to the 50th quantile (Q50). We suggest to make a decision based on how the quantiles look. SNP-dense datasets will benefit from higher quantiles to separate real LOH blocks from background noise (for example Q50). More SNP-depleted datasets will, instead, benefit from a lower quantile (e.g. Q10) so to not exclude many true positives.

Details on this rationale are available [here](../docs/QUANTILES.md). 

The final values to be used in `jloh extract` in the `--min-snps-kbp` parameter are two: one corresponding to heterozygous SNPs/Kbp, the other to homozygous SNPs/kbp. When using two parental genomes, you will generate two sets of these values (one set per genome). We suggest you average them and use the two averages as final parameter setting in `jloh extract`. 

For example: 
```
Parent A (Q25): 20 heterozygous, 10 homozygous 
Parent B (Q25): 14 heterozygous, 6 homozygous 

Final parameter setting: 
(20 + 14) / 2 = 17
(10 + 6) / 2  = 8

jloh extract --min-snps-kbp 17,8 <other_options>

```

Finally, we call LOH blocks using **JLOH extract**, limiting the output to regions showing divergence in genome-to-genome mapping.

```
jloh extract \
--threads 4 \
--vcfs A.ff.vcf B.ff.vcf \
--min-snps-kbp <Het>,<Homo> \
--bams A.bam B.bam \
--refs parent_A.fasta parent_B.fasta \
--regions A_and_B.regions.bed \
--sample <STRING> \
--output-dir <PATH>
```

The `<sample>.LOH_blocks.{A,B}.tsv` files contained in `--output-dir` will contain all *bona fide* blocks found with this approach in either A or B genome. The `<sample>.LOH_candidates.{A,B}.tsv` files, instead, contain all the blocks that were called from the program, regardless of the regions specified in `--regions`.
