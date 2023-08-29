# A guide to interpret quantiles produced by JLOH stats 

The **stats** module of JLOH produces an output that describes the quantiles of SNP density obtained for: 
- All SNPs 
- Heterozygous SNPs 
- Homozygous SNPs 

These quantiles are modelled on the SNP density distributions from genome regions of size `--window-size`. If Q10 (quantile 10) says 5 heterozygous SNPs and 2 homozygous SNPs, it means that 90% of the genome has *at least* 5 het SNPs/kbp and 2 homo SNPs/kbp. 

If we set the thresholds to Q10, and Q10 is 5 het SNPs/kbp and 2 homo SNPs/kbp (`--min-snps-kbp 5,2`), it means that we are considering true LOH blocks only those that have less than 5 heterozygous SNPs/kbp, and we are assigning to the alternative allele only those that have at least 2 homozygous SNPs/kbp.