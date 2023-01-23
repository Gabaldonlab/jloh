# A guide to interpret quantiles produced by JLOH stats 

The **stats** module of JLOH produces an output that describes the quantiles of SNP density obtained for: 
- All SNPs 
- Heterozygous SNPs 
- Homozygous SNPs 

These quantiles define the percentage of genome that has that many SNPs/kbp or more. For example, if Q10 (quantile 10) says 5 heterozygous SNPs and 2 homozygous SNPs, it means that 90% of the genome (100 - 10) has *at least* 5 het SNPs/kbp and 2 homo SNPs/kbp. 

Another way to look at SNP density quantiles is to think of them as the percentage (%) of the genome that *you would cut away* by setting the parameter to that value. If we set the thresholds to Q10, and Q10 is 5 het SNPs/kbp and 2 homo SNPs/kbp (`--min-snps-kbp 5,2`), it means that we are subdividing the genome in two groups: one below and one above the threshold, the first corresponding to *noise* and the second corresponding to *signal*.