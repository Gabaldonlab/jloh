.. _one-parent:

Hybrids with only one known parental
====================================

If you're working in a hybrid setup but only have one parental genome at hand, *you can still use JLOH*. There's only a few adjustments you shall take into account. 

To generate the necessary BAM and VCF files, follow the guidelines described in :ref:`hybrid-wf`, doing each step only once (i.e. with the only parental genome you have). You won't be able to run :ref:`jloh-g2g`, since you don't have two parental genomes. 

When running :ref:`jloh-stats`, you won't have any information regarding the homozygosity derived from the missing parental. Hence, you must go forward with the quantile values you choose from the distribution of your only known parental. 

When running :ref:`jloh-extract`, simply run it in default mode (i.e. *without* ``--assign-blocks``) and pass the parental genome you have as ``--ref``. Blocks annotated as ``REF`` did not have many homozygous SNPs and are assumed to carry the allele of the parental genome you have in your possession. Those annotated as ``ALT`` are likely carrying the allele of the other parental.

When running :ref:`jloh-plot`, use it in ``--one-ref`` mode as you would do with a normal (i.e. non-hybrid) run. 