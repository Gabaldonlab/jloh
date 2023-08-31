.. _model-snp-density:

Modelling SNP density 
=====================

.. image:: ../images/snp_density.png


How to do it
------------

A crucial aspect of using JLOH is modelling SNP density. This is done with :ref:`jloh-stats`, which profiles the density of heterozygous and homozygous SNPs over the genome in windows of adjustable size (``--window-size`` and ``--step-size``). This step is aimed at choosing an appropriate SNP density threshold to infer LOH blocks when running :ref:`jloh-extract`.

To infer the density distribution of the SNPs over the genome, one must run: 

.. code-block:: bash 

    jloh stats --vcf my_variants.vcf


This step will produce an output that looks like this::

    -- SNPs/Kbp Statistics --

    S       Gen   Het   Homo
    Median  5.0   4.0   4.0
    Mean    8.5   8.5   5.6
    Max     103   103   110
    Min     1     1     1

    -- SNPs/Kbp Quantiles --

    Q     Gen    Het    Homo
    5%    1.0    1.0    1.0
    10%   1.0    1.0    1.0
    15%   1.0    1.0    2.0
    50%   5.0    4.0    4.0
    85%   15.0   15.4   9.0
    90%   20.0   21.0   11.0
    95%   30.0   30.0   15.0


The user must choose a **Het** and a **Homo** value from a quantile (e.g. Q10) or a descriptor (e.g. median). The values will be used as thresholds to separate blocks into candidate LOH blocks and not. The two chosen values (e.g. 8 and 5) will become the arguments of the ``--min-snps-kbp`` parameter of :ref:`jloh-extract` (see below). 

.. note:: 

    Higher quantiles are stricter in the detection of LOH blocks, they increase the precision but reduce the recall. 

.. tip:: 

    This step tells you how heterozygous and how homozygous your dataset is. Depending on that, you may choose the quantiles appropriately to minimize false positives. 


Explanation
-----------

A quantile cuts a distribution in two. For example, "quantile 10" means that 10% of the items in the distribution have a value that is lower than the one desccribed by Q10. 

If we're talking about SNP densities, a Q10 = 5 means that 10% of the windows in the analysed genome have less than 5 SNPs/kbp. 

In a genomic scenario, we can expect the majority of the genome to float around a certain average het/homo SNPs/kbp value. Then, there will be regions of very low heterozygosity (candidate LOH blocks), and regions of high homozygosity (candidate LOH blocks with an alternative allele). These are the ones we want to select for. 

This is why :ref:`jloh-stats` shows the **mean** and the **median** of the distribution together with the quantiles. The easiest approach is to take the mean as threshold: any window that has less than that many SNPs/kbp will be a candidate LOH block.

However, for greater precision, one can choose to use a quantile (Q15, Q10, Q5). these will be lower than the average and still include 15%, 10%, or 5% of the total windows of the genome. 

Displayed higher quantiles (Q85, Q90, Q95) are for extreme cases. Hypothetically, if about 90% of a genome has undergone LOH, one could set the threshold at Q90 to make sure that nothing is left behind. Note that using high quantiles in case studies with "normal" or low heterozygosity will massively introduce false positives. 