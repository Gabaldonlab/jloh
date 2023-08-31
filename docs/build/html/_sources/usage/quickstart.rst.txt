.. _quickstart:

Quick start
===========

If you're in a hurry and can't read too much documentation, here's a quick way to use JLOH with commonly available files derived from NGS data analysis. If, instead, you're working with a hybrid organism, make sure you read this section: :ref:`hybrid-wf`.

Requirements
------------

JLOH is built to work as an addition to any common **variant calling** pipeline. These pipelines are based on two steps: 1) read mapping and 2) variant calling. Regardless of the tool used to map the reads, the most common alignment format is the SAM format (or its binary version BAM). In the calling step, regardless of the tool used to call variants the output is most often in VCF format.

These are the files you need as input to infer LOH blocks. 

+-----------+------------------------------------------------------------------------------------------------------------------+
| File type | Description                                                                                                      |
+===========+==================================================================================================================+
| FASTA     | The reference genome sequence where you map your genomic reads.                                                  |
+-----------+------------------------------------------------------------------------------------------------------------------+
| BAM       | The output of the mapping step.                                                                                  |
+-----------+------------------------------------------------------------------------------------------------------------------+
| VCF       | The file containing all the single-nucleotide polymorphisms (SNPs) called from the BAM file onto the FASTA file. |
+-----------+------------------------------------------------------------------------------------------------------------------+

These three files will highlight the positions in which the genotype represented from the **reads** has lost heterozygosity when compared to the genotype represented by the **reference**. 


Calculating SNP density
-----------------------

This is done with :ref:`jloh-stats`, for details see :ref:`model-snp-density`. Run the command:

.. code-block:: bash 

    jloh stats --vcf my_variants.vcf 


And choose the thresholds of SNP density for heterozygous and homozygous SNPs. 


Extracting blocks
-----------------

The second step is the inference of LOH blocks. This step is done with ``jloh extract``. At the very minimum, the program requires these parameters:

.. code-block:: bash

    jloh extract --vcf my_variants.vcf --bam my_mappings.bam --ref my_reference_genome.fasta


Besides the three input parameters, we encourage you to set at least two other parameters, even though they have default values: 

  *  ``--min-snps-kbp <N,N>``: heterozygous/homozygous minimum SNP/kbp densities to label a region as heterozygous/homozygous 
  *  ``--threads``: number of parallel operations 


Among the output files, there is also a table with the inferred LOH blocks which will end with ``*LOH_blocks.tsv``. This is the main output of the program.