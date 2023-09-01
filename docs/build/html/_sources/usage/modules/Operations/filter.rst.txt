.. _jloh-filter:

jloh filter
===========

Description
-----------

Filter produced LOH blocks based on different criteria. 

Usage
-----

.. code-block:: bash 

    jloh filter --loh <LOH.tsv> [options]

Parameters
----------

Input/Output
^^^^^^^^^^^^

.. function:: --loh <PATH>

    A \*.LOH_blocks.tsv file produced by :ref:`jloh-extract`.

Basic filters
^^^^^^^^^^^^^

.. function:: --length <INT>

    Minimum length (bp) of blocks to keep. 

.. function:: --alleles ["REF"|"ALT"|"NA"]

    Allele annotation of blocks to keep. Multiple choices possible if space-separated (e.g. "REF NA").

.. function:: --zygosity ["homo"|"hemi"|"NA"]

    Zygosity of blocks to keep. Multiple choices possible if space-separated (e.g. "hemi NA").

.. function:: --coverage <INT>

    Minimum coverage of blocks to keep. 

.. function:: --region <STR>

    Region where blocks must be kept, declared as *chr:start-end*, or simply *chr*. 

SNP filters 
^^^^^^^^^^^

.. function:: --snps <INT>

    Minimum number of SNPs found in block in order to keep it. 

.. function:: --het-snps <INT>

    Minimum number of heterozygous SNPs found in block in order to keep it. 

.. function:: --homo-snps <INT>

    Minimum number of homozygous SNPs found in block in order to keep it. 

.. function:: --snps-kbp <INT>

    Minimum SNPs/kbp density found in block in order to keep it. 

.. function:: --het-snps-kbp <INT>

    Minimum heterozygous SNPs/kbp density found in block in order to keep it. 

.. function:: --homo-snps-kbp <INT>

    Minimum homozygous SNPs/kbp density found in block in order to keep it. 