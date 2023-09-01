.. _jloh-stats:

jloh stats
==========

Description
-----------

Calculate heterozygous and homozygous SNP density distributions from a VCF file

Usage
-----

.. code-block:: bash 

    jloh stats --vcf <VCF>

Parameters
----------

Input / Output 
^^^^^^^^^^^^^^

.. function:: --vcf <PATH>

    Input VCF file to produce statistics on. 

Other options
^^^^^^^^^^^^^

.. function:: --threads <INT>

    Number of parallel operations. 

.. function:: --window-size <INT>

    Distribution of SNP density will be calculated on genomic sliding windows of this size. 

.. function:: --step-size <INT>

    Sliding windows used in distribution will increase by this step size.