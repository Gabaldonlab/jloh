.. _jloh-extract:

jloh extract
============

Description
-----------

Extract candidate LOH blocks from single-nucleotide polymorphisms (SNPs) called from reads mapped onto a reference genome. Alternatively, from reads derived from a hybrid mapped onto its parental genomes (at least one). 

Usage
-----

.. code-block:: bash 

    jloh extract --vcf <VCF> --ref <FASTA> --bam <BAM> [options]


Or, if using ``--assign-blocks``:

.. code-block:: bash 

    jloh extract --assign-blocks --vcfs [<PATH_1> <PATH_2>] --refs [<PATH_1> <PATH_2>] --bams [<PATH_1> <PATH_2>] [options]

Parameters
----------

Default mode
^^^^^^^^^^^^

.. function:: --vcf <PATH>

    VCF file containing single-nucleotide polymorphisms (SNPs).

.. function:: --bam <PATH>

    BAM file containing read mapping records. 

.. function:: --ref <PATH>

    FASTA file where reads were mapped. 

Assign-blocks mode (hybrids)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: --assign-blocks 

    Activate block assignment mode. 

.. function:: --vcfs [<PATH_1> <PATH_2>]

    VCF file containing single-nucleotide polymorphisms (SNPs).

.. function:: --bams [<PATH_1> <PATH_2>]

    BAM file containing read mapping records. 

.. function:: --refs [<PATH_1> <PATH_2>]

    FASTA file where reads were mapped. 

Common Parameters
^^^^^^^^^^^^^^^^^

Variants
********

.. function:: --min-snps-kbp [<INT>,<INT>]

    Comma-separated set of two integer values defining heterozygous and homozygous minimum SNPs/kbp densities. See details at :ref:`model-snp-density`.

.. function:: --filter-mode ["all"|"pass"]

    Either "all" or "pass". Whether to select only VCF entries that have the ``PASS`` annotation or not.

.. function:: --min-af <FLOAT>

    Minimum allele frequency to consider a SNP heterozygous. Useful when working with polyploid species. 

.. function:: --max-af <FLOAT>

    Maximum allele frequency to consider a SNP heterozygous. Useful when working with polyploid species. 

Blocks
******

.. function:: --min-length <INT>

    Minimum length of accepted candidate LOH blocks. 

.. function:: --coarseness <INT>

    Minimum length of initial building blocks that are used to define LOH blocks (i.e. nothing shorter than this will be considered an interesting interval).

.. function:: --min-frac-cov <FLOAT>

    Minimum fraction of positions of a candidate LOH block to include it in the final list. 

.. function:: --hemi <FLOAT>

    Threshold of coverage ratio between candidate block and surrounding up/downstream regions, below which a block is considered hemizygous (i.e. carrying only one copy).

.. function:: --overhang <INT>

    Size of the up/downstream region checked to define zygosity (see ``--hemi``).

.. function:: --min-overhang <FLOAT>

    Fraction of the ``--overhang`` that must be present to infer zygosity (e.g. at the beginning of a chromosome).

.. function:: --merge-uncov <INT>

    Number of uncovered positions (bp) separating two blocks that are ignored, producing a merged block. 

Misc
****

.. function:: --sample <STR>

    Sample name to include in output files.

.. function:: --output-dir <PATH>

    Path to the output directory. 

.. function:: --threads <INT>

   Number of parallel operations performed. 

.. function:: --regions <PATH>

    BED file containing regions where blocks shall be searched in. This BED file may be created via :ref:`jloh-g2g` or it may be a custom BED file with regions of interest.

.. function:: --os-scratch

    Use the system's temporary directory instead of creating temporary files in the output path.