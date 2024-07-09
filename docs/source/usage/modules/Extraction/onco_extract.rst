.. _jloh-onco_extract:

jloh onco_extract
=================

Description
-----------

Extract candidate LOH blocks from single-nucleotide polymorphisms (SNPs) called from reads mapped onto a reference genome using cancer data. Data may be a single set of BAM+VCF or a matched tumor/normal pair (recommended). 

Usage
-----

.. code-block:: bash 

    jloh onco_extract --vcfs <VCF_control> <VCF_tumor> --ref <FASTA> --bams <BAM_control> <BAM_tumor> [options]


Or, if using ``--single-mode``:

.. code-block:: bash 

    jloh onco_extract --vcf <VCF> --ref <FASTA> --bam <BAM> [options]


Parameters
----------

Default mode
^^^^^^^^^^^^

.. function:: --vcfs [<PATH_1> <PATH_2>]

    VCF files (space-separated) from control & tumor, in this order. 

.. function:: --bams [<PATH_1> <PATH_2>]

    BAM files from used to call the --vcfs (space-separated).

.. function:: --ref <PATH>

    FASTA file where reads were mapped. 

Single sample mode
^^^^^^^^^^^^^^^^^^

.. function:: --single-mode 

    Activate block assignment mode. 

.. function:: --vcf <PATH_1>

    VCF file containing single-nucleotide polymorphisms (SNPs).

.. function:: --bam <PATH_1>

    BAM file containing read mapping records. 

.. function:: --refs <PATH_1>

    FASTA file where reads were mapped. 

Common Parameters
^^^^^^^^^^^^^^^^^

Input/Output
************

.. function:: --sample <STR>

    Sample name for output files.

.. function:: --output-dir <PATH>

    Path to an output directory (created if not existing)

.. function:: --regions <PATH>

    Path to a BED file with regions of interest. The BED file **must** contain 4 columns: chromosome, start position, end position, and annotation. Annotation may be anything (gene name, transcript name, exon name, locus) as long as it is an alphanumeric string.

Variants
********

.. function:: --max-dist <INT>

    Maximum distance allowed between SNPs for them to still be retained within the same homozygous / heterozygous block.

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

.. function:: --min-snps <INT>

    Minimum number of homozygous SNPs to consider a block in the final results. Homozygous SNPs are an indicator of LOH when a paired normal sample is present. 

.. function:: --min-snps-het <INT>

    Minimum number of heterozygous SNPs to discard a block in the final results. 

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