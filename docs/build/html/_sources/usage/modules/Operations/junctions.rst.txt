.. _jloh-junctions:

jloh junctions
==============

Description
-----------

Calculate number of block-to-block junctions over the genome.

Usage
-----

.. code-block:: bash 

    jloh junctions --blocks [<TSV> ...] [options]

Parameters
----------

Input / Output
^^^^^^^^^^^^^^

.. function:: --blocks <PATH> | [<PATH_1> <PATH_2>]

    This input parameter works in two different ways, both with \*.LOH_blocks.tsv files produced by :ref:`jloh-extract`.

    - **one input file**: Junctions are searched within the file between blocks of different alleles (REF, ALT).
    - **two input files**: Junctions are searched between the two files. 

.. function:: --gff <PATH>

    GFF file containing gene models that will be used to produce more in-depth statistics in the output. 

.. function:: --genome <PATH>

    FASTA file with the genome sequence where blocks were called. 

Other options
^^^^^^^^^^^^^

.. function:: --max-dist <INT>

    Maximum distance between LOH blocks to consider them part of a candidate junction. 