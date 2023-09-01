.. _jloh-chimeric:

jloh chimeric
=============

Description
-----------

This module extracts genes from a genome annotation file which contain traces of LOH towards two different haplotypes (REF and ALT). 

Biologically, these are likely the result of an LOH breakpoint placed within an intron of a gene, so that part of the gene has the REF haplotype, and part of the gene has the ALT one. 

Usage
-----

.. code-block:: bash 

    jloh chimeric --blocks-A <TSV> --blocks-B <TSV> --het <BED> --gff <GFF> [options]

Parameters
----------

.. function:: --blocks-A <PATH>

    A TSV file produced by :ref:`jloh-extract`. 
    
.. function:: --blocks-B <PATH>

    Another TSV file produced by another run of :ref:`jloh-extract`, profiling a different genotype from ``--blocks-A``. 
    
.. function:: --het <PATH>

    Heterozygous regions identified by :ref:`jloh-extract`. 
    
    If the extraction was done with ``--assign-blocks`` it has produced two files with Heterozygous blocks (one per parental genome). In that case, concatenate them together. 

.. function:: --gff <PATH>

    GFF3 file containing gene models from the reference FASTA file used to infer LOH blocks. 

.. function:: --out-prefix <PATH/STR>

    Pre-pend this prefix to each output file (path allowed). 

.. function:: --quiet <BOOL>

    Suppress warnings.

.. function:: --feature <STR>

    What feature to look for (see GFF format, 3rd column). Examples are "gene", "exon", "CDS", "mRNA".

.. function:: --min-overlap <FLOAT>

    Minimum fraction [0.0 - 1.0] of a feature that must overlap an LOH block to create an instance of candidate chimeric feature. 

.. function:: --id-attr <STR>

    Attribute to search for to get the names of the features analysed (GFF3 last column, e.g. "ID").