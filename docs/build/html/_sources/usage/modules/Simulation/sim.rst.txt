.. _jloh-sim:

jloh sim
========

Description
-----------

Simulate a divergent copy of a genome, introducing mutations and loss-of-heterozygosity (LOH).

Usage
-----

.. code-block:: bash 

    jloh sim --fasta <FASTA> [options]

Parameters
----------

Input / Output 
^^^^^^^^^^^^^^

.. function:: --fasta <PATH>

    Input FASTA file to generate mutations and LOH on.

.. function:: --out-fasta <PATH>

    Path to save the output mutated FASTA file. 

.. function:: --out-haplotypes <PATH>

    Path to save the output table with the coordinates of the produced haplotypes. 

Other options
^^^^^^^^^^^^^

.. function:: --threads <INT>

    Number of parallel operations. 

.. function:: --mean-haplotype-size <INT>

    Generated haplotypes will have this average size. 

.. function:: --min-haplotype-length <INT>

    Minimum length of generated haplotypes. 

.. function:: --divergence <FLOAT>

    Average divergence to produce in the output mutated genome [0.0 - 1.0].

.. function:: --loh <FLOAT>
    
    Apply this level of LOH to the mutated genome [0.0 - 1.0].

.. function:: --chrom-name-replace [<STR> <STR>]

    Two space-separated strings that define pattern and replacement to search and substitute in the chromosome names. The output FASTA file will have chromosome names that contain the replacement string instead of the pattern string. 