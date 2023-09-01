.. _jloh-g2g:

jloh g2g
========

Description
-----------

Align two genome sequences to find regions where they diverge. Useful when later running :ref:`jloh-extract` in ``--assign-blocks`` mode.

Usage
-----

.. code-block:: bash 

    jloh g2g --target <FASTA> --query <FASTA> [options]

Parameters
----------

Input / Output
^^^^^^^^^^^^^^

.. function:: --ref-A <PATH>

    First reference genome (FASTA) to use in the alignment.

.. function:: --ref-B <PATH>

    Second reference genome (FASTA) to use in the alignment.

.. function:: --output-dir <PATH>

    Path to the output directory. 

Mapping presets
^^^^^^^^^^^^^^^

For more information on the nucmer parameters meaning, consult `their manual <https://mummer.sourceforge.net/manual/#nucmer>`_.

.. function:: --default <BOOL>

    nucmer parameters: -c 65 -b 200 -l 20

.. function:: --sensitive <BOOL>

    nucmer parameters: -c 100 -b 50 -l 50 --mum

.. function:: --relaxed <BOOL>

    nucmer parameters: -c 50 -b 500 -l 20

Other
^^^^^

.. function:: --est-divergence <FLOAT>

    Estimated divergence between the two aligned genomes. Produced alignments will be retained if diverging up to *twice* the declared divergence. I.e. if the user declares ``--est-divergence 0.05``, alignments will be kept down to 90% identity.

.. function:: --min-length <INT>

    Minimum length of the final alignments produced by nucmer that are kept for further analyses.

.. function:: --all2vcf-exe <PATH>

    Path to the `all2vcf <https://github.com/MatteoSchiavinato/all2vcf>`_ executable. We provide a version of the tool within the jloh installation directory (``src/all2vcf``). 