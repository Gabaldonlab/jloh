.. _jloh-intersect:

jloh intersect
==============

Description
-----------

Perform intersection/removal operations with output files of :ref:`jloh-extract`.

Usage
-----

.. code-block:: bash 

    jloh intersect --loh-A <LOH_A.tsv> --loh-B <LOH_B.tsv> [options]

Parameters
----------

Input / Output
^^^^^^^^^^^^^^

.. function:: --loh-A <PATH>

    \*.LOH_blocks.tsv output file produced by :ref:`jloh-extract`. 

.. function:: --loh-B <PATH>

    Another \*.LOH_blocks.tsv output file produced by :ref:`jloh-extract`. 

Other options
^^^^^^^^^^^^^

.. function:: --mode <STR>

    What to do with the two input files. Choices are: 

    - **intersection** (default): keep LOH blocks in common between the two files 
    - **complement**: keep lines unique to ``--loh-A``
    - **unique**: exclude common lines between the two files 

.. function:: --min-ovl <FLOAT>

    Minimum overlap [0.0 - 1.0] between blocks to validate operation requested in ``--mode``. 