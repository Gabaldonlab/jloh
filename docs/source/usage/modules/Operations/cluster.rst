.. _jloh-cluster:

jloh cluster
============

Description
-----------

Cluster multiple :ref:`jloh-extract` runs by overlap of their LOH blocks. 

Usage
-----

.. code-block:: bash 

   jloh cluster [options] --loh <LOH_A.tsv> ... <LOH_N.tsv> 

Parameters
----------

.. function:: --loh <PATH_1> ... <PATH_n>

   Space-separated list of \*.LOH_blocks.tsv files produced by :ref:`jloh-extract`. 

.. function:: --out-prefix <PATH|STR>

   Pre-pend this prefix to each output file (path allowed). 

.. function:: --max-dist <FLOAT>

   Maximum distance [0.0 - 1.0] between elements in a cluster. 

.. function:: --threads <INT>

   Number of parallel operations performed. 