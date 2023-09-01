.. _jloh-plot:

jloh plot
=========

Description
-----------

Plot LOH propensity over all or specific chromosomes from multiple JLOH extract TSV output files.

Usage
-----

.. code-block:: bash 

    jloh plot [--one-ref | --two-ref] --loh <TSV_1> ... <TSV_n> [options]

Parameters
----------

Modes
^^^^^

.. function:: --one-ref <BOOL>

    Plot block propensity towards LOH or heterozygosity. Must be used in combination with ``--het`` parameter.
    This mode compares the regions with candidate LOH blocks with the regions with heterozygous blocks. 

.. function:: --two-ref <BOOL>

    Plot LOH propensity towards genome A or B, with results obtained from :ref:`jloh-extract` ran in ``--assign-blocks`` mode.
    This mode compares the blocks assigned to "REF" with those assigned to "ALT". 

Input / Output 
^^^^^^^^^^^^^^

Input
*****

Choose one of these three options: 

.. function:: --loh [<PATH_1> ... <PATH_n>]

    \*.LOH_blocks.tsv output files produced by :ref:`jloh-extract`. 

.. function:: --loh-files <PATH>

    File of file names (.fofn) containing paths to \*.LOH_blocks.tsv files produced by :ref:`jloh-extract`, one per line. This option is an alternative to ``--loh``.

.. function:: --names-file <PATH>

    TAB-separated file containing paths to \*.LOH_blocks.tsv files produced by :ref:`jloh-extract`, and an associated *sample name*. This option is an alternative to ``--loh``.

.. function:: --het <PATH>

    \*.exp.het_blocks.tsv output file produced by :ref:`jloh-extract`. This file contains annotated heterozygous blocks and is needed only in ``--one-ref`` mode.

Output
******

.. function:: --output-dir <PATH>

    Path to the output directory. 

.. function:: --prefix <STR>

    Prefix to include in output files. 

Plot construction options 
^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: --ref-name <STR>

    Name to use when labelling the reference (REF) allele. 

.. function:: --alt-name <STR>

    Name to use when labelling the alternative (ALT) allele. 

.. function:: --by-sample <BOOL>

    Produce a plot for each sample instead of a plot for each chromosome. 

.. function:: --merge <BOOL>

    The output plot will combine all chromosomes into one. 

.. function:: --clusters <PATH>

    Pass the output file produced by :ref:`jloh-cluster` to sort the samples in the output plot according to the cluster.

.. function:: --threads <INT>

    Number of parallel workers to use in the plot table construction.

.. function:: --chr <STR>

    Restrict the analysis to this specific chromosome. 

.. function:: --chr-file <PATH> 

    Restrict the analysis to these specific chromosomes (file with names, one per line).

.. function:: --window-size <INT>

    Size of the plotting chromosomal windows. Smaller windows are more precise but slower to build. 

.. function:: --contrast <STR>

    Increase plot contrast for samples with low LOH rate. Choose one of: 

    - **off** (default): leave plot contrast untouched. 
    - **low**: set contrast to "low". This scales the LOH propensity to the 0 - 0.25 range.
    - **mid**: set contrast to "mid". This scales the LOH propensity to the 0 - 0.50 range.
    - **high**: set contrast to "high". This scales the LOH propensity to the 0 - 0.75 range.
    - **max**: set contrast to "max". This scales the LOH propensity to the 0 - 1.00 range.

    At the moment, this function works only with the ``--one-ref`` mode.

R / ggplot2 options 
^^^^^^^^^^^^^^^^^^^

More info on some of these options can be found in the `tidyverse official ggplot2 theme() <https://ggplot2.tidyverse.org/reference/theme.html>`_ manual.
Colors are represented in HEX notation, find yours `here <https://htmlcolorcodes.com/>`_.

.. function:: --r-exec <PATH>

    Path to the ``Rscript`` executable.

.. function:: --aspect-ratio <FLOAT>

    y/x length ratio.

.. function:: --width <INT>

    Plot width (px).

.. function:: --height <INT>

    Plot height (px). 

.. function:: --res <INT>

    Plot resolution. 

.. function:: --het-color <STR>

    Color to use to represent heterozygous blocks (``#...`` HEX code).

.. function:: --ref-color <STR>

    Color to use to represent reference LOH blocks (``#...`` HEX code).

.. function:: --alt-color <STR>

    Color to use to represent alternative LOH blocks (``#...`` HEX code).

.. function:: --midpoint-color <STR>

    Color to use to represent the midpoint between ALT and REF (``#...`` HEX code).