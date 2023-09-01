.. _home:

.. jloh documentation master file, created by
   sphinx-quickstart on Tue Aug 29 11:07:04 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the JLOH documentation
=================================

**What is JLOH?**

**JLOH** (`/jay'lo/`) is a tool to extract, filter, and analyse loss of heterozygosity (LOH) blocks based on single-nucleotide polymorphisms (VCF), read mapping information (BAM), and a reference genome sequence (FASTA). 

A publication is currently being reviewed. The preprint can be found in bioRxiv (`Schiavinato et al., 2023 <https://doi.org/10.1101/2023.05.04.539368>`_).

.. image:: images/graphical_abstract.jpg
 
.. toctree::
   :maxdepth: 2
   :caption: Get started

   usage/installation
   usage/quickstart
   usage/run_test_data
   usage/model_snp_density
   usage/hybrid_wf
   usage/one_parent
   
.. toctree::
   :maxdepth: 2
   :caption: Extraction tools

   usage/modules/Extraction/stats
   usage/modules/Extraction/g2g 
   usage/modules/Extraction/extract 

.. toctree:: 
   :maxdepth: 2
   :caption: Operations tools 

   usage/modules/Operations/filter 
   usage/modules/Operations/intersect 
   usage/modules/Operations/chimeric
   usage/modules/Operations/junctions  
   usage/modules/Operations/cluster 

.. toctree::
   :maxdepth: 2
   :caption: Visualization tools 

   usage/modules/Visualization/plot 

.. toctree::
   :maxdepth: 2
   :caption: Simulation tools 

   usage/modules/Simulation/sim 

