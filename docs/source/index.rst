.. jloh documentation master file, created by
   sphinx-quickstart on Tue Aug 29 11:07:04 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the JLOH documentation
=================================

**What is JLOH?**


**JLOH** (`/jay'lo/`) is a tool to extract, filter, and analyse loss of heterozygosity (LOH) blocks based on single-nucleotide polymorphisms (VCF), read mapping information (BAM), and a reference genome sequence (FASTA). 

.. image:: graphical_abstract.jpg
 
.. toctree::
   :maxdepth: 1
   :caption: Get started:

   usage/installation
   usage/quickstart
   usage/run_test_data
   usage/default_wf
   usage/hybrid_wf

.. toctree::
   :maxdepth: 1
   :caption: Modules description:

   usage/modules/stats
   usage/modules/g2g
   usage/modules/extract
   usage/modules/filter
   usage/modules/intersect 
   usage/modules/cluster 
   usage/modules/chimeric 
   usage/modules/junctions 
   usage/modules/plot 
   usage/modules/sim 