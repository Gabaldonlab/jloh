Installation
============

Docker
------

The latest stable version of JLOH is available as a docker image in `Docker Hub <https://hub.docker.com/repository/docker/cgenomics/jloh>`_. To obtain it, make sure you have a functioning `Docker Engine <https://docs.docker.com/engine/install/>`_ in your system. If that's the case, simply run this command from your terminal:

.. code-block:: bash 

    docker pull cgenomics/jloh

First run
^^^^^^^^^

If everything went smoothly, you should be able to run JLOH simply by running: 

.. code-block:: bash

    docker run --rm -t -i cgenomics/jloh --help



.. note:: 
    The docker image downloaded will contain all the dependencies and will function as is (no need to read further).


GitHub
------

The latest development version of JLOH is available in `GitHub <https://github.com/Gabaldonlab/jloh.git>`_. Below find the instructions on how to proceed with the installation from GitHub.

Latest release
**************

To get the latest *stable* release, open the `JLOH GitHub page <https://github.com/Gabaldonlab/jloh.git>`_ and navigate to the releases section on the right. 

.. note:: 
    Once you have downloaded it, make sure you install all the dependencies listed below before using it.

From source
***********

To get the latest *unstable* (development) release, make sure you have a functioning `Git Engine <>`_ in your system. If that's the case, simply clone it in the folder you prefer: 

.. code-block:: bash

    git clone https://github.com/Gabaldonlab/jloh.git

.. note:: 
    Once you have downloaded it, make sure you install all the dependencies listed below before using it.


Install dependencies
^^^^^^^^^^^^^^^^^^^^

These are the dependencies that **jloh** requires in order to function properly.

+------------------+------------------+------------+------------------------------------------------+
| Program          | Type             | Version    | Links                                          |
+==================+==================+============+================================================+
| all2vcf          | Program          | 0.7.3      | https://github.com/MatteoSchiavinato/all2vcf   |
+------------------+------------------+------------+------------------------------------------------+
| bedtools         | Program          | 2.30       | https://bedtools.readthedocs.io/en/latest      |
+------------------+------------------+------------+------------------------------------------------+
| Biopython        | Module           | 1.79       | https://biopython.org/                         |
+------------------+------------------+------------+------------------------------------------------+
| MUMmer           | Program          | 3.1        | https://mummer4.github.io/install/install.html |
+------------------+------------------+------------+------------------------------------------------+
| numpy            | Module           | 1.21.4     |https://numpy.org/                              |
+------------------+------------------+------------+------------------------------------------------+
| pandas           | Module           | 1.3.5      | https://pandas.pydata.org/                     |
+------------------+------------------+------------+------------------------------------------------+
| pybedtools       | Module           | 0.8.2      | https://daler.github.io/pybedtools/            |
+------------------+------------------+------------+------------------------------------------------+
| pandarallel      | Module           | 1.6.1      | https://pypi.org/project/pandarallel/          |
+------------------+------------------+------------+------------------------------------------------+
| pysam            | Module           | 0.1.7      | https://pypi.org/project/pysam/                |
+------------------+------------------+------------+------------------------------------------------+
| Python           | Interpreter      | 3.6.1      | https://www.python.org/downloads/              |
+------------------+------------------+------------+------------------------------------------------+
| samtools         | Program          | 1.13       | http://www.htslib.org/download/                |
+------------------+------------------+------------+------------------------------------------------+

.. Note:: 
    The **pybedtools** python module will look for the `bedtools <https://bedtools.readthedocs.io/en/latest/>`_ executable in the ``$PATH``, while the **pysam** python module will look for the `samtools <http://www.htslib.org/download/>`_ executable in the ``$PATH``. Moreover, **MUMmer** and its tools (e.g. ``nucmer``) must be in the ``$PATH`` as well. The installation of MUMmer can be a bit cumbersome in certain systems, hence we recommend you proceed with conda (``conda install -c bioconda mummer``). This will place in the ``$PATH`` all the toolkit from the MUMmer arsenal, in particular the tools needed by JLOH to run: ``nucmer``, ``delta-filter``, ``show-snps``. 


First run
^^^^^^^^^

If everything went smoothly, you should be able to run JLOH simply by running: 

.. code-block:: bash 

    ./jloh --help 

