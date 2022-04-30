.. image:: /ONT_logo.png
  :width: 800
  :alt: [Oxford Nanopore Technologies]
  :target: https://nanoporetech.com/

******************

Sockeye
"""""""""

Sockeye is a research Snakemake pipeline designed to identify the cell barcode
and UMI sequences present in nanopore sequencing reads generated from 10X
single-cell libraries.

The inputs are raw nanopore reads (FASTQ) generated from the sequencing
instrument and reference files that can be downloaded from `10X
<https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>`_.
The pipeline outputs a gene x cell expression matrix, as well as a BAM file of
aligned reads tagged with cell barcode and UMI information.

Detailed documentation for all Sockeye commands and algorithms can be found on
the `Sockeye documentation page <https://nanoporetech.github.io/sockeye/>`_.

Prerequisites
-------------

``conda`` must be installed in order to create the base environment where the
Sockeye snakemake pipeline will run. Installation instructions can be found in
the conda `documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

Installation
------------

The project source code must first be cloned from the Oxford Nanopore repository
on GitHub:

::

   git clone git@github.com:nanoporetech/sockeye.git
   cd sockeye

Next you must create and activate the ``conda`` environment (named ``sockeye``)
that contains the necessary packages for calling the Snakemake pipeline:

::

   conda env create -f environment.yml
   conda activate sockeye

Getting Started
---------------

Prior to demultiplexing any nanopore reads, sample sheet information and
relevant pipeline configurations must be specified:

Editing the sample sheet
^^^^^^^^^^^^
``config/runs.csv``

Setting up the pipeline
^^^^^^^^^^^^^^^^^^^^^^
``config/config.yml``
