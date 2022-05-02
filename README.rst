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

Setting up the pipeline
^^^^^^^^^^^^^^^^^^^^^^

The pipeline configurations are described in a YAML file

``config/config.yml``

::

   SAMPLE_SHEET: "./config/samples.csv"

   OUTPUT_BASE: /PATH/TO/OUTPUT/BASE/DIRECTORY

   ################################################################################
   # 10x SUPPORTING FILES                                                         #
   ################################################################################
   # Reference files can be downloaded from the 10x website using either curl or wget:
   # For the human GRCh38 reference, the commands would be:
   # curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
   # or
   # wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz

   ######### REF_GENOME_DIR #########
   # REF_GENOME_DIR refers the path to reference directory as downloaded from 10x,
   # e.g. PATH/TO/10X/DOWNLOADS/refdata-gex-GRCh38-2020-A for a human reference.
   REF_GENOME_DIR: PATH/TO/10X/DOWNLOADS/refdata-gex-GRCh38-2020-A

   ######### BC_SUPERLIST #########
   # The 10x cell barcode full whitelist (BC_SUPERLIST) can be downloaded from:
   # wget https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz
   # BC_SUPERLIST path can point to either the .txt.gz or .txt file
   BC_SUPERLIST: PATH/TO/10X/DOWNLOADS/3M-february-2018.txt
   ################################################################################

   MAX_THREADS: 4

   READ_STRUCTURE_BATCH_SIZE: 40000
   READ_STRUCTURE_BARCODE_LENGTH: 16
   READ_STRUCTURE_UMI_LENGTH: 12
   READ_STRUCTURE_READ1: CTACACGACGCTCTTCCGATCT
   READ_STRUCTURE_TSO: ATGTACTCTGCGTTGATACCACTGCTT
   READ_STRUCTURE_FLAGS: ""

   BARCODE_READ1_SUFF_LENGTH: 10
   BARCODE_KNEEPLOT_FLAGS: ""
   BARCODE_MAX_ED: 2
   BARCODE_MIN_ED_DIFF: 2

   GENE_ASSIGNS_MINQV: 60

   UMI_GENOMIC_INTERVAL: 1000
   UMI_CELL_GENE_MAX_READS: 20000
   UMI_CLUSTER_MAX_THREADS: 4

   MATRIX_MIN_GENES: 100
   MATRIX_MIN_CELLS: 3
   MATRIX_MAX_MITO: 5
   MATRIX_NORM_COUNT: 10000

   # Using a comma-separated list, specify which genes should be annotated in the
   # UMAP plots (e.g. CD19,PAX5,XBP1)
   UMAP_PLOT_GENES: CD19,CD24,CD27,CD38,CD79A,CD79B,PAX5,XBP1

   # Set the maximum resources to devote to the minimap2 alignment step
   RESOURCES_MM2_MEM_GB: 50
   RESOURCES_MM2_MAX_THREADS: 4

Editing the sample sheet
^^^^^^^^^^^^
The path to the sample sheet is defined in the ``config.yml`` file described above. This sample sheet contains details about the input run IDs and ONT read directory. The input read directory specified in the sample sheet can contain multiple ``*.fastq``, ``*.fq``, ``*.fastq.gz`` or ``*.fq.gz`` files, but all file extensions must be the same. A mixture of file extensions is not supported.

``config/samples.csv``
