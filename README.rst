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

The pipeline configurations are described in the YAML file ``config/config.yml``:

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

Most of the parameters defined in the ``config/config.yml`` file can normally remain unchanged. However, certain fields require editing, such as:

::

   OUTPUT_BASE     # Base directory where run_id-specific output folders will be written
   REF_GENOME_DIR  # Path to the downloaded 10X reference data
   BC_SUPERLIST    # Path to the downloaded 10X cell barcode whitelist (i.e. 3M-february-2018.txt.gz)
   MAX_THREADS     # Maximum number of threads to use for various steps in the pipeline
   UMAP_PLOT_GENES # Genes to annotate in UMAP plots

Editing the sample sheet
^^^^^^^^^^^^
The path to the sample sheet is defined by the ``SAMPLE_SHEET`` variable in the ``config.yml`` file described above (set to ``./config/samples.csv`` by default). This sample sheet contains details about the input run IDs and ONT read directory. Sockeye can launch analyses of multiple runs simultaneously, which is useful especially when submitting the analyses to a compute cluster.

The input read directory specified in the sample sheet can contain multiple ``*.fastq``, ``*.fq``, ``*.fastq.gz`` or ``*.fq.gz`` files, but all file extensions must be the same. A mixture of file extensions is not supported.

``config/samples.csv``

::

   run_id,path
   run1,/PATH/TO/ONT/READS1.fq.gz
   run2,/PATH/TO/ONT/READS2.fq.gz
   run3,/PATH/TO/ONT/READS3.fq.gz

Launching Sockeye
^^^^^^^^^^^^^^^^^

Once the Sockeye environment has been created and activated (see Installation above) and both the ``config.yml`` and ``samples.csv`` files have been edited, the Sockeye pipeline is ready to be launched.

Launch Sockeye locally from the Sockeye repository using:

::

   snakemake --use-conda --configfile config/config.yml -pr all

If your cluster system supports Distributed Resource Management Application API (DRMAA), you can submit the Sockeye pipeline to your job scheduler using:
::

   snakemake --configfile config/config.yml --latency-wait 300 --drmaa ' -V -cwd -P applications -l m_mem_free={resources.mem}G -pe mt {threads} ' --default-resources mem=1 --jobs 1000 --use-conda --drmaa-log-dir ./drmaa_logs -pr all

More details on cluster execution for various systems can be found `here <https://snakemake.readthedocs.io/en/stable/executing/cluster.html>`_.

Pipeline output
---------------

The pipeline output will be written to a directory defined by ``OUTPUT_BASE`` in the ``config/config.yml`` file. For instance, using the example ``config/config.yml`` and ``config/sample_sheet.csv`` files shown above, the pipeline output would be written to three separate directories, one for each ``run_id``:

::

   /PATH/TO/OUTPUT/BASE/DIRECTORY/run1
   /PATH/TO/OUTPUT/BASE/DIRECTORY/run2
   /PATH/TO/OUTPUT/BASE/DIRECTORY/run3

Each run_id-specific output folder will contain the following subdirectories:

::

   /PATH/TO/OUTPUT/BASE/DIRECTORY/run1
   |
   |-- adapters   # contains output from the characterization of read structure based on adapters
   |-- align      # output from the alignment to the reference
   |-- demux      # demultiplexing results, primarily in the tagged.sorted.bam file
   |-- matrix     # gene expression matrix and UMAP outputs
   \-- saturation # plots describing the library sequencing saturation

The most useful outputs of the pipeline are likely:

* ``adapters/configs.stats.json``: provides a summary of sequencing statistics and observed read configurations, such as

  - ``n_reads``: number of total reads in the input fastq(s)
  - ``rl_mean``: mean read length
  - ``n_fl``: total number of reads with the read1-->TSO or TSO'-->read1' adapter configuration (i.e. full-length reads)
  - ``n_plus``: number of reads with the read1-->TSO configuration
  - ``n_minus``: number of reads with the TSO'-->read1' configuration

* ``demux/tagged.sorted.bam``: BAM file of alignments to the reference where each alignment contains the following sequence tags

  - CB: corrected cell barcode sequence
  - CR: uncorrected cell barcode sequence
  - CY: Phred quality scores of the uncorrected cell barcode sequence
  - UB: corrected UMI sequence
  - UR: uncorrected UMI sequence
  - UY: Phred quality scores of the uncorrected UMI sequence

* ``matrix/gene_expression.processed.tsv``: TSV containing the gene (rows) x cell (columns) expression matrix, processed and normalized according to the parameters defined in the ``config/config.yml`` file:

  - ``MATRIX_MIN_GENES``: cells with fewer than this number of expressed genes will be removed
  - ``MATRIX_MIN_CELLS``: genes present in fewer than this number of cells will be removed
  - ``MATRIX_MAX_MITO``: cells with more than this percentage of counts belonging to mitochondrial genes will be removed
  - ``MATRIX_NORM_COUNT``: normalize all cells to this number of total counts per cell

License and Copyright
---------------------

|copy| 2020-22 Oxford Nanopore Technologies Ltd.

.. |copy| unicode:: 0xA9 .. copyright sign

Sockeye is distributed under the terms of the Oxford Nanopore
Technologies, Ltd.  Public License, v. 1.0.  If a copy of the License
was not distributed with this file, You can obtain one at
http://nanoporetech.com

Research Release
----------------

Research releases are provided as technology demonstrators to provide early access to features or stimulate Community development of tools. Support for this software will be minimal and is only provided directly by the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the developers may have limited resource for support of this software. Research releases may be unstable and subject to rapid iteration by Oxford Nanopore Technologies.
