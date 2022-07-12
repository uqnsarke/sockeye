import pathlib
import os
import sys
import gzip
import shutil
import requests
import pandas as pd
from glob import glob
from itertools import product
from snakemake.utils import validate, min_version


##### load config and point to scripts #####
# configfile: "config/config.yml"
if "--configfile" in sys.argv:
    i = sys.argv.index("--configfile")
    config_path = sys.argv[i + 1]

    configfile: config_path


else:
    config_path = "config/config.yml"

    configfile: config_path


SCRIPT_DIR = srcdir("scripts")


###############################################
# Download 10X barcode longlists if necessary #
###############################################
urls = {
    "3M-february-2018.txt.gz": "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz",
    "737K-august-2016.txt": "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt",
}

ll_dir = config["BC_LONGLIST_DIR"]

if not os.path.exists(ll_dir):
    os.mkdir(ll_dir)

for fn in urls.keys():
    local_fn = os.path.join(ll_dir, fn)

    if not os.path.exists(local_fn):
        print(f"Downloading {fn} to {ll_dir}")
        with open(local_fn, "wb") as f:
            chunkSize = 1024
            r = requests.get(urls[fn], stream=True)
            with open(local_fn, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunkSize):
                    if chunk:  # filter out keep-alive new chunks
                        f.write(chunk)

############################
# Validate kit_configs.csv #
############################
kit_df = pd.read_csv("./config/kit_configs.csv", sep=",", comment="#")


#######################
# Validate config.yml #
#######################
if not config.get("SAMPLE_SHEET"):
    raise Exception(f"Please define SAMPLE_SHEET in the {config_path}")
elif not os.path.exists(config["SAMPLE_SHEET"]):
    raise Exception(f"Path specified for SAMPLE_SHEET in {config_path} not found!")
sample_sheet = pd.read_csv(
    config.get("SAMPLE_SHEET", "./config/samples.csv"), sep=",", comment="#"
).set_index("run_id", drop=True)
RUN_IDS = sample_sheet.index

if not config.get("BC_LONGLIST_DIR"):
    raise Exception(f"Please define BC_LONGLIST_DIR in the {config_path}")

if not config.get("REF_GENOME_DIR"):
    raise Exception(f"Please define REF_GENOME_DIR in the {config_path}")
elif not os.path.exists(config["REF_GENOME_DIR"]):
    raise Exception(f"Path specified for REF_GENOME_DIR in {config_path} not found!")
REF_GENOME_DIR = pathlib.Path(config["REF_GENOME_DIR"])
REF_GENOME_FASTA = REF_GENOME_DIR / "fasta/genome.fa"
REF_GENES_GTF = REF_GENOME_DIR / "genes/genes.gtf"

PLOT_GENES = list(map(lambda x: x.upper(), config.get("UMAP_PLOT_GENES", None).split(",")))

BC_LONGLIST_DIR = Path(config["BC_LONGLIST_DIR"])
BC_LONGLIST_3PRIME = BC_LONGLIST_DIR / "3M-february-2018.txt.gz"
BC_LONGLIST_5PRIME = BC_LONGLIST_DIR / "737K-august-2016.txt"

##### Set output location #####
OUTPUT_BASE = pathlib.Path(config.get("OUTPUT_BASE", "./output"))


###################
# READ ADAPTER SCAN FILES
###################
INGEST_DIR = OUTPUT_BASE / "{run_id}" / "ingest"
FOFN = str(INGEST_DIR / "fofn.txt")
CHUNKED_FASTQ_DIR = INGEST_DIR / "chunked_fastqs"
CHUNKED_FASTQS = str(CHUNKED_FASTQ_DIR / "{batch_id}.fastq")
ADAPTERS_DIR = OUTPUT_BASE / "{run_id}" / "adapters"
READ_CONFIG_CHUNKED = str(ADAPTERS_DIR / "{batch_id}.configs.tsv")
STRANDED_FQ_CHUNKED = str(ADAPTERS_DIR / "{batch_id}.stranded.fastq")
READ_CONFIG = str(ADAPTERS_DIR / "configs.tsv")
STRANDED_FQ = str(ADAPTERS_DIR / "reads.stranded.fastq.gz")
CONFIG_STATS = str(ADAPTERS_DIR / "configs.stats.json")


###################
# REF ALIGN FILES
###################
REF_GENES_DIR = REF_GENES_GTF.parents[0]
REF_GENES_BED = str(REF_GENES_DIR / f'{REF_GENES_GTF.name.replace(".gtf", ".bed")}')
ALIGN_DIR = OUTPUT_BASE / "{run_id}" / "align"
REF_CHROM_SIZES = str(ALIGN_DIR / "chrom_sizes.tsv")
SAM_TMP = str(ALIGN_DIR / "tmp.sam")
BAM_UNSORT_TMP = str(ALIGN_DIR / "tmp.unsort.bam")
BAM_SORT = str(ALIGN_DIR / "sorted.bam")
BAM_SORT_BAI = str(ALIGN_DIR / "sorted.bam.bai")


###################
# BARCODE + UMI DEMUX FILES
###################
DEMUX_DIR = OUTPUT_BASE / "{run_id}" / "demux"
# All chrom files
BAM_BC_UNCORR_TMP = str(DEMUX_DIR / "bc_extract.tmp.sorted.bam")
BAM_BC_UNCORR = str(DEMUX_DIR / "bc_extract.sorted.bam")
BAM_BC_UNCORR_BAI = str(DEMUX_DIR / "bc_extract.sorted.bam.bai")
BARCODE_COUNTS = str(DEMUX_DIR / "uncorrected_bc_counts.tsv")
BARCODE_WHITELIST = str(DEMUX_DIR / "whitelist.tsv")
BARCODE_KNEEPLOT = str(DEMUX_DIR / "kneeplot.png")

# Chrom-specific files
SPLIT_DIR = DEMUX_DIR / "chroms"
CHROM_BAM_BC_UNCORR = str(SPLIT_DIR / "{chrom}.sorted.bam")
CHROM_BAM_BC_UNCORR_BAI = str(SPLIT_DIR / "{chrom}.sorted.bam.bai")
CHROM_BAM_BC_TMP = str(SPLIT_DIR / "{chrom}.bc_assign.tmp.bam")
CHROM_BAM_BC = str(SPLIT_DIR / "{chrom}.bc_assign.bam")
CHROM_BAM_BC_BAI = str(SPLIT_DIR / "{chrom}.bc_assign.bam.bai")
CHROM_ASSIGNED_BARCODE_COUNTS = str(SPLIT_DIR / "{chrom}.bc_assign_counts.tsv")
CHROM_BED_BC = str(SPLIT_DIR / "{chrom}.bc_assign.bed")
SPLIT_ANNOT_DIR = DEMUX_DIR / "refs"
CHROM_GTF = str(SPLIT_ANNOT_DIR / "{chrom}.genes.gtf")
CHROM_TSV_GENE_ASSIGNS = str(SPLIT_DIR / "{chrom}.read.gene_assigns.tsv")
CHROM_BAM_BC_GENE_TMP = str(SPLIT_DIR / "{chrom}.bc_assign.gene.tmp.bam")
CHROM_BAM_BC_GENE = str(SPLIT_DIR / "{chrom}.bc_assign.gene.bam")
CHROM_BAM_BC_GENE_BAI = str(SPLIT_DIR / "{chrom}.bc_assign.gene.bam.bai")
CHROM_BAM_FULLY_TAGGED_TMP = str(SPLIT_DIR / "{chrom}.tagged.tmp.bam")
CHROM_BAM_FULLY_TAGGED = str(SPLIT_DIR / "{chrom}.tagged.bam")
CHROM_BAM_FULLY_TAGGED_BAI = str(SPLIT_DIR / "{chrom}.tagged.bam.bai")

# Merged files
BAM_FULLY_TAGGED = str(DEMUX_DIR / "tagged.sorted.bam")
BAM_FULLY_TAGGED_BAI = str(DEMUX_DIR / "tagged.sorted.bam.bai")
ASSIGNED_BARCODE_COUNTS = str(DEMUX_DIR / "assigned_bc_counts.tsv")
CELL_UMI_GENE_TSV = str(DEMUX_DIR / "cell_umi_gene.tsv")


###################
# GENE EXPRESSION FILES
###################
MATRIX_DIR = OUTPUT_BASE / "{run_id}" / "matrix"
MATRIX_COUNTS_TSV = str(MATRIX_DIR / "gene_expression.counts.tsv")
MATRIX_PROCESSED_TSV = str(MATRIX_DIR / "gene_expression.processed.tsv")
MATRIX_UMAP_TSV = str(MATRIX_DIR / "gene_expression.umap.tsv")
MATRIX_UMAP_PLOT_GENE = str(MATRIX_DIR / "umap.gene.{plot_gene}.png")
MATRIX_UMAP_PLOT_MITO = str(MATRIX_DIR / "umap.mitochondrial.png")
MATRIX_UMAP_PLOT_TOTAL = str(MATRIX_DIR / "umap.total.png")


###################
# LIBRARY SATURATION FILES
###################
SAT_DIR = OUTPUT_BASE / "{run_id}" / "saturation"
SAT_PLOT = str(SAT_DIR / "saturation_curves.png")


##### include rules #####
include: "rules/stranding.smk"
include: "rules/align.smk"
include: "rules/process_bams.smk"


wildcard_constraints:
    chrom=r"\w+(\.\d+)?",


##### target rules #####
rule all:
    input:
        expand(CONFIG_STATS, run_id=RUN_IDS),
        expand(MATRIX_UMAP_TSV, run_id=RUN_IDS),
        expand(SAT_PLOT, run_id=RUN_IDS),
        expand(MATRIX_UMAP_PLOT_TOTAL, run_id=RUN_IDS),
        expand(MATRIX_UMAP_PLOT_MITO, run_id=RUN_IDS),
        expand(MATRIX_UMAP_PLOT_GENE, run_id=RUN_IDS, plot_gene=PLOT_GENES),
