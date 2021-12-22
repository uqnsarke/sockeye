import argparse
import logging
import os
import subprocess
import sys

import pandas as pd
import pysam
import seaborn as sns
from matplotlib import pyplot as plt
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Aligned BAM file with gene, barcode, and UMI \
        tags",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output plot file with UMI saturation curve [output.png]",
        type=str,
        default="output.png",
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def run_subprocess(cmd):
    """
    Run OS command and return stdout & stderr
    """
    p = subprocess.Popen(
        cmd,
        shell=True,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )
    stdout, stderr = p.communicate()
    return str(stdout), str(stderr)


def plot_saturation_curves(res, args):
    """
    Output a single file with two subplots:
    1. Median number of unique genes per cell
    2. Median number of unique UMIs per cell
    """
    fig = plt.figure(figsize=[10, 5])

    ax1 = fig.add_subplot(1, 2, 1)
    sns.pointplot(x="reads_pc", y="genes_pc", data=res, ax=ax1, color="orange")

    ax1.set_xlabel("Median reads per cell")
    ax1.set_ylabel("Genes per cell")
    ax1.set_title("Gene saturation")

    # ax.xaxis.set_major_locator(MultipleLocator(2000))
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax.xaxis.set_minor_locator(MultipleLocator(500))

    # ax.yaxis.set_major_locator(MultipleLocator(2000))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax.yaxis.set_minor_locator(MultipleLocator(500))
    #
    # ax.set_xlim([-5,20000])
    # ax.set_ylim([-5,20000])
    # custom_legend(ax)
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax1.grid()

    ax2 = fig.add_subplot(1, 2, 2)
    sns.pointplot(x="reads_pc", y="umis_pc", data=res, ax=ax2, color="purple")

    ax2.set_xlabel("Median reads per cell")
    ax2.set_ylabel("UMIs per cell")
    ax2.set_title("UMI saturation")

    # ax.xaxis.set_major_locator(MultipleLocator(2000))
    # ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax.xaxis.set_minor_locator(MultipleLocator(500))

    # ax.yaxis.set_major_locator(MultipleLocator(2000))
    # ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax.yaxis.set_minor_locator(MultipleLocator(500))
    #
    # ax.set_xlim([-5,20000])
    # ax.set_ylim([-5,20000])
    # custom_legend(ax)
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax2.grid()

    fig.tight_layout()
    fig.savefig(args.output)


def load_bam_entries(n_aligns, args):
    """
    Load read_id, gene, corrected barcode, and corrected UMI values from BAM
    """
    bam = pysam.AlignmentFile(args.bam, "rb")

    records = []
    for i, align in enumerate(tqdm(bam.fetch(), total=n_aligns)):
        read_id = align.query_name

        # Annotated gene name = GN:Z
        gene = align.get_tag("GN")
        # Corrected cell barcode = CB:Z
        bc = align.get_tag("CB")
        # Corrected UMI = UB:Z
        umi = align.get_tag("UB")

        if gene != "NA":
            records.append((read_id, gene, bc, umi))

    bam.close()

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "barcode", "umi"]
    )

    n_barcodes = df["barcode"].unique().shape[0]
    return df, n_barcodes


def count_bam_entries(args):
    """
    Call `samtools view -c -F 260 <BAM>` to get number of entries
    """
    stdout, stderr = run_subprocess(f"samtools view -c -F 260 {args.bam}")
    return int(stdout.strip())


def downsample_reads(df, n_aligns):
    """
    Downsample dataframe of BAM entries and tabulate genes and UMIs per cell
    """
    fractions = [
        0.01,
        0.02,
        0.03,
        0.04,
        0.05,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]

    records = []
    for fraction in tqdm(fractions, total=len(fractions)):
        df_ = df.sample(frac=fraction)
        downsamp_reads = df_.shape[0]
        # Get the unique number of reads, genes and UMIs per cell barcode
        reads_per_cell = df_.groupby("barcode")["read_id"].nunique().median()
        genes_per_cell = df_.groupby("barcode")["gene"].nunique().median()
        umis_per_cell = df_.groupby("barcode")["umi"].nunique().median()
        records.append(
            (fraction, downsamp_reads, reads_per_cell, genes_per_cell, umis_per_cell)
        )

    res = pd.DataFrame.from_records(
        records,
        columns=["downsamp_frac", "downsamp_reads", "reads_pc", "genes_pc", "umis_pc"],
    )
    return res


def main(args):
    init_logger(args)

    logger.info(f"Counting alignments in {args.bam}")
    n_aligns = count_bam_entries(args)

    logger.info(f"Reading read_ids, genes, barcodes, and UMIs from {args.bam}")
    df, n_barcodes = load_bam_entries(n_aligns, args)

    logger.info("Downsampling BAM entries for per-cell gene and UMI saturation curves")
    res = downsample_reads(df, n_aligns)

    logger.info(f"Plotting saturation curves in {args.output}")
    plot_saturation_curves(res, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
