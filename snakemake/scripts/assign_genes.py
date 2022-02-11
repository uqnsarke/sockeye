import argparse
import collections
import logging
import os
import pathlib
import shutil
import sys

import bioframe as bf
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bed",
        help="BED file of alignments intervals",
        type=str,
    )

    parser.add_argument(
        "gtf",
        help="GTF file of gene annotations",
        type=str,
    )

    # parser.add_argument(
    #     "unmapped_reads",
    #     help="List of unmapped reads from the BAM",
    #     type=str,
    # )

    # Optional arguments
    parser.add_argument(
        "-q",
        "--mapq",
        help="Minimum mapping quality to use for feature assignment [60]",
        type=int,
        default=60,
    )

    parser.add_argument(
        "--output",
        help="Output file [./read_annotations.tsv]",
        type=str,
        default="./read_annotations.tsv",
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


def load_gtf(args):
    """ """
    cols = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    df = pd.read_csv(args.gtf, sep="\t", comment="#", header=None, names=cols)
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "GTF file not loading as a valid dataframe!"

        # Restrict annotations to the gene level and isolate gene name from attributes
        df = df[df["feature"] == "gene"]
        df["attribute"] = df["attribute"].str.split(";", expand=True).iloc[:, 3]
        df["attribute"] = df["attribute"].str.split(" ", expand=True).iloc[:, -1]
        df["attribute"] = df["attribute"].str.replace('"', "")
    return df


def load_bed(args):
    """
    Maybe use bioframe.io.read_bigbed if these files are too big.
    """
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]
    df = pd.read_csv(args.bed, sep="\t", header=None, names=cols)
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "BED file not loading as a valid dataframe!"

    df["aln_len"] = df["end"] - df["start"]
    return df


def load_unmapped_reads(args):
    """ """
    cols = ["row_read", "chrom", "start"]
    df = pd.read_csv(
        args.unmapped_reads, sep="\t", header=None, usecols=["row_read"], names=cols
    )
    df[["index_bed", "read"]] = df["row_read"].str.split(":", expand=True)
    df = df.drop("row_read", axis=1)
    df["status"] = "Unassigned_unmapped"
    df["score"] = 0
    df["gene"] = "NA"
    # Switch from 1-index to 0-index
    df["index_bed"] = df["index_bed"].astype(int) - 1
    return df


def assign_status_low_mapq(df, args):
    """
    Assign all reads with ambiguous mapping as "Unassigned_mapq"
    """
    df.loc[df["score"] < args.mapq, "status"] = "Unassigned_mapq"
    df.loc[df["score"] < args.mapq, "gene"] = "NA"
    return df


def assign_status_ambiguous_overlap(df):
    """
    Assign all reads with equal overlap with multiple features as
    "Unassigned_ambiguous"
    """
    is_ambiguous = df[df["status"] == "Unknown"].duplicated(
        subset=["index_bed", "overlap_bp"], keep=False
    )
    ambiguous_idx = is_ambiguous.index[is_ambiguous]
    df.loc[ambiguous_idx, "status"] = "Unassigned_ambiguous"
    df.loc[ambiguous_idx, "gene"] = "NA"
    df = df.drop_duplicates(subset=["index_bed", "overlap_bp", "status"])

    return df


def assign_status_no_features(df):
    """
    Assign all reads without any features as "Unassigned_no_features"
    """
    df.loc[df["gene"] == 0, "status"] = "Unassigned_no_features"
    df.loc[df["gene"] == 0, "gene"] = "NA"
    return df


def assign_status_unmapped(df):
    """ """

    return df


def find_largest_overlap(df):
    """
    Find the largest overlap when multiple overlaps are present.
    """
    # Find the indices of the largest overlap for each alignment entry
    max_ovlp_idx = df.groupby(["index_bed"])["overlap_bp"].idxmax().sort_values().values

    # Only keep the entry with the largest overlap
    df = df.loc[max_ovlp_idx, :]
    unknown_idx = [i for i in max_ovlp_idx if df.loc[i, "status"] == "Unknown"]
    df.loc[unknown_idx, "status"] = "Assigned"
    return df


def get_overlaps(bed, gtf):
    """ """
    df = bf.overlap(
        bed,
        gtf,
        how="left",
        suffixes=("_bed", "_gtf"),
        return_overlap=True,
        return_index=True,
    )

    # !!!Now keeping index_bed to keep track of the input order. We'll want to
    # get the largest overlap for each index_bed value, which will properly
    # handle the issue of supplementary alignments in the BAM.

    # Only keep relevant columns
    df = df[
        [
            "index_bed",
            "name_bed",
            "chrom_bed",
            "score_bed",
            "strand_gtf",
            "attribute_gtf",
            "overlap_start",
            "overlap_end",
        ]
    ].fillna(0)

    df = df.rename(
        columns={
            "name_bed": "read",
            "chrom_bed": "chrom",
            "score_bed": "score",
            "strand_gtf": "strand",
            "attribute_gtf": "gene",
        }
    )
    df["score"] = df["score"].astype(int)
    df["overlap_bp"] = df["overlap_end"] - df["overlap_start"]

    return df


def main(args):
    pd.set_option("display.max_rows", None)
    gtf = load_gtf(args)
    bed = load_bed(args)

    if (bed.shape[0] > 0) & (gtf.shape[0] > 0):
        # The bed file has alignments and the chromosome has annotations,
        # so process the overlaps
        df = get_overlaps(bed, gtf)

        df["status"] = "Unknown"
        df = assign_status_low_mapq(df, args)
        df = assign_status_ambiguous_overlap(df)
        df = assign_status_no_features(df)
        # df = assign_status_unmapped(df)
        df = find_largest_overlap(df)

        df = df[["read", "status", "score", "gene", "index_bed"]]

        # Now add the unmapped reads from the end of the BAM
        # unmapped = load_unmapped_reads(args)
        # df = pd.concat([df, unmapped], axis=0)

        df = df.reset_index(drop=True)
        # print(df.head(500))
        # print(df.tail(500))
        print(df.groupby("status")["index_bed"].nunique())
        df = df.drop(["index_bed"], axis=1)

        df.to_csv(args.output, sep="\t", index=False, header=False)

    else:
        # The bed file contained no alignments or the chromosome does not
        # have any annotations, so output empty file
        f = open(args.output, "w")
        f.close()


if __name__ == "__main__":
    args = parse_args()

    main(args)
