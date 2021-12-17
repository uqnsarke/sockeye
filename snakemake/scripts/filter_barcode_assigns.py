import argparse
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "table", help="FASTA file of putative barcode sequences", type=str
    )

    # Optional arguments
    parser.add_argument(
        "--max_ed",
        help="Max edit distance between putative barcode \
                        and the matching whitelist barcode [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--min_ed_diff",
        help="Min difference in edit distance between the \
                        (1) putative barcode vs top hit and (2) putative \
                        barcode vs runner-up hit [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--output_bc",
        help="Output file of filtered barcode assignments \
                        [bc_assigns.filtered.tsv]",
        type=str,
        default="bc_assigns.filtered.tsv",
    )

    parser.add_argument(
        "--output_counts",
        help="Output file of filtered barcode assignment read \
                        counts per barcode [bc_assigns.filtered.counts.tsv]",
        type=str,
        default="bc_assigns.filtered.counts.tsv",
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
    # logger.info(f'Generated {len(df_bcs)} barcodes of length {n}')
    # logger.warning('!! Failures detected !!')


def main(args):
    init_logger(args)

    logger.info(f"Reading demux table {args.table}")
    df = pd.read_csv(args.table, sep="\t")

    # Output file of filtered barcode assignments
    df = df[
        (df["match_ed"] <= args.max_ed)
        & (df["match_runner_up_diff"] >= args.min_ed_diff)
    ]

    df = df.rename(
        columns={"found_bc": "bc_uncorr", "match_bc": "bc_corr", "qscore": "bc_qv"}
    )
    logger.info(f"Writing reads per barcode to {args.output_bc}")
    df[["read_id", "bc_corr", "bc_uncorr", "bc_qv"]].to_csv(
        args.output_bc, sep="\t", index=False
    )

    # Output file of read counts per barcode
    bc_read_counts = (
        df.groupby("bc_corr")[["read_id"]]
        .agg("count")
        .reset_index()
        .rename(columns={"bc_corr": "barcode", "read_id": "reads"})
        .sort_values("reads", ascending=False)
    )
    logger.info(f"Writing reads per barcode to {args.output_counts}")
    bc_read_counts.to_csv(args.output_counts, sep="\t", index=False)


if __name__ == "__main__":
    args = parse_args()

    main(args)
