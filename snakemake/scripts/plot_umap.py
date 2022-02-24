import argparse
import logging
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "umap",
        help="File containing the 2D UMAP projection of cell barcodes.",
        type=str,
    )

    parser.add_argument(
        "full_matrix",
        help="File containing the full expression matrix that was used for \
        the UMAP projection.",
        type=str,
        default=None,
    )

    # Optional arguments
    parser.add_argument(
        "--output_prefix",
        help="Output file prefix to use when saving annotated UMAP plots \
            (<output_prefix>.FEATURE.png) [umap]",
        type=str,
        default="umap",
    )

    parser.add_argument(
        "-g",
        "--genes",
        help="List of genes to annotate in UMAP plots, separated by commas \
        (e.g. --genes=CD19,CD38,CD27) [None]",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-t",
        "--target_cells",
        help="List of cells to highlight in UMAP plots [None]",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-s", "--size", help="Size of markers [15]", type=int, default=15
    )

    parser.add_argument(
        "-a", "--alpha", help="Transpancy of markers [0.4]", type=float, default=0.4
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


def remove_top_right_axes(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="both", direction="in")
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def scatterplot(df, values, args):
    """ """
    fig = plt.figure(figsize=[8, 8])
    ax = fig.add_axes([0.08, 0.08, 0.85, 0.85])

    if values.name == "total":
        cmap = cm.jet
    else:
        cmap = cm.cool

    plot = ax.scatter(
        df["D1"],
        df["D2"],
        s=args.size,
        edgecolor="None",
        c=values,
        cmap=cmap,
        alpha=args.alpha,
    )

    remove_top_right_axes(ax)

    plt.xlim([df["D1"].min() - 1, df["D1"].max() + 1])
    plt.ylim([df["D2"].min() - 1, df["D2"].max() + 1])

    n_cells = df.shape[0]
    if values.name == "highlight":
        title = f"{n_cells} cells: highlighted cells from {args.target_cells}"
    else:
        title = f"{n_cells} cells: {values.name}"
        plt.colorbar(plot)

    ax.set_title(title)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")

    fn = f"{args.output_prefix}.{values.name}.png"
    plt.savefig(fn)


def get_expression(args):
    """ """
    df_f = (
        pd.read_csv(args.full_matrix, delimiter="\t")
        .rename(columns={"gene": "cell"})
        .set_index("cell")
    )
    df_f = df_f.transpose()

    if args.genes:
        genes = args.genes.split(",")

        # Make sure requested genes are in the matrix
        for g in genes:
            if g not in df_f.columns:
                logger.info(f"WARNING: gene {g} not found in expression matrix!")

    # Create annotation dataframe with requested features
    df_annot = pd.DataFrame()
    df_annot["total"] = np.exp(df_f).sum(axis=1) - 1
    for g in genes:
        df_annot[g] = df_f.loc[:, g]

    return df_annot


def main(args):
    init_logger(args)

    df = pd.read_csv(args.umap, delimiter="\t")

    df_annot = get_expression(args)

    for feature in df_annot.columns:
        logger.info(f"Plotting UMAP with {feature} annotation")

        values = df_annot.loc[:, feature]
        scatterplot(df, values, args)

    if args.target_cells:
        logger.info(f"Plotting UMAP with highlighted cells from {args.target_cells}")

        df_highlight = pd.read_csv(args.target_cells, header=None, names=["barcode"])
        df_highlight["barcode"] = df_highlight["barcode"].str.split(
            "-", n=0, expand=True
        )
        df_highlight["highlight"] = "red"
        df_highlight = df_highlight.set_index("barcode")
        df_highlight = pd.merge(df, df_highlight, on="barcode", how="left").fillna(
            "lightgray"
        )
        values = df_highlight.loc[:, "highlight"]
        scatterplot(df, values, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
