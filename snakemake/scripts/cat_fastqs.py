import argparse
import collections
import logging
import multiprocessing
import os
import pathlib
import shutil

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "fofn",
        help="File of filenames containing paths to all input FASTQs",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "-t",
        "--threads",
        help="Threads to use. This is also the number of \
        chunks that will be output [4]",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--output_dir",
        help="Output directory for chunked FASTQ files [./chunked_fastqs]",
        type=str,
        default="./chunked_fastqs",
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


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def load_fofn(fn):
    """ """
    input_fastqs = []
    for line in open(fn, "r"):
        path = line.strip()
        input_fastqs.append(path)
    return input_fastqs


def cat_files(chunk_fns, i, args):
    """ """
    out_fn = os.path.join(args.output_dir, f"proc.{i}.fastq.gz")
    with open(out_fn, "wb") as outfile:
        for fname in chunk_fns:
            with open(fname, "rb") as infile:
                outfile.write(infile.read())


def main(args):
    input_fastqs = load_fofn(args.fofn)

    # Create output directory
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.mkdir(args.output_dir)

    chunk_size = int((len(input_fastqs) / args.threads)) + 1
    for i, chunk_fns in enumerate(chunks(input_fastqs, chunk_size)):
        cat_files(chunk_fns, i + 1, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
