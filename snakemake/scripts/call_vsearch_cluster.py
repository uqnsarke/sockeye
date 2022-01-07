import argparse
import logging
import os
import subprocess
import sys

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    # parser.add_argument(
    #     "bc_fasta",
    #     help="FASTA containing all extracted barcodes to \
    #                     cluster",
    #     type=str,
    # )

    parser.add_argument(
        "hq_bc_fasta",
        help="FASTA containing all high-quality (i.e. filtered \
                        against the barcode superlist) extracted barcodes to \
                        cluster",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "-o",
        "--output",
        help="Output FASTA file containing cluster details and \
                        consensus sequences [output.vsearch_clustering.fasta]",
        default="output.vsearch_clustering.fasta",
        type=str,
    )

    parser.add_argument(
        "-t",
        "--threads",
        help="Threads to use for VSEARCH clustering [4]",
        default=4,
        type=int,
    )

    parser.add_argument(
        "-b",
        "--bc_len",
        help="Length of expected cell barcodes [16]",
        default=16,
        type=int,
    )

    parser.add_argument(
        "-i",
        "--min_id",
        help="Minimum sequence identity for clusters [0.94]",
        default=0.94,
        type=float,
    )

    parser.add_argument(
        "--verbosity",
        help="Logging level: <=2 logs info, <=3 logs warnings",
        default=2,
        type=int,
    )

    args = parser.parse_args()

    return args


def run_subprocess(cmd, pipe=True):
    """
    Run OS command and return stdout & stderr
    """
    if pipe:
        p = subprocess.Popen(
            cmd,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
        )
    else:
        p = subprocess.Popen(cmd, shell=True, universal_newlines=True)
    stdout, stderr = p.communicate()
    return str(stdout), str(stderr)


def check_vsearch():
    stdout, stderr = run_subprocess("vsearch --quiet -h")
    if stderr.find("vsearch: command not found") > -1:
        logging.error("Could not load find VSEARCH -- check installation")
        sys.exit(1)


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def main(args):
    init_logger(args)
    check_vsearch()

    # if os.path.getsize(args.hq_bc_fasta) != 0:
    #     # We have high-quality barcodes to cluster
    #     input_fasta = args.hq_bc_fasta
    # else:
    #     # Just cluster all found barcodes
    input_fasta = args.hq_bc_fasta

    logger.info("Running VSEARCH clustering")
    vsearch_cmd = f"vsearch --clusterout_id --consout {args.output} \
                  --minseqlength {args.bc_len} \
                  --maxseqlength {args.bc_len} \
                  --threads {args.threads} \
                  -id {args.min_id} \
                  --clusterout_sort \
                  --qmask none \
                  --cluster_fast {input_fasta}"

    stdout, stderr = run_subprocess(vsearch_cmd, pipe=False)
    logger.info(f"Wrote VSEARCH clusters to {args.output}")


if __name__ == "__main__":
    args = parse_args()

    main(args)
