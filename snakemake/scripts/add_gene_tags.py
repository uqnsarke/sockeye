import argparse
import gzip
import logging
import os
import pathlib
import subprocess
import sys

import numpy as np
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("bam", help="Sorted BAM file", type=str)

    parser.add_argument(
        "fc",
        help="FeatureCounts read/gene assignments file. \
        IMPORTANT: featureCounts must have been run on the same sorted BAM file \
        that serves as the other input argument.",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output BAM file containing aligned reads with gene name tags (GN) \
        [gene.sorted.bam]",
        type=str,
        default="gene.sorted.bam",
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


def get_bam_info(bam):
    """
    Use `samtools idxstat` to get number of alignments and names of all contigs
    in the reference.

    :param bam: Path to sorted BAM file
    :type bame: str
    :return: Sum of all alignments in the BAM index file and list of all chroms
    :rtype: int,list
    """
    bam = pysam.AlignmentFile(bam, "rb")
    stats = bam.get_index_statistics()
    n_aligns = int(np.sum([contig.mapped for contig in stats]))
    chroms = [contig.contig for contig in stats]
    bam.close()
    return n_aligns, chroms


def process_bam_entries(args):
    """
    Iterate simultaneously throught the BAM and the featureCounts TSV file of
    read/gene assignments. These files must be sorted identically, otherwise it
    will throw an exception. The gene names found in the featureCounts TSV are
    added to the GN tag for each alignment in the output BAM.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
    n_reads, chroms = get_bam_info(args.bam)

    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_out = pysam.AlignmentFile(args.output, "wb", template=bam)

    logger.info(f"Adding gene tags (GN) to {args.bam}")
    with open(args.fc, "r") as fc:
        for align in tqdm(bam.fetch(), total=n_reads):
            line = fc.readline().strip()
            fc_id = line.split("\t")[0]
            fc_gene = line.split("\t")[3]

            assert fc_id == align.query_name, "BAM and featureCounts reads not ordered"

            # Annotated gene name = GN:Z
            align.set_tag("GN", fc_gene, value_type="Z")

            bam_out.write(align)

    bam.close()
    bam_out.close()


def main(args):
    init_logger(args)

    process_bam_entries(args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
