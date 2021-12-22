import argparse
import gzip
import logging
import os
import pathlib
import subprocess
import sys

import pysam
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("bam", help="Aligned BAM file", type=str)

    parser.add_argument(
        "fastq", help="FASTQ file of reads with tag values in read headers", type=str
    )

    parser.add_argument("fc", help="FeatureCounts read/gene assignments file", type=str)

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output barcode/umi tagged BAM file name \
                        [sorted.bc_umi.bam]",
        type=str,
        default="sorted.bc_umi.bam",
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


class Read:
    def __init__(self, read_id):
        self.read_id = read_id

    def add_barcode_info(self, read_name):
        # Read key=value pairs from read header into a dictionary
        kv = dict([kv.split("=") for kv in read_name.split(" ") if kv.find("=") > -1])
        self.bc_uncorr = kv["bc_uncorr"]
        self.bc_corr = kv["bc_corr"]
        self.bc_qv = kv["bc_qv"]
        self.umi_uncorr = kv["umi_uncorr"]
        # self.umi_corr = ###
        self.umi_qv = kv["umi_qv"]

    def add_gene_info(self, gene):
        self.gene = gene


def count_fastq_reads(fastq):
    logger.info("Counting reads")
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in tqdm(f, unit_scale=0.25, unit=" reads"):
            number_lines += 1
    return number_lines / 4


def open_fastq(fastq):
    if fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq, "r+")
    return f


def read_fastq(args):
    """
    Read cell barcode and UMI values from FASTQ read headers
    """
    n_reads = int(count_fastq_reads(args.fastq))
    f = open_fastq(args.fastq)
    record_iter = FastqGeneralIterator(f)
    read_tags = {}
    logging.info("Reading barcode / UMI tag info from FASTQ read headers")
    for (read_name, seq, qv) in tqdm(record_iter, total=n_reads):
        read_id = read_name.split(" ")[0]
        r = Read(read_id)
        r.add_barcode_info(read_name)
        read_tags[read_id] = r

    return read_tags


def read_fc(args, read_tags):
    """
    Read aligned gene names from featureCounts output file
    """
    for i, line in enumerate(open(args.fc, "r")):
        read_id = line.strip().split("\t")[0]
        gene = line.strip().split("\t")[3]

        # Only some of the featureCounts reads will have barcodes / UMIs
        if read_tags.get(read_id):
            read_tags[read_id].add_gene_info(gene)

    return read_tags


def process_bam_entries(read_tags, args):
    """ """
    logger.info(f"Counting alignments in {args.bam}")
    stdout, stderr = run_subprocess(f"samtools view -c -F 260 {args.bam}")
    n_aligns = int(stdout.strip())

    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_tagged = pysam.AlignmentFile(args.output, "wb", template=bam)

    logger.info(f"Adding barcode / UMI tags to {args.bam}")
    for i, align in enumerate(tqdm(bam.fetch(), total=n_aligns)):
        read_id = align.query_name

        # If alignment has assigned cell barcode / UMI
        if read_tags.get(read_id):
            # Corrected cell barcode = CB:Z
            align.set_tag("CB", read_tags[read_id].bc_corr, value_type="Z")
            # Uncorrected cell barcode = CR:Z
            align.set_tag("CR", read_tags[read_id].bc_uncorr, value_type="Z")
            # Cell barcode quality score = CY:Z
            align.set_tag("CY", read_tags[read_id].bc_qv, value_type="Z")
            # Uncorrected UMI = UR:Z
            align.set_tag("UR", read_tags[read_id].umi_uncorr, value_type="Z")
            # UMI quality score = UY:Z
            align.set_tag("UY", read_tags[read_id].umi_qv, value_type="Z")
            # Annotated gene name = GN:Z
            align.set_tag("GN", read_tags[read_id].gene, value_type="Z")
        else:
            # Corrected cell barcode = CB:Z
            align.set_tag("CB", "XXXXXXXXX", value_type="Z")
            # Uncorrected cell barcode = CR:Z
            align.set_tag("CR", "XXXXXXXXX", value_type="Z")
            # Cell barcode quality score = CY:Z
            align.set_tag("CY", "0.0", value_type="Z")
            # Uncorrected UMI = UR:Z
            align.set_tag("UR", "XXXXXXXXX", value_type="Z")
            # UMI quality score = UY:Z
            align.set_tag("UY", "0.0", value_type="Z")
            # Annotated gene name = GN:Z
            align.set_tag("GN", "NA", value_type="Z")

        bam_tagged.write(align)

    bam.close()
    bam_tagged.close()


def main(args):
    init_logger(args)

    read_tags = read_fastq(args)

    read_tags = read_fc(args, read_tags)

    process_bam_entries(read_tags, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
