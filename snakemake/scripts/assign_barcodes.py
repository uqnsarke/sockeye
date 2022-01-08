import argparse
import gzip
import logging
import math
import multiprocessing
import os
import pathlib
import re
import shutil
import tempfile

import editdistance as ed
import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "fastq",
        help="FASTQ file of reads with putative barcode sequences in header",
        type=str,
    )

    parser.add_argument(
        "whitelist", help="File containing list of expected cell barcodes", type=str
    )

    # Optional arguments
    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "-k", help="Kmer size to use for whitelist filtering [5]", type=int, default=5
    )

    parser.add_argument(
        "--output_all",
        help="Output file containing all reads with unfiltered barcode assignments \
        [reads.bc_corr.fastq]",
        type=str,
        default="reads.bc_corr.fastq",
    )

    parser.add_argument(
        "--output_filtered",
        help="Output file containing only reads with high-confidence barcode \
        assignments [reads.bc_corr.filt.fastq]",
        type=str,
        default="reads.bc_corr.filt.fastq",
    )

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
        "--barcode_length", help="Cell barcode length [16]", type=int, default=16
    )

    parser.add_argument(
        "-b",
        "--batch_size",
        help="Number of reads to load at once [20000]",
        type=int,
        default=20000,
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    # Create temp dir and add that to the args object
    p = pathlib.Path(args.output_all)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)
    # progress = tqdm if logging_level <= 20 else None
    # logger.info(f'Generated {len(df_bcs)} barcodes of length {n}')
    # logger.warning('!! Failures detected !!')


class Barcode:
    def __init__(self, read_id, found_bc):
        """
        Define basic characteristic of the barcode that we
        want to classify
        """
        self.read_id = read_id
        self.found_bc = found_bc

    def calc_ed_with_whitelist(self, whitelist):
        """
        Find minimum and runner-up barcode edit distance by
        iterating through the whitelist of expected barcodes
        """
        self.match_bc = "X" * len(self.found_bc)
        self.match_ed = len(self.found_bc)
        self.runner_up_bc = "X" * len(self.found_bc)
        self.runner_up_ed = len(self.found_bc)
        for wl_bc in whitelist:
            d = ed.eval(self.found_bc, wl_bc)
            if d < self.match_ed:
                self.runner_up_ed = self.match_ed
                self.runner_up_bc = self.match_bc
                self.match_ed = d
                self.match_bc = wl_bc
        self.match_runner_up_diff = self.runner_up_ed - self.match_ed


def get_val_from_read_header(REGEX, header):
    """
    Use regular expression to extract value from read header
    """
    m = re.search(REGEX, header)
    return m.group(1)


def calc_ed(tup):
    """
    Filter the whitelist by kmers and find least edit distance match
    """
    read = tup[0]
    kmer_to_bc_index = tup[1]
    whitelist = tup[2]
    args = tup[3]

    REGEX_bc_uncorr = r"bc_uncorr=([ACGT]{{{}}})".format(args.barcode_length)
    barcode = get_val_from_read_header(REGEX_bc_uncorr, read.description)

    REGEX_bc_qv = r"bc_qv=(\d+\.\d+)"
    bc_qv = get_val_from_read_header(REGEX_bc_qv, read.description)

    kmers = split_barcode_into_kmers(barcode, args.k)
    filt_whitelist = filter_whitelist_by_kmers(whitelist, kmers, kmer_to_bc_index)

    # Instantiate barcode object
    read_bc = Barcode(read.id, barcode)
    # Calc edit distances and store min & runner up
    read_bc.calc_ed_with_whitelist(filt_whitelist)

    # Construct the updated read header
    read.description = f"{read.id}"
    read.description += f" bc_uncorr={barcode}"
    read.description += f" bc_corr={read_bc.match_bc}"
    read.description += f" bc_qv={bc_qv}"
    read.description += f" ed={read_bc.match_ed}"
    read.description += f" diff={read_bc.match_runner_up_diff}"

    # Check barcode match edit distance and difference to runner-up edit distance
    if (read_bc.match_ed <= args.max_ed) and (
        read_bc.match_runner_up_diff >= args.min_ed_diff
    ):
        # This is a high-confidence barcode assignment
        read.features.append({"hi-conf": True})
    else:
        # This is NOT a high-confidence barcode assignment
        read.features.append({"hi-conf": False})

    return read


def launch_pool(procs, funct, args):
    p = multiprocessing.Pool(processes=procs)
    try:
        results = p.map(funct, args)
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


def filter_whitelist_by_kmers(wl, kmers, kmer_to_bc_index):
    """
    Given a list of whitelisted barcodes, return just the
    subset that contain any of the kmers contained in the
    query barcode.
    """
    # collect sets of indices that each kmer points to
    id_sets = [
        kmer_to_bc_index[kmer] for kmer in kmers if kmer in kmer_to_bc_index.keys()
    ]

    # retain all barcodes that have at least one kmer match with the query barcode
    all_filt_indices = list(set().union(*id_sets))
    filt_wl = [wl[i] for i in all_filt_indices]
    return filt_wl


def launch_alignment_pool(batch_reads, whitelist, kmer_to_bc_index, args):
    func_args = []

    for read in batch_reads:
        # Read barcode from read header
        func_args.append((read, kmer_to_bc_index, whitelist, args))

    reads = launch_pool(args.threads, calc_ed, func_args)
    return reads


def split_barcode_into_kmers(bc, k):
    kmers = []
    for i in range(0, len(bc) - k + 1):
        kmer = bc[i : i + k]
        kmers.append(kmer)
    return kmers


def load_whitelist(args):
    wl = []
    with open(args.whitelist) as file:
        for line in file:
            bc = line.strip().split("-")[0]
            wl.append(bc)

    wl.sort()
    kmer_to_bc_index = {}
    for index, bc in enumerate(wl):
        bc_kmers = split_barcode_into_kmers(bc, args.k)
        for bc_kmer in bc_kmers:
            if bc_kmer not in kmer_to_bc_index.keys():
                kmer_to_bc_index[bc_kmer] = set([index])
            else:
                kmer_to_bc_index[bc_kmer].add(index)
    return wl, kmer_to_bc_index


def open_fastq(fastq):
    if fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq, "r+")
    return f


def check_input_format(fastq):
    f = open_fastq(fastq)
    line = f.readline()
    if line[0] == "@":
        pass
    else:
        raise ("Unexpected file type! FASTQ only!")
    return


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def count_barcodes(fastq):
    logging.info("Counting putative barcodes")
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in tqdm(f, unit_scale=0.5, unit=" barcodes"):
            number_lines += 1
    return int(number_lines / 4)


def write_tmp_fastq(fastq_entries, args):
    tmp_read_fastq = tempfile.NamedTemporaryFile(
        prefix="tmp.read.bc_corr.", suffix=".fastq", dir=args.tempdir, delete=False
    )
    tmp_filt_read_fastq = tempfile.NamedTemporaryFile(
        prefix="tmp.read.bc_corr.filt.", suffix=".fastq", dir=args.tempdir, delete=False
    )
    fastq_filt_entries = [r for r in fastq_entries if r.features[0]["hi-conf"]]
    SeqIO.write(fastq_entries, tmp_read_fastq.name, "fastq")
    SeqIO.write(fastq_filt_entries, tmp_filt_read_fastq.name, "fastq")
    return tmp_read_fastq.name, tmp_filt_read_fastq.name


def main(args):
    init_logger(args)
    whitelist, kmer_to_bc_index = load_whitelist(args)
    check_input_format(args.fastq)
    n_barcodes = count_barcodes(args.fastq)
    n_batches = int(math.ceil(n_barcodes / args.batch_size))

    f = open_fastq(args.fastq)
    record_iter = SeqIO.parse(f, "fastq")

    logger.info(
        "Processing {n} barcodes in {b} batches".format(n=n_barcodes, b=n_batches)
    )
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)
    tmp_fastqs = []
    tmp_filt_fastqs = []
    for i, batch_reads in enumerate(
        tqdm(batch_iterator(record_iter, args.batch_size), total=n_batches)
    ):

        fastq_entries = launch_alignment_pool(
            batch_reads, whitelist, kmer_to_bc_index, args
        )
        tmp_fastq, tmp_filt_fastq = write_tmp_fastq(fastq_entries, args)
        tmp_fastqs.append(tmp_fastq)
        tmp_filt_fastqs.append(tmp_filt_fastq)

    logger.info(f"Writing all reads with bc_corr to {args.output_all}")
    with open(args.output_all, "wb") as f_out:
        for tmp_fastq in tmp_fastqs:
            with open(tmp_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info(f"Writing filtered reads with bc_corr to {args.output_filtered}")
    with open(args.output_filtered, "wb") as f_out:
        for tmp_filt_fastq in tmp_filt_fastqs:
            with open(tmp_filt_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info("Cleaning up")
    [os.remove(fn) for fn in tmp_fastqs]
    [os.remove(fn) for fn in tmp_filt_fastqs]
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
