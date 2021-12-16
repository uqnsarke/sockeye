import argparse
import gzip
import logging
import math
import multiprocessing
import os
import pathlib
import shutil
import tempfile

import editdistance as ed
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "fasta", help="FASTA file of putative barcode sequences", type=str
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
        "-o",
        "--output",
        help="Output file [bc_assigns.tsv]",
        type=str,
        default="bc_assigns.tsv",
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
    p = pathlib.Path(args.output)
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
    def __init__(self, read_name, found_bc):
        """
        Define basic characteristic of the barcode that we
        want to classify
        """
        self.found_bc = found_bc
        self.read_id = read_name.split(" ")[0]
        self.rlen = int(read_name.split(" ")[1].split("=")[1])
        self.read1_ed = int(read_name.split(" ")[2].split("=")[1])
        self.qscore = float(read_name.split(" ")[3].split("=")[1])
        self.umi = read_name.split(" ")[4].split("=")[1]

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


def calc_ed(tup):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment
    """
    read_name = tup[0]
    seq = tup[1]
    kmer_to_bc_index = tup[2]
    whitelist = tup[3]
    args = tup[4]

    kmers = split_barcode_into_kmers(seq, args.k)
    filt_whitelist = filter_whitelist_by_kmers(whitelist, kmers, kmer_to_bc_index)
    # filt_whitelist = whitelist

    # Instantiate barcode object
    ReadBC = Barcode(read_name, seq)
    # Calc edit distances and store min & runner up
    ReadBC.calc_ed_with_whitelist(filt_whitelist)
    return ReadBC


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


def launch_alignment_pool(batch_barcodes, whitelist, kmer_to_bc_index, args):
    func_args = []

    for read_num, (read_name, barcode) in enumerate(batch_barcodes):
        # filt_whitelist = whitelist
        func_args.append((read_name, barcode, kmer_to_bc_index, whitelist, args))

    fasta_entries = launch_pool(args.threads, calc_ed, func_args)
    return fasta_entries


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
    # for kmer,indexes in kmer_to_bc_index.items():
    #     print(kmer, indexes)
    return wl, kmer_to_bc_index


def open_fasta(fasta):
    if fasta.split(".")[-1] == "gz":
        f = gzip.open(fasta, "rt")
    else:
        f = open(fasta, "r+")
    return f


def check_input_format(fasta):
    f = open_fasta(fasta)
    line = f.readline()
    if line[0] == ">":
        pass
    else:
        raise ("Unexpected file type! FASTA only!")
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


def count_barcodes(fasta):
    print("Counting putative barcodes")
    number_lines = 0
    with open_fasta(fasta) as f:
        for line in tqdm(f, unit_scale=0.5, unit=" barcodes"):
            number_lines += 1
    return int(number_lines / 2)


def write_tmp_table(barcodes):
    tmp_table = tempfile.NamedTemporaryFile(
        prefix="tmp.assigns.", suffix=".tsv", dir=args.tempdir, delete=False
    )
    d = [vars(ReadBC) for ReadBC in barcodes]

    df = pd.DataFrame.from_dict(d).set_index("read_id")
    df.to_csv(tmp_table.name, sep="\t", index=True)
    return tmp_table.name


def main(args):
    init_logger(args)
    whitelist, kmer_to_bc_index = load_whitelist(args)
    check_input_format(args.fasta)
    n_barcodes = count_barcodes(args.fasta)
    n_batches = int(math.ceil(n_barcodes / args.batch_size))

    f = open_fasta(args.fasta)
    record_iter = SimpleFastaParser(f)

    logger.info(
        "Processing {n} barcodes in {b} batches".format(n=n_barcodes, b=n_batches)
    )
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)
    tmp_tables = []
    for i, batch_barcodes in enumerate(
        tqdm(batch_iterator(record_iter, args.batch_size), total=n_batches)
    ):

        barcodes = launch_alignment_pool(
            batch_barcodes, whitelist, kmer_to_bc_index, args
        )
        tmp_table = write_tmp_table(barcodes)
        tmp_tables.append(tmp_table)

    # Merge temp demux tables then clean up
    demux_table_fn = args.output
    logger.info(f"Writing read demux table to {demux_table_fn}")
    df = pd.concat([pd.read_csv(fn, sep="\t") for fn in tmp_tables], axis=0)
    df = df.drop(["runner_up_bc", "runner_up_ed"], axis=1)
    df.to_csv(demux_table_fn, sep="\t", index=False)

    logger.info("Cleaning up")
    [os.remove(fn) for fn in tmp_tables]
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
