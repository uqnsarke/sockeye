import argparse
import gzip
import logging
import mmap
import multiprocessing
import os
import pathlib
import shutil
import sys
import tempfile

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

try:
    import parasail
except ImportError:
    logging.error(
        "Could not load parasail library. Please try reinstalling "
        "parasail using pip."
    )
    sys.exit(1)

parasail_alg = parasail.sw_trace

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "fastq", help="Stranded sequencing reads (fastq or fasta)", type=str
    )

    parser.add_argument(
        "barcodes", help="Table of filtered barcode assignments (tsv)", type=str
    )

    # Optional arguments
    parser.add_argument(
        "-r",
        "--read1_adapter",
        help="Read1 adapter sequence to use in the alignment query [CTACACGACGCTCTTCCGATCT]",
        default="CTACACGACGCTCTTCCGATCT",
        type=str,
    )

    parser.add_argument(
        "--read1_suff_length",
        help="Use this many suffix bases from read1 sequence \
                        in the alignment query. For example, specifying 12 \
                        would mean that the last 12 bases of the specified \
                        read1 sequence will be included in the probe sequence \
                        [10]",
        default=10,
        type=int,
    )

    parser.add_argument(
        "-T",
        "--polyT_length",
        help="Length of polyT sequence to use in the alignment query [10]",
        type=int,
        default=10,
    )

    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "--barcode_length", help="Cell barcode length [16]", type=int, default=16
    )

    parser.add_argument("--umi_length", help="UMI length [12]", type=int, default=12)

    parser.add_argument(
        "-o", "--gap_open", help="Gap open penalty [2]", type=int, default=2
    )

    parser.add_argument(
        "-e", "--gap_extend", help="Gap extend penalty [4]", type=int, default=4
    )

    parser.add_argument("-m", "--match", help="Match score [5]", type=int, default=5)

    parser.add_argument(
        "-x", "--mismatch", help="Mismatch score [-1]", type=int, default=-1
    )

    parser.add_argument(
        "-n",
        "--acg_to_n_match",
        help="Score for A/C/G<-->N match [1]",
        type=int,
        default=1,
    )

    parser.add_argument(
        "-s", "--t_to_n_match", help="Score for T<-->N match [1]", type=int, default=1
    )

    parser.add_argument(
        "-w",
        "--window",
        help="Number of bases to query at start of read [100]",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--output",
        help="Output file containing stranded read, uncorrected cell barcode, \
                    corrected cell barcode, cell barcode QV, uncorrected UMI, \
                    and UMI QV entries [bc_extracted_umi.fastq]",
        type=str,
        default="bc_extracted_umi.fastq",
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


def update_matrix(args):
    """
    Create new parasail scoring matrix. 'N' is used as wildcard character
    for barcodes and has its own match parameter (0 per default).
    :return: args: arguments
    :rtype args, args
    """
    matrix = parasail.matrix_create("ACGTN", args.match, args.mismatch)

    ############################
    # SCORING MATRIX POSITIONS #
    ############################
    #     A   C   G   T   N
    # A   0   1   2   3   4   5
    # C   6   7   8   9   10  11
    # G   12  13  14  15  16  17
    # T   18  19  20  21  22  23
    # N   24  25  26  27  28  29

    # Update scoring matrix so that N matches A/C/G/N
    # pointers = [4, 9, 14, 20, 21, 22, 24]
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        matrix.pointer[0].matrix[i] = args.acg_to_n_match

    # Update N to T matching score so that we start
    # aligning dT sequence instead of Ns.
    matrix.pointer[0].matrix[22] = args.t_to_n_match
    matrix.pointer[0].matrix[27] = args.t_to_n_match
    return matrix


def launch_pool(procs, funct, args):
    p = multiprocessing.Pool(processes=procs)
    try:
        results = p.map(funct, args)
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


def find(target, myList):
    for i in range(len(myList)):
        if myList[i] == target:
            yield i


def align_query(tup):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment
    """
    read_id = tup[0]
    read_seq = tup[1]
    read_qv = tup[2]
    bc_corr = tup[3]
    bc_uncorr = tup[4]
    bc_qv = tup[5]
    args = tup[6]

    prefix_seq = read_seq[: args.window]
    prefix_qv = read_qv[: args.window]

    matrix = update_matrix(args)

    # Use only the specified suffix length of the read1 adapter
    read1_probe_seq = args.read1_adapter[-args.read1_suff_length :]
    # Compile the actual query sequence of <read1_suffix><bc_corr>NNN...N<TTTTT....>
    probe_seq = "{r}{bc}{umi}{pT}".format(
        r=read1_probe_seq,
        bc=bc_corr,
        umi="N" * args.umi_length,
        pT="T" * args.polyT_length,
    )

    alignment = parasail_alg(
        s1=prefix_seq,
        s2=probe_seq,
        open=args.gap_open,
        extend=args.gap_extend,
        matrix=matrix,
    )

    idxs = list(find("N", alignment.traceback.ref))
    if len(idxs) > 0:
        umi_start = min(idxs)
        umi_end = max(idxs)
        umi = alignment.traceback.query[umi_start : umi_end + 1]

        umi = umi.strip("-")
        start_idx = prefix_seq.find(umi)
        umi_qv = prefix_qv[start_idx : (start_idx + len(umi))]
        # print("{} bc_corr={} bc_uncorr={} bc_qv={:.1f} umi_uncorr={} umi_qv={:.1f}".format(read_id, bc_corr, bc_uncorr, bc_qv, umi, np.mean(umi_qv)))
        # print(alignment.traceback.ref)
        # print(alignment.traceback.comp)
        # print(alignment.traceback.query)
        # print()
    else:
        # No Ns in the alignment -- ignore
        umi = "XXXXXXXXX"
        umi_qv = 0.0

    fastq_entry = SeqRecord(
        Seq(read_seq),
        id=read_id,
        description="bc_corr={} bc_uncorr={} bc_qv={:.1f} umi_uncorr={} umi_qv={:.1f}".format(
            bc_corr, bc_uncorr, bc_qv, umi, np.mean(umi_qv)
        ),
    )
    fastq_entry.letter_annotations["phred_quality"] = read_qv

    return fastq_entry


def launch_alignment_pool(batch_reads, args):
    func_args = []

    for r in batch_reads:
        func_args.append(
            (
                r.read_id,
                r.seq,
                list(map(lambda x: ord(x) - 33, r.qv)),
                r.bc_corr,
                r.bc_uncorr,
                r.bc_qv,
                args,
            )
        )

    read_entries = launch_pool(args.threads, align_query, func_args)
    return read_entries


def check_input_format(fastq):
    f = open_fastq(fastq)
    line = f.readline()

    if line[0] == "@":
        pass
    else:
        raise ("Unexpected file type! Only *.fastq, and *.fq recognized.")
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


def open_fastq(fastq):
    if fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq, "r+")
    return f


def count_reads(fastq):
    logger.info("Counting reads")
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in tqdm(f, unit_scale=0.25, unit=" reads"):
            number_lines += 1
    return number_lines / 4


def write_tmp_files(fastq_entries, args):
    tmp_fastq_file = tempfile.NamedTemporaryFile(
        prefix="tmp.bc.", suffix=".fastq", dir=args.tempdir, delete=False
    )
    SeqIO.write(fastq_entries, tmp_fastq_file.name, "fastq")
    return tmp_fastq_file.name


class Read:
    def __init__(self, read_id, seq, qv):
        self.read_id = read_id
        self.seq = seq
        self.qv = qv

    def add_barcode_info(self, bc_uncorr, bc_corr, bc_qv):
        self.bc_uncorr = bc_uncorr
        self.bc_corr = bc_corr
        self.bc_qv = bc_qv


def get_read_batch_barcodes(batch_entries, args):
    # Load all assigned barcode information
    df = pd.read_csv(args.barcodes, sep="\t").set_index("read_id")

    batch_reads = []
    # Create Read objects for just those reads in the batch
    for (read_name, seq, qv) in batch_entries:
        read_id = read_name.split(" ")[0]
        r = Read(read_id, seq, qv)
        if read_id in df.index:
            r.add_barcode_info(
                df.loc[read_id, "bc_uncorr"],
                df.loc[read_id, "bc_corr"],
                df.loc[read_id, "bc_qv"],
            )
            batch_reads.append(r)
    return batch_reads


def main(args):
    init_logger(args)
    check_input_format(args.fastq)
    n_reads = int(count_reads(args.fastq))
    n_batches = int(n_reads / args.batch_size)

    f = open_fastq(args.fastq)
    record_iter = FastqGeneralIterator(f)

    logger.info("Processing {n} reads in {b} batches".format(n=n_reads, b=n_batches))
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)

    tmp_files = []
    for i, batch_entries in enumerate(
        tqdm(batch_iterator(record_iter, args.batch_size), total=n_batches)
    ):

        batch_reads = get_read_batch_barcodes(batch_entries, args)
        read_entries = launch_alignment_pool(batch_reads, args)
        tmp_batch_files = write_tmp_files(read_entries, args)
        tmp_files.append(tmp_batch_files)

    logger.info(f"Writing all barcoded reads to {args.output}")
    with open(args.output, "wb") as f_out:
        for tmp_read in tmp_files:
            with open(tmp_read, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info("Cleaning up")
    [os.remove(fn) for fn in tmp_files]
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
