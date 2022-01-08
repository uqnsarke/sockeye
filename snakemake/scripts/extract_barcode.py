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

import editdistance as ed
import numpy as np
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

    parser.add_argument("--umi_length", help="UMI length [10]", type=int, default=10)

    parser.add_argument(
        "--bc_superlist",
        help="If specified, also outputs barcodes that are in \
                        the supplied barcode superlist. This option will \
                        write bc_umi sequences to the output file \
                        specified by --output_hq_bc [None]",
        type=str,
        default=None,
    )

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
        "--max_read1_ed",
        help="Max edit distance with read1 adapter sequence (upstream of cell barcode) [3]",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--output_reads",
        help="Output FASTQ file containing stranded reads with uncorrected barcodes \
        in the read header [reads.bc_uncorr.fastq]",
        type=str,
        default="reads.bc_uncorr.fastq",
    )

    parser.add_argument(
        "--output_hq_bc",
        help="Output FASTA file containing high-quality barcode entries [bc_uncorr.hq.fasta]",
        type=str,
        default="bc_uncorr.hq.fasta",
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
    p = pathlib.Path(args.output_reads)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    if args.bc_superlist == "None":
        args.bc_superlist = None

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


def edit_distance(query, target):
    return ed.eval(query, target)


def align_adapter(tup):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment
    """
    fastq_entry = tup[0]
    args = tup[1]

    read_id = fastq_entry.id
    read_seq = str(fastq_entry.seq)
    read_qv = fastq_entry.letter_annotations["phred_quality"]

    def create_fasta_entry(barcode_only, read_id, read1_ed, bc_qv):
        fasta_entry = SeqRecord(
            Seq(barcode_only),
            id=read_id,
            description="read1_ed={} bc_qv={:.1f}".format(read1_ed, np.mean(bc_qv)),
        )
        return fasta_entry

    matrix = update_matrix(args)

    # Use only the specified suffix length of the read1 adapter
    read1_probe_seq = args.read1_adapter[-args.read1_suff_length :]
    # Compile the actual query sequence of <read1_suffix>NNN...NNN<TTTTT....>
    probe_seq = "{r}{bc}{umi}{pT}".format(
        r=read1_probe_seq,
        bc="N" * args.barcode_length,
        umi="N" * args.umi_length,
        pT="T" * args.polyT_length,
    )

    prefix_seq = read_seq[: args.window]
    prefix_qv = read_qv[: args.window]

    alignment = parasail_alg(
        s1=prefix_seq,
        s2=probe_seq,
        open=args.gap_open,
        extend=args.gap_extend,
        matrix=matrix,
    )

    # Find the position of the Ns in the alignment. These correspond
    # to the cell barcode + UMI sequences bound by the read1 and polyT
    idxs = list(find("N", alignment.traceback.ref))
    if len(idxs) > 0:
        # The Ns in the probe successfully aligned to sequence
        bc_start = min(idxs)
        umi_end = max(idxs)

        # The read1 adapter comprises the first part of the alignment
        read1 = alignment.traceback.query[0:bc_start]
        read1_ed = edit_distance(read1, read1_probe_seq)

        # The barcode + UMI sequences in the read correspond to the
        # positions of the aligned Ns in the probe sequence
        barcode_umi = alignment.traceback.query[bc_start : umi_end + 1]
        barcode_only = alignment.traceback.query[
            bc_start : (bc_start + args.barcode_length)
        ]
    else:
        # No Ns in the probe successfully aligned -- ignore this read
        read1_ed = len(read1_probe_seq)
        barcode_umi = ""
        barcode_only = ""

    barcode_umi = barcode_umi.strip("-")

    # Require minimal read1 edit distance and require perfect barcode length
    condition1 = read1_ed <= args.max_read1_ed
    ##########################
    # TO-DO: can we be more flexible about barcode length?
    ##########################
    condition2 = len(barcode_only.strip("-")) == args.barcode_length

    if condition1 and condition2:
        start_idx = prefix_seq.find(barcode_umi)
        bc_qv = prefix_qv[start_idx : (start_idx + args.barcode_length + 1)]

        # print(alignment.traceback.ref)
        # print(alignment.traceback.comp)
        # print(alignment.traceback.query)
        # print()

        fasta_entry = create_fasta_entry(barcode_only, read_id, read1_ed, bc_qv)
        fastq_entry.description = "{} bc_uncorr={} bc_qv={:.2f}".format(
            fastq_entry.id, barcode_only, np.mean(bc_qv)
        )
    else:
        fasta_entry = None
        fastq_entry = None

    return fasta_entry, fastq_entry


def launch_alignment_pool(batch_entries, args):
    func_args = []

    for r in batch_entries:
        if len(str(r.seq)) > 0:
            func_args.append((r, args))

    results = launch_pool(args.threads, align_adapter, func_args)
    fasta_entries, fastq_entries = list(zip(*results))
    return fasta_entries, fastq_entries


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


def load_superlist(args):
    wl = []
    with open(args.bc_superlist) as file:
        for line in tqdm(file, total=get_num_lines(args.bc_superlist)):
            wl.append(line.strip())
    return set(np.array(wl))


def open_fastq(fastq):
    if fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq, "r+")
    return f


def mmap_line_count(f):
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    while buf.readline():
        lines += 1
    return lines


def get_num_lines(file_path):
    f = open(file_path, "r+")
    lines = mmap_line_count(f)
    return lines


def count_reads(fastq):
    logger.info("Counting reads")
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in tqdm(f, unit_scale=0.25, unit=" reads"):
            number_lines += 1
    return number_lines / 4


def write_tmp_files(fasta_entries, fastq_entries, wl, args):
    hq_fastas = []
    for fasta_entry in fasta_entries:
        if fasta_entry is not None:
            if args.bc_superlist is not None:
                if fasta_entry.seq in wl:
                    hq_fastas.append(fasta_entry)

    hq_tmp_fasta = tempfile.NamedTemporaryFile(
        prefix="tmp.hq.bc.", suffix=".fasta", dir=args.tempdir, delete=False
    )
    tmp_read_fastq = tempfile.NamedTemporaryFile(
        prefix="tmp.read.bc_uncorr.", suffix=".fastq", dir=args.tempdir, delete=False
    )
    SeqIO.write(hq_fastas, hq_tmp_fasta.name, "fasta")
    valid_fastq_entries = [entry for entry in fastq_entries if entry is not None]
    SeqIO.write(valid_fastq_entries, tmp_read_fastq.name, "fastq")
    return hq_tmp_fasta.name, tmp_read_fastq.name


def main(args):
    init_logger(args)
    check_input_format(args.fastq)
    n_reads = int(count_reads(args.fastq))
    n_batches = int(n_reads / args.batch_size)

    if args.bc_superlist is not None:
        logger.info("Loading barcode superlist")
        wl = load_superlist(args)
    else:
        wl = None

    f = open_fastq(args.fastq)
    record_iter = SeqIO.parse(f, "fastq")

    logger.info("Processing {n} reads in {b} batches".format(n=n_reads, b=n_batches))
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)

    tmp_read_fastqs = []
    hq_tmp_fastas = []
    for i, batch_entries in enumerate(
        tqdm(batch_iterator(record_iter, args.batch_size), total=n_batches)
    ):

        fasta_entries, fastq_entries = launch_alignment_pool(batch_entries, args)
        hq_tmp_fasta, tmp_read_fastq = write_tmp_files(
            fasta_entries, fastq_entries, wl, args
        )
        tmp_read_fastqs.append(tmp_read_fastq)
        hq_tmp_fastas.append(hq_tmp_fasta)

    logger.info(
        f"Writing reads with putative barcodes in header to {args.output_reads}"
    )
    with open(args.output_reads, "wb") as f_out:
        for tmp_fastq in tmp_read_fastqs:
            with open(tmp_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info(f"Writing superlist-filtered putative barcodes to {args.output_hq_bc}")
    with open(args.output_hq_bc, "wb") as f_out:
        for tmp_fasta in hq_tmp_fastas:
            with open(tmp_fasta, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info("Cleaning up")
    [os.remove(fn) for fn in tmp_read_fastqs]
    [os.remove(fn) for fn in hq_tmp_fastas]
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
