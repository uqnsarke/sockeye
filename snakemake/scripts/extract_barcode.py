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
import parasail
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

# Make logger globally accessible to all functions
logger = logging.getLogger(__name__)


def parse_args():
    """
    Parse the command line arguments

    :return args: object containing all supplied arguments
    :rtype args: class argparse.Namespace
    """
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "fastq",
        help="FASTQ file of stranded sequencing reads (gzipped ending in *.fastq.gz \
        or *.fq.gz supported)",
        type=str,
    )

    parser.add_argument(
        "superlist",
        help="Comprehensive whitelist of all possible cell barcodes. For example, \
        the file 3M-february-2018.txt.gz can be downloaded at https://github.com/\
        10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation\
        /3M-february-2018.txt.gz",
        type=str,
        default=None,
    )

    # Optional arguments
    parser.add_argument(
        "-r",
        "--read1_adapter",
        help="Read1 adapter sequence to use in the alignment query \
            [CTACACGACGCTCTTCCGATCT]",
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
        "-t",
        "--threads",
        help="Threads to use [4]",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--barcode_length",
        help="Cell barcode length [16]",
        type=int,
        default=16,
    )

    parser.add_argument(
        "--umi_length",
        help="UMI length [10]",
        type=int,
        default=10,
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
        help="Max edit distance with read1 adapter sequence (upstream of cell \
        barcode) [3]",
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
        help="Output FASTA file containing high-quality barcode entries \
        [bc_uncorr.hq.fasta]",
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
    output_dir = p.parents[0]
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=output_dir)
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    """
    Initialize the logger using the specified verbosity level.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
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

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return: matrix: custom parasail alignment matrix
    :rtype matrix: parasail.bindings_v2.Matrix
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
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        matrix.pointer[0].matrix[i] = args.acg_to_n_match

    # Update N to T matching score so that we start
    # aligning dT sequence instead of Ns.
    matrix.pointer[0].matrix[22] = args.t_to_n_match
    matrix.pointer[0].matrix[27] = args.t_to_n_match
    return matrix


def launch_pool(func, func_args, procs=1):
    """
    Use multiprocessing library to create pool and map function calls to
    that pool

    :param procs: Number of processes to use for pool
    :type procs: int, optional
    :param func: Function to exececute in the pool
    :type func: function
    :param func_args: List containing arguments for each call to function <funct>
    :type func_args: list
    :return results: List of results returned by each call to function <funct>
    :rtype results: list
    """
    p = multiprocessing.Pool(processes=procs)
    try:
        results = p.map(func, func_args)
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


def find(char, string):
    """
    Return iterator of indices for positions in a string
    corresponding to a target character

    :param char: Target character whose positions you want to locate in the string
    :type char: str
    :param string: String to search for the target character positions
    :type string: str
    :return i: Indices in the string corresponding to the target character
    :rtype i: iterator
    """
    for i in range(len(string)):
        if string[i] == char:
            yield i


def edit_distance(query, target):
    """
    Return Levenshtein distance between the two supplied strings

    :param query: Query string to compare against the target
    :type query: str
    :param target: Target string to compare against the query
    :type target: str
    :return d: Calculated Levenshtein distance between query and target
    :rtype d: int
    """
    d = ed.eval(query, target)
    return d


def align_adapter(tup):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment

    :param tup: Tuple containing the function arguments
    :type tup: tup
    :return fasta_entry: Biopython FASTA Seq record containing identified (uncorrected)
        cell barcode sequence
    :rtype fasta_entry: class 'Bio.SeqRecord.SeqRecord'
    :return fastq_entry: Biopython FASTQ Seq record containing the read with the
        identified cell barcode (bc_uncorr) and barcode mean quality value (bc_qv)
        added to the read header
    :rtype fastq_entry: class 'Bio.SeqRecord.SeqRecord'
    """
    fastq_entry = tup[0]
    args = tup[1]

    read_id = fastq_entry.id
    read_seq = str(fastq_entry.seq)
    read_qv = fastq_entry.letter_annotations["phred_quality"]

    def create_fasta_entry(barcode, read_id, read1_ed, bc_qv):
        """
        Create a Biopython FASTA SeqRecord using the supplied sequence info

        :param barcode: Nucleotide sequence of the identified (uncorrected) cell barcode
        :type barcode: str
        :param read_id: Read ID
        :type read_id: str
        :param read1_ed: Levenshtein distance between expected and observed read1
            adapter sequence
        :type read1_ed: int
        :param bc_qv: Mean quality value for positions corresponding to the identified
            cell barcode
        :type bc_qv: float
        :return fasta_entry: Biopython FASTA Seq record containing the identified
            (uncorrected) cell barcode sequence
        :rtype fasta_entry: class 'Bio.SeqRecord.SeqRecord'
        """
        fasta_entry = SeqRecord(
            Seq(barcode),
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

    parasail_alg = parasail.sw_trace

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

        print(read_id)
        print(alignment.traceback.ref)
        print(alignment.traceback.comp)
        print(alignment.traceback.query)
        print()

        fasta_entry = create_fasta_entry(barcode_only, read_id, read1_ed, bc_qv)
        fastq_entry.description = "{} bc_uncorr={} bc_qv={:.2f}".format(
            fastq_entry.id, barcode_only, np.mean(bc_qv)
        )
    else:
        fasta_entry = None
        fastq_entry = None

    return fasta_entry, fastq_entry


def launch_alignment_pool(batch, args):
    """
    First builds a list of tuples (<func_args>), where each tuple contains the
    arguments required by the function that is being executed by the
    multiprocessing pool. Then launch the pool and return the results.

    :param batch: List of FASTQ SeqRecord entries batched from the input FASTQ file
    :type batch: list
    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return fasta_entries: List of Biopython FASTA Seq records containing
        identified (uncorrected) cell barcode sequences
    :rtype fasta_entries: list
    :return fastq_entries: List of Biopython FASTQ Seq records containing the
        reads with the identified cell barcode (bc_uncorr) and barcode mean
        quality value (bc_qv) added to the read headers
    :rtype fastq_entries: list
    """
    func_args = []

    for r in batch:
        if len(str(r.seq)) > 0:
            func_args.append((r, args))

    results = launch_pool(align_adapter, func_args, args.threads)
    fasta_entries, fastq_entries = list(zip(*results))
    return fasta_entries, fastq_entries


def check_input_format(input_file):
    """
    Check that the input is a FASTQ file by just reading the first line.

    :param input_file: Filename to check
    :type input_file: str
    """
    f = open_fastq(input_file)
    line = f.readline()

    if line[0] == "@":
        pass
    else:
        raise ("Unexpected file type! Check your FASTQ file.")


def batch_iterator(iterator, batch_size):
    """
    Returns lists of length batch_size.

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


def load_superlist(superlist):
    """
    Read contents of the file containing all possible cell barcode sequences.
    File can be uncompressed or gzipped.

    :param superlist: Path to file containing all possible cell barcodes, e.g.
        3M-february-2018.txt
    :type superlist: str
    :return wl: Set of all possible cell barcodes
    :rtype wl: set
    """
    ext = pathlib.Path(superlist).suffix
    fn = pathlib.Path(superlist).name
    wl = []
    if ext == ".gz":
        with gzip.open(superlist, "rt") as file:
            for line in tqdm(file, desc=f"Loading barcodes in {fn}", unit=" barcodes"):
                wl.append(line.strip())
    elif ext == ".txt":
        with open(superlist) as file:
            for line in tqdm(file, desc=f"Loading barcodes in {fn}", unit=" barcodes"):
                wl.append(line.strip())
    wl = set(wl)
    return wl


def open_fastq(fastq):
    """
    Open the supplied FASTQ file, depending on whether it is gzipped or uncompressed

    :param fastq: Path to supplied FASTQ file
    :type fastq: str
    :return f: File handle for FASTQ file
    :rtype f: class '_io.TextIOWrapper'
    """
    ext = pathlib.Path(fastq).suffix
    if ext == ".gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq, "r+")
    return f


def count_reads(fastq):
    """
    Get the number of reads from a FASTQ file. Does this by simply dividing the
    number of lines in the file by four.

    :param fastq: Path to supplied FASTQ file
    :type fastq: str
    :return n_reads: Number of reads in FASTQ file
    :rtype n_reads: int
    """
    n_lines = 0
    fn = pathlib.Path(fastq).name
    with open_fastq(fastq) as f:
        for line in tqdm(
            f, desc=f"Counting reads in {fn}", unit_scale=0.25, unit=" reads"
        ):
            n_lines += 1
    n_reads = int(n_lines / 4)
    return n_reads


def write_tmp_files(fasta_entries, fastq_entries, wl, args):
    """
    Write temporary files for each batch of reads. Each batch will generate
    (1) a temporary FASTQ of reads with the barcode info in the read headers,
    (2) a temporary FASTA of the uncorrected cell barcode found in each read.

    :param fasta_entries: List of Biopython FASTA Seq records containing
        identified (uncorrected) cell barcode sequences
    :type fasta_entries: list
    :param fastq_entries: List of Biopython FASTQ Seq records containing the
        reads with the identified cell barcode (bc_uncorr) and barcode mean
        quality value (bc_qv) added to the read headers
    :type fastq_entries: list
    :param wl:
    :type wl: set
    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return hq: Path to temporary FASTA file of cell barcodes found in each read
    :rtype hq: str
    :return read: Path to temporary FASTQ file of reads with cell barcode info
        in the read headers
    :rtype read: str
    """
    hq_fastas = []
    for fasta_entry in fasta_entries:
        if fasta_entry is not None:
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
    hq = hq_tmp_fasta.name
    read = tmp_read_fastq.name
    return hq, read


def main(args):
    init_logger(args)
    check_input_format(args.fastq)
    n_reads = count_reads(args.fastq)
    n_batches = int(np.ceil(n_reads / args.batch_size))

    # logger.info("Loading barcode superlist")
    wl = load_superlist(args.superlist)

    f = open_fastq(args.fastq)
    record_iter = SeqIO.parse(f, "fastq")

    logger.info(
        f"Processing {n_reads} total reads in {n_batches} batches of "
        f"{args.batch_size} reads"
    )

    # Create temporary directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)

    # Process reads from the FASTQ file in batches to manage memory usage
    tmp_read_fastqs = []
    tmp_hq_fastas = []
    batches = batch_iterator(record_iter, args.batch_size)
    for batch in tqdm(batches, total=n_batches):

        fasta_entries, fastq_entries = launch_alignment_pool(batch, args)
        hq_tmp_fasta, tmp_read_fastq = write_tmp_files(
            fasta_entries, fastq_entries, wl, args
        )
        tmp_read_fastqs.append(tmp_read_fastq)
        tmp_hq_fastas.append(hq_tmp_fasta)

    logger.info(
        f"Writing reads with putative barcodes in header to {args.output_reads}"
    )
    with open(args.output_reads, "wb") as f_out:
        for tmp_fastq in tmp_read_fastqs:
            with open(tmp_fastq, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info(f"Writing superlist-filtered putative barcodes to {args.output_hq_bc}")
    with open(args.output_hq_bc, "wb") as f_out:
        for tmp_fasta in tmp_hq_fastas:
            with open(tmp_fasta, "rb") as f_:
                shutil.copyfileobj(f_, f_out)

    logger.info("Cleaning up temporary files")
    [os.remove(fn) for fn in tmp_read_fastqs]
    [os.remove(fn) for fn in tmp_hq_fastas]
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
