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
import numpy as np
import parasail
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Sorted BAM file of stranded sequencing reads aligned to a reference. \
            Alignments must have the CR and CY tags.",
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
        "--output",
        help="Output BAM file containing aligned reads with tags for uncorrected \
        barcodes (CR), corrected barcodes (CB), barcode QVs (CY), uncorrected \
        UMIs (UR), and UMI QVs (UY) [bc_corr.umi_uncorr.sorted.bam]",
        type=str,
        default="bc_corr.umi_uncorr.sorted.bam",
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
        "--barcode_length", help="Cell barcode length [16]", type=int, default=16
    )

    parser.add_argument("--umi_length", help="UMI length [12]", type=int, default=12)

    parser.add_argument(
        "-w",
        "--window",
        help="Number of bases to query at start of read [100]",
        type=int,
        default=100,
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


def find(target, myList):
    for i in range(len(myList)):
        if myList[i] == target:
            yield i


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


def calc_ed_with_whitelist(bc_uncorr, whitelist):
    """
    Find minimum and runner-up barcode edit distance by
    iterating through the whitelist of expected barcodes
    """
    bc_match = "X" * len(bc_uncorr)
    bc_match_ed = len(bc_uncorr)
    next_bc_match_ed = len(bc_uncorr)
    for wl_bc in whitelist:
        d = ed.eval(bc_uncorr, wl_bc)
        if d < bc_match_ed:
            next_bc_match_ed = bc_match_ed
            bc_match_ed = d
            bc_match = wl_bc
    next_match_diff = next_bc_match_ed - bc_match_ed

    return bc_match, bc_match_ed, next_match_diff


def parse_probe_alignment(p_alignment, align, prefix_seq, prefix_qv):
    """
    Parse a parasail alignment alignment and add uncorrected UMI and UMI QV
    values as tags to the BAM alignment.

    :param p_alignment: parasail alignment
    :type p_alignment:
    :param align: pysam BAM alignment
    :type align:
    :return: pysam BAM alignment with the UR and UY tags added
    :rtype:
    """
    # Find the position of the Ns in the parasail alignment. These correspond
    # to the UMI sequences bound by the cell barcode and polyT
    idxs = list(find("N", p_alignment.traceback.ref))
    if len(idxs) > 0:
        umi = p_alignment.traceback.query[min(idxs) : max(idxs) + 1]

        umi = umi.strip("-")
        start_idx = prefix_seq.find(umi)
        umi_qv = np.mean(prefix_qv[start_idx : (start_idx + len(umi))])

        # print(alignment.traceback.ref)
        # print(alignment.traceback.comp)
        # print(alignment.traceback.query)
        # print()

        # Uncorrected UMI = UR:Z
        align.set_tag("UR", umi, value_type="Z")
        # UMI quality score = UY:Z
        align.set_tag("UY", "{:.2f}".format(umi_qv), value_type="Z")

    return align


def get_uncorrected_umi(align, args):
    """
    Aligns a probe sequence containing the read1+corrected_barcode+Ns+polyT to
    the read. Bases aligning to the Ns in the probe sequence correspond to the
    UMI positions. Extract those bases and consider those to be the uncorrected
    UMI sequence.

    :param align:
    :type align:
    :param args:
    :type args:
    :return:
    :rtype:
    """
    prefix_seq = align.get_forward_sequence()[: args.window]
    prefix_qv = align.get_forward_qualities()[: args.window]

    # Use only the specified suffix length of the read1 adapter
    read1_probe_seq = args.read1_adapter[-args.read1_suff_length :]
    # Compile the actual query sequence of <read1_suffix><bc_corr>NNN...N<TTTTT....>
    probe_seq = "{r}{bc}{umi}{pT}".format(
        r=read1_probe_seq,
        bc=align.get_tag("CB"),
        umi="N" * args.umi_length,
        pT="T" * args.polyT_length,
    )

    matrix = update_matrix(args)
    parasail_alg = parasail.sw_trace

    p_alignment = parasail_alg(
        s1=prefix_seq,
        s2=probe_seq,
        open=args.gap_open,
        extend=args.gap_extend,
        matrix=matrix,
    )

    align = parse_probe_alignment(p_alignment, align, prefix_seq, prefix_qv)

    return align


def process_bam_records(tup):
    """
    Process BAM records to assign each read a corrected cell barcode and an
    uncorrected UMI. Do this by loading and processing the barcode whitelist
    then iterating over alignments.
    For each alignment:
    1. Calculate edit distance between uncorrected barcode and barcodes in whitelist
    2.

    :param tup:
    :type tup: tup
    :return:
    :rtype:
    """
    input_bam = tup[0]
    chrom = tup[1]
    args = tup[2]

    # Load barcode whitelist and map kmers to indices in whitelist for faster
    # barcode matching
    whitelist, kmer_to_bc_index = load_whitelist(args.whitelist, args.k)

    # Open input BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    # Open temporary output BAM file for writing
    suff = f".{chrom}.bam"
    chrom_bam = tempfile.NamedTemporaryFile(
        prefix="tmp.align.", suffix=suff, dir=args.tempdir, delete=False
    )
    bam_out = pysam.AlignmentFile(chrom_bam.name, "wb", template=bam)

    for align in bam.fetch(contig=chrom):
        # Make sure each alignment in this BAM has an uncorrected baracode and
        # barcode QV
        assert align.has_tag("CR") and align.has_tag("CY"), "CR or CY tags not found"

        bc_uncorr = align.get_tag("CR")

        # Decompose uncorrected barcode into N k-mers
        bc_uncorr_kmers = split_barcode_into_kmers(bc_uncorr, args.k)
        # Filter the whitelist to only those with at least one of the k-mers
        # from the uncorrected barcode
        filt_whitelist = filter_whitelist_by_kmers(
            whitelist, bc_uncorr_kmers, kmer_to_bc_index
        )

        # Calc edit distances between uncorrected barcode and the filtered
        # whitelist barcodes
        bc_match, bc_match_ed, next_match_diff = calc_ed_with_whitelist(
            bc_uncorr, filt_whitelist
        )

        # Check barcode match edit distance and difference to runner-up edit distance
        condition1 = bc_match_ed <= args.max_ed
        condition2 = next_match_diff >= args.min_ed_diff
        if condition1 and condition2:
            # Add corrected cell barcode = CB:Z
            align.set_tag("CB", bc_match, value_type="Z")

            # Add corrected barcode to probe sequence to fish out uncorrected UMI
            align = get_uncorrected_umi(align, args)

            # Only write BAM entry in output file if we've assigned a corrected
            # barcode and an uncorrected UMI
            bam_out.write(align)

    bam.close()
    bam_out.close()

    return chrom_bam.name


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
    :return: List of results returned by each call to function <funct>
    :rtype: list
    """
    p = multiprocessing.Pool(processes=procs)
    try:
        results = list(tqdm(p.imap(func, func_args), total=len(func_args)))
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


def split_barcode_into_kmers(bc, k):
    kmers = []
    for i in range(0, len(bc) - k + 1):
        kmer = bc[i : i + k]
        kmers.append(kmer)
    return kmers


def load_whitelist(whitelist, k=5):
    wl = []
    with open(whitelist) as file:
        for line in file:
            bc = line.strip().split("-")[0]
            wl.append(bc)

    wl.sort()
    kmer_to_bc_index = {}
    for index, bc in enumerate(wl):
        bc_kmers = split_barcode_into_kmers(bc, k)
        for bc_kmer in bc_kmers:
            if bc_kmer not in kmer_to_bc_index.keys():
                kmer_to_bc_index[bc_kmer] = set([index])
            else:
                kmer_to_bc_index[bc_kmer].add(index)
    return wl, kmer_to_bc_index


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


def main(args):
    init_logger(args)
    n_reads, chroms = get_bam_info(args.bam)

    # Create temporary directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)

    # Process BAM alignments from each chrom separately
    func_args = []
    for chrom in chroms:
        func_args.append((args.bam, chrom, args))

    chrom_bam_fns = launch_pool(process_bam_records, func_args, args.threads)

    logger.info(f"Writing BAM with CR, CB, CY, UR, and UY tags to {args.output}")
    tmp_bam = tempfile.NamedTemporaryFile(
        prefix="tmp.align.", suffix=".unsorted.bam", dir=args.tempdir, delete=False
    )
    merge_parameters = ["-f", tmp_bam.name] + list(chrom_bam_fns)
    pysam.merge(*merge_parameters)

    pysam.sort("-@", str(args.threads), "-o", args.output, tmp_bam.name)
    pysam.index(args.output)

    logger.info("Cleaning up temporary files")
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
