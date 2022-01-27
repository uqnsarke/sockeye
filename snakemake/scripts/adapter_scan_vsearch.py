import argparse
import gzip
import logging
import multiprocessing
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
from glob import glob

import numpy as np
import pandas as pd
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("fastq", help="FASTQ of ONT reads", type=str)

    # Optional arguments
    parser.add_argument(
        "--output_fastq",
        help="Output file name for stranded FASTQ entries \
                        [stranded.fastq]",
        type=str,
        default="stranded.fastq",
    )

    parser.add_argument(
        "--output_tsv",
        help="Output file name for adapter configurations \
                        [adapters.tsv]",
        type=str,
        default="adapters.tsv",
    )

    parser.add_argument(
        "--output_vsearch",
        help="If specified, write VSEARCH output to \
                        <adapters.vsearch.tsv> [None]",
        type=str,
        default=None,
    )

    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "-b",
        "--batch_size",
        help="Number of reads per batch [100000]",
        type=int,
        default=100000,
    )

    parser.add_argument(
        "-r1",
        "--read1",
        help="Read1 adapter sequence \
                        [CTACACGACGCTCTTCCGATCT]",
        type=str,
        default="CTACACGACGCTCTTCCGATCT",
    )

    parser.add_argument(
        "-s",
        "--TSO",
        help="TSO sequence [ATGTACTCTGCGTTGATACCACTGCTT]",
        type=str,
        default="ATGTACTCTGCGTTGATACCACTGCTT",
    )

    parser.add_argument(
        "-i",
        "--min_adapter_id",
        help="Minimum adapter alignment identity for VSEARCH \
                        [0.7]",
        type=float,
        default=0.7,
    )

    parser.add_argument(
        "--only_strand_full_length",
        help="Do not try to strand-orient reads where either \
                        just a single read1 or single TSO adapter was found \
                        [False]",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "-a",
        "--adapters_fasta",
        help="Filename for adapter query sequences \
                        [read1_TSO.fasta]",
        type=str,
        default="read1_TSO.fasta",
    )

    parser.add_argument(
        "--no_fastq",
        help="Do not write the stranded fastqs \
                        (specified by --output_fastq) [False]",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "--dolomite",
        help="Reads are from Dolomite instead of 10x [False]",
        action="store_true",
        default=False,
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
    p = pathlib.Path(args.output_tsv)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    # Add tempdir path to the adapters_fasta
    args.adapters_fasta = os.path.join(args.tempdir, args.adapters_fasta)

    return args


# complement translation table with support for regex punctuation
COMPLEMENT_TRANS = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn"
)


def run_subprocess(cmd):
    """
    Run OS command and return stdout & stderr
    """
    p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    return str(stdout), str(stderr)


def check_vsearch():
    stdout, stderr = run_subprocess("vsearch --quiet -h")
    if stderr.find("vsearch: command not found") > -1:
        logging.error("Could not load find VSEARCH -- check installation")
        sys.exit(1)


def write_adapters_fasta(args):
    adapters = []
    for adapter, seq in {
        "read1_f": args.read1,
        "read1_r": args.read1[::-1].translate(COMPLEMENT_TRANS),
        "TSO_f": args.TSO,
        "TSO_r": args.TSO[::-1].translate(COMPLEMENT_TRANS),
    }.items():
        entry = SeqRecord(Seq(seq), id=adapter, name="", description="")

        adapters.append(entry)
    SeqIO.write(adapters, args.adapters_fasta, "fasta")


def write_tmp_fasta(batch_reads, args):
    tmp_fasta = tempfile.NamedTemporaryFile(
        prefix="tmp.reads.", suffix=".fasta", dir=args.tempdir, delete=False
    )

    with open(tmp_fasta.name, mode="w") as f_out:
        for r in batch_reads:
            f_out.write(">" + str(r.name) + "\n")
            f_out.write(str(r.sequence) + "\n")
    return tmp_fasta.name


def call_vsearch(tmp_fasta, args):
    """ """
    tmp_vsearch = tmp_fasta.replace(".fasta", ".vsearch.tsv")

    vsearch_cmd = "vsearch --usearch_global {fasta} --db {adapters} \
    --threads 1 --minseqlength 20 --maxaccepts 5 --id {id} --strand plus \
    --wordlength 3 --minwordmatches 10 --output_no_hits --userfields \
    'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
    --userout {output}".format(
        fasta=tmp_fasta,
        id=args.min_adapter_id,
        adapters=args.adapters_fasta,
        output=tmp_vsearch,
    )
    stdout, stderr = run_subprocess(vsearch_cmd)
    os.remove(tmp_fasta)
    return tmp_vsearch


def get_valid_adapter_pair_positions_in_read(read):
    valid_pairs_n = 0
    fl_pairs = []

    def write_valid_pair_dict(read, adapter_1_idx, adapter_2_idx, valid_pairs_n):
        read_id = read["query"].iloc[0]
        pair_str = (
            f"{read.iloc[adapter_1_idx]['target']}-{read.iloc[adapter_2_idx]['target']}"
        )
        fl_pair = {
            "read_id": "{}_{}".format(read_id, valid_pairs_n),
            "config": pair_str,
            "start": read.iloc[adapter_1_idx]["qilo"],
            "end": read.iloc[adapter_2_idx]["qihi"],
        }
        return fl_pair, valid_pairs_n

    compat_adapters = {"read1_f": "TSO_f", "TSO_r": "read1_r"}

    # First find all instances of read1_f
    for adapter1 in ["read1_f", "TSO_r"]:
        adapter_1_idxs = read.index[read["target"] == adapter1]
        for adapter_1_idx in adapter_1_idxs:
            # For each found read1_f, examine next found adapter
            adapter_2_idx = adapter_1_idx + 1
            # Make sure there are enough alignments to allow this indexing
            if adapter_2_idx < read.shape[0]:
                # Is the next found adapter a TSO_f?
                if read.iloc[adapter_2_idx]["target"] == compat_adapters[adapter1]:
                    # This is a valid adapter pairing (read1_f-TSO_f)
                    fl_pair, valid_pairs_n = write_valid_pair_dict(
                        read, adapter_1_idx, adapter_2_idx, valid_pairs_n
                    )
                    if adapter1 == "read1_f":
                        fl_pair["strand"] = "+"
                    else:
                        fl_pair["strand"] = "-"
                    valid_pairs_n += 1
                    fl_pairs.append(fl_pair)
    return fl_pairs


def filter_dolomite(df):
    """
    The read1/TSO adapter alignments are overlapping due
    to their largely complementary sequences. We need
    to identify such overlapping alignments and keep
    only the highest quality alignment.
    """
    to_concat = []
    for read_id, df_ in df.groupby("query"):
        for field in ["qilo", "qihi"]:
            df_ = df_.sort_values([field, "id"]).drop_duplicates(
                subset=field, keep="last"
            )
        to_concat.append(df_)

    df = pd.concat([df_ for df_ in to_concat], axis=0)
    return df


def add_entry_to_read_info(
    read_info,
    orig_read_id,
    new_read_id,
    readlen,
    start,
    end,
    fl,
    stranded,
    strand,
    orig_adapter_config,
    adapter_config,
    lab,
):
    read_info[orig_read_id][new_read_id] = {
        "readlen": readlen,
        "start": start,
        "end": end,
        "fl": fl,
        "stranded": stranded,
        "orig_strand": strand,
        "orig_adapter_config": orig_adapter_config,
        "adapter_config": adapter_config,
        "lab": lab,
    }
    return read_info


def parse_vsearch(tmp_vsearch, args):
    colnames = [
        "query",
        "target",
        "id",
        "alnlen",
        "mism",
        "opens",
        "qilo",
        "qihi",
        "qstrand",
        "tilo",
        "tihi",
        "ql",
        "tl",
    ]

    df = pd.read_csv(tmp_vsearch, sep="\t", header=None, names=colnames)

    if args.dolomite:
        df = filter_dolomite(df)

    read_info = {}
    for read_id, read in df.groupby("query"):
        # Sort aligned adapters by their position in the read
        read = read.sort_values("qilo").reset_index()

        # Get the original adapter config
        orig_adapter_config = "-".join(read["target"])

        orig_read_id = read["query"].iloc[0]
        read_info[orig_read_id] = {}

        # Look for reads with valid consecutive features.
        # For example, a read1_f followed immediately by
        # a TSO_f, or a TSO_r followed immediately by
        # a read1_r.
        fl_pairs = get_valid_adapter_pair_positions_in_read(read)
        if len(fl_pairs) > 0:
            for fl_pair in fl_pairs:
                fl = True
                stranded = True
                read_id = fl_pair["read_id"]
                strand = fl_pair["strand"]
                readlen = fl_pair["end"] - fl_pair["start"]
                start = fl_pair["start"]
                end = fl_pair["end"]
                adapter_config = fl_pair["config"]
                lab = "full_len"
                read_info = add_entry_to_read_info(
                    read_info,
                    orig_read_id,
                    read_id,
                    readlen,
                    start,
                    end,
                    fl,
                    stranded,
                    strand,
                    orig_adapter_config,
                    adapter_config,
                    lab,
                )
        else:
            # No valid adapter pairs found. Either single adapter,
            # no adapter, or weird artifact read
            read_id = "{}_0".format(read["query"].iloc[0])
            fl = False
            stranded = False
            strand = "*"
            readlen = read["ql"].iloc[0]
            start = 0
            end = readlen - 1
            adapter_config = "-".join(list(read["target"].values))
            if adapter_config in ["TSO_r-TSO_f", "TSO_f-TSO_r"]:
                lab = "double_tso"
            elif adapter_config in ["read1_r-read1_f", "read1_f-read1_r"]:
                lab = "double_read1"
            elif adapter_config in ["TSO_f", "TSO_r"]:
                lab = "single_tso"
                if not args.only_strand_full_length:
                    # We want to strand and trim reads where we only have
                    # a TSO sequence but no read1. These MIGHT contain
                    # the cell barcode, UMI, polyT and cDNA sequence, but
                    # since the read1 is low-quality, these might be of
                    # of dubious value. We can only trim the TSO end of the
                    # read based on TSO aligned positions, but won't touch
                    # the putative read1 end.
                    stranded = True
                    if adapter_config == "TSO_f":
                        strand = "+"
                        start = 0
                        end = read.iloc[0]["qihi"]
                        readlen = end - start
                    elif adapter_config == "TSO_r":
                        strand = "-"
                        start = read.iloc[0]["qilo"]
                        end = read.iloc[0]["ql"] - 1
                        readlen = end - start
                    else:
                        raise Exception("Shouldn't be here!")
            elif adapter_config in ["read1_f", "read1_r"]:
                lab = "single_read1"
                if not args.only_strand_full_length:
                    # We want to strand and trim reads where we only have
                    # a read1 sequence but no TSO. These should contain
                    # the cell barcode, UMI, polyT and cDNA sequence, so
                    # should still have significant value. We can only trim
                    # the read1 end of the read, but won't touch the putative
                    # TSO end.
                    stranded = True
                    if adapter_config == "read1_f":
                        strand = "+"
                        start = read.iloc[0]["qilo"]
                        end = read.iloc[0]["ql"] - 1
                        readlen = end - start
                    elif adapter_config == "read1_r":
                        strand = "-"
                        start = 0
                        end = read.iloc[0]["qihi"]
                        readlen = end - start
                    else:
                        raise Exception("Shouldn't be here!")
            elif adapter_config == "*":
                lab = "no_adapters"
            else:
                lab = "other"
            read_info = add_entry_to_read_info(
                read_info,
                orig_read_id,
                read_id,
                readlen,
                start,
                end,
                fl,
                stranded,
                strand,
                orig_adapter_config,
                adapter_config,
                lab,
            )
    return read_info, colnames


def revcomp_adapter_config(adapters_string):
    d = {
        "read1_f": "read1_r",
        "read1_r": "read1_f",
        "TSO_f": "TSO_r",
        "TSO_r": "TSO_f",
    }

    "TSO_r-read1_r"
    rc_string = "-".join([d[a] for a in adapters_string.split("-")[::-1]])
    return rc_string


def write_stranded_fastq(batch_reads, read_info, args):
    """ """
    tmp_fastq = tempfile.NamedTemporaryFile(
        prefix="tmp.stranded.", suffix=".fastq", dir=args.tempdir, delete=False
    )

    with open(tmp_fastq.name, mode="w") as f_out:
        for r in batch_reads:
            read_id = r.name.split(" ")[0]
            for subread_id in read_info[read_id].keys():
                d = read_info[read_id][subread_id]
                # subread_info = read_info[read_id][subread_id]
                subread_seq = str(r.sequence[d["start"] : d["end"]])
                subread_quals = r.quality[d["start"] : d["end"]]
                if d["orig_strand"] == "-":
                    rc_config = revcomp_adapter_config(d["adapter_config"])
                    d["adapter_config"] = rc_config
                    subread_seq = subread_seq[::-1].translate(COMPLEMENT_TRANS)
                    subread_quals = subread_quals[::-1]
                f_out.write(f"@{subread_id}\n")
                f_out.write(f"{subread_seq}\n")
                f_out.write("+\n")
                f_out.write(f"{subread_quals}\n")

    return tmp_fastq.name


def open_fastq(fastq):
    if args.fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq)
    return f


def count_reads(fastq):
    logging.info("Counting reads")
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in tqdm(f, unit_scale=0.25, unit=" reads"):
            number_lines += 1
    return number_lines / 4


def batch_iterator(iterator, args):
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
        while len(batch) < args.batch_size:
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch, args


def get_subread_info(read_info):
    subread_info = []
    for read_id, subread_d in read_info.items():
        for subread_id, attr_d in subread_d.items():
            attr_d["read_id"] = subread_id
            subread_info.append(attr_d)
    return subread_info


def write_tmp_table(tmp_fastq, subread_info):
    df = pd.DataFrame.from_records(subread_info)
    tmp_table = tmp_fastq.replace(".fastq", ".info.tsv")
    df.to_csv(tmp_table, sep="\t", index=False)
    return tmp_table


def process_batch(tup):
    batch_reads = tup[0]
    args = tup[1]

    # VSEARCH needs a FASTA
    tmp_fasta = write_tmp_fasta(batch_reads, args)

    tmp_vsearch = call_vsearch(tmp_fasta, args)
    read_info, vsearch_cols = parse_vsearch(tmp_vsearch, args)
    if not args.no_fastq:
        stranded_tmp_fastq = write_stranded_fastq(batch_reads, read_info, args)
    else:
        stranded_tmp_fastq = None
    subread_info = get_subread_info(read_info)
    tmp_table = write_tmp_table(stranded_tmp_fastq, subread_info)
    return stranded_tmp_fastq, tmp_table, vsearch_cols


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

    # If specified batch size is > total number of reads, reduce batch size
    n_reads = count_reads(args.fastq)
    args.batch_size = min(n_reads, args.batch_size)
    n_batches = int(np.ceil(n_reads / args.batch_size))

    # Open pysam handle to input FASTQ
    fastq = pysam.FastxFile(args.fastq)

    # Create temp directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)

    # Write a FASTA file containing the adapter sequences for VSEARCH to use
    write_adapters_fasta(args)

    logging.info("Processing {} batches of {} reads".format(n_batches, args.batch_size))
    tmp_tables = []
    with multiprocessing.Pool(args.threads) as p:
        r = list(
            tqdm(
                p.imap(process_batch, batch_iterator(fastq, args)),
                total=n_batches,
            )
        )
    tmp_fastqs, tmp_tables, vsearch_cols = zip(*r)
    vsearch_cols = vsearch_cols[0]

    # Merge temp tables and fastqs then clean up
    logging.info(f"Writing output table to {args.output_tsv}")
    if len(tmp_tables) > 1:
        pd.concat([pd.read_csv(d, sep="\t") for d in tmp_tables], axis=0).to_csv(
            args.output_tsv, sep="\t", index=False
        )
    else:
        shutil.copy(tmp_tables[0], args.output_tsv)

    if not args.no_fastq:
        logging.debug(f"Writing stranded fastq to {args.output_fastq}")
        with open(args.output_fastq, "wb") as f_out:
            for tmp_fastq in tmp_fastqs:
                with open(tmp_fastq, "rb") as f_:
                    shutil.copyfileobj(f_, f_out)

    if args.output_vsearch is not None:
        # Merge temp VSEARCH tables
        logging.debug(f"Writing VSEARCH output to {args.output_vsearch}")
        glob_str = os.path.join(args.tempdir, "*.vsearch.tsv")
        pd.concat(
            [pd.read_csv(d, sep="\t", header=None) for d in glob(glob_str)], axis=0
        ).to_csv(args.output_vsearch, sep="\t", index=False, header=vsearch_cols)

    logging.debug("Cleaning up")
    os.remove(args.adapters_fasta)
    [os.remove(fn) for fn in tmp_tables]
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
