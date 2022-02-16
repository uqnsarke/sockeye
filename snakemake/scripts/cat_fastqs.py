import argparse
import collections
import gzip
import logging
import multiprocessing
import os
import pathlib
import shutil

import numpy as np
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "fofn",
        help="File of filenames containing paths to all input FASTQs",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "-t",
        "--threads",
        help="Threads to use. This is also the number of \
        chunks that will be output [4]",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--output_dir",
        help="Output directory for chunked FASTQ files [./chunked_fastqs]",
        type=str,
        default="./chunked_fastqs",
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


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


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


def load_fofn(fn):
    """ """
    input_fastqs = []
    for line in open(fn, "r"):
        path = line.strip()
        input_fastqs.append(path)
    return input_fastqs


def open_fastq(fastq):
    if fastq.split(".")[-1] == "gz":
        f = gzip.open(fastq, "rt")
    else:
        f = open(fastq)
    return f


def count_reads(fastq):
    number_lines = 0
    with open_fastq(fastq) as f:
        for line in tqdm(f, unit_scale=0.25, unit=" reads"):
            number_lines += 1
    f.close()
    return int(number_lines / 4)


def cat_files(tup):
    """ """
    chunk_fns = tup[0]
    i = tup[1]
    ext = tup[2]
    args = tup[3]

    out_fn = os.path.join(args.output_dir, f"proc.{i}.{ext}")
    if ext.split(".")[-1] == "gz":
        with open(out_fn, "wb") as outfile:
            for fname in chunk_fns:
                with open(fname, "rb") as infile:
                    outfile.write(infile.read())
    else:
        with open(out_fn, "w") as outfile:
            for fname in chunk_fns:
                with open(fname, "r") as infile:
                    outfile.write(infile.read())


def get_input_file_ext(input_fastqs, args):
    """
    Determine the extension of the FASTQ files listed in
    the input fofn. We will maintain this extension in
    the catted output files.

    :param input: List of filenames that will be catted together
    :type input: list
    :param args: object containing all supplied arguments
    :type args: class 'argparse.Namespace'
    :return: Extension (e.g. filename.<ext>)
    :rtype: str
    """
    if len(input_fastqs) > 1:
        ext = set([p.split(".")[-1] for p in input_fastqs])
        assert len(ext) == 1, f"Unexpected mixture of file extensions in {args.fofn}"

        ext = list(ext)[0]

        if ext == "gz":
            subext = set([p.split(".")[-2] for p in input_fastqs])
            assert (
                len(subext) == 1
            ), f"Unexpected mixture of file extensions in {args.fofn}"

            subext = list(subext)[0]
            ext = f"{subext}.{ext}"

    elif len(input_fastqs) == 1:
        ext = input_fastqs[0].split(".")[-1]
        if ext == "gz":
            subext = input_fastqs[0].split(".")[-2]
            ext = f"{subext}.{ext}"

    return ext


def main(args):
    init_logger(args)

    input_fastqs = load_fofn(args.fofn)

    # Determine file extension (e.g. ".fastq.gz", ".fastq", ".fq.gz", or ".fq")
    ext = get_input_file_ext(input_fastqs, args)

    # Create output directory
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.mkdir(args.output_dir)

    if len(input_fastqs) > args.threads:
        # There are many batched fastqs, so combine them into N=<args.threads>
        # files for efficient downstream processing. Makes use of multiprocessing
        chunk_size = int((len(input_fastqs) / args.threads)) + 1
        func_args = []
        for i, chunk_fns in enumerate(chunks(input_fastqs, chunk_size)):
            func_args.append((chunk_fns, i + 1, ext, args))

        p = multiprocessing.Pool(processes=args.threads)
        try:
            p.imap(cat_files, func_args)
            p.close()
            p.join()
        except KeyboardInterrupt:
            p.terminate()

    elif len(input_fastqs) < args.threads:
        # There are less input files than available threads, so split them up
        # into N=<args.threads> files for efficient downstream processing
        # 1. Cat the files together
        # 2. Count the reads
        # 3. Split combined file into N=args.threads chunks
        if len(input_fastqs) == 1:
            tmp_combined_fn = input_fastqs[0]
        else:
            tmp_combined_fn = os.path.join(args.output_dir, "tmp.combined.fastq")
            if ext.split(".")[-1] == "gz":
                with open(tmp_combined_fn, "wb") as outfile:
                    for fname in input_fastqs:
                        with gzip.open(fname, "rb") as infile:
                            for line in infile:
                                outfile.write(line)
            else:
                with open(tmp_combined_fn, "w") as outfile:
                    for fname in input_fastqs:
                        with open(fname, "r") as infile:
                            for line in infile:
                                outfile.write(line)

        logger.info("Counting total reads")
        n_reads = count_reads(tmp_combined_fn)
        chunk_size = int(np.ceil(n_reads / args.threads))

        read_iterator = pysam.FastxFile(tmp_combined_fn)

        for i, chunk_reads in enumerate(batch_iterator(read_iterator, chunk_size)):
            out_fn = os.path.join(args.output_dir, f"proc.{i}.fastq")
            with open(out_fn, "w") as outfile:
                for read in chunk_reads:
                    outfile.write(f"@{read.name}\n")
                    outfile.write(f"{read.sequence}\n")
                    outfile.write("+\n")
                    outfile.write(f"{read.quality}\n")
        os.remove(tmp_combined_fn)

    else:
        # There are already N=<args.threads> input files, so just copy them into
        # the output directory for downstream processing
        for i, fn in enumerate(input_fastqs):
            dest_fn = os.path.join(args.output_dir, f"proc.{i}.{ext}")
            shutil.copy(fn, dest_fn)


if __name__ == "__main__":
    args = parse_args()

    main(args)
