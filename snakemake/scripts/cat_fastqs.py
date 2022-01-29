import argparse
import collections
import logging
import multiprocessing
import os
import pathlib
import shutil

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


def load_fofn(fn):
    """ """
    input_fastqs = []
    for line in open(fn, "r"):
        path = line.strip()
        input_fastqs.append(path)
    return input_fastqs


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
    ext = set([p.split(".")[-1] for p in input_fastqs])
    assert len(ext) == 1, f"Unexpected mixture of file extensions in {args.fofn}"

    ext = list(ext)[0]

    if ext == "gz":
        subext = set([p.split(".")[-2] for p in input_fastqs])
        assert len(subext) == 1, f"Unexpected mixture of file extensions in {args.fofn}"

        subext = list(subext)[0]
        ext = f"{subext}.{ext}"

    return ext


def main(args):
    input_fastqs = load_fofn(args.fofn)

    # Determine file extension (e.g. ".fastq.gz", ".fastq", ".fq.gz", or ".fq")
    ext = get_input_file_ext(input_fastqs, args)

    # Create output directory
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.mkdir(args.output_dir)

    chunk_size = int((len(input_fastqs) / args.threads)) + 1
    func_args = []
    for i, chunk_fns in enumerate(chunks(input_fastqs, chunk_size)):
        func_args.append((chunk_fns, i + 1, ext, args))

    p = multiprocessing.Pool(processes=args.threads)
    try:
        results = p.imap(cat_files, func_args)
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


if __name__ == "__main__":
    args = parse_args()

    main(args)
