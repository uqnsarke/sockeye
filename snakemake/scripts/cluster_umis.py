import argparse
import itertools
import logging
import multiprocessing
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam
from editdistance import eval as edit_distance
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="BAM file of alignments with tags for gene (GN), corrected barcode \
        (CB) and uncorrected UMI (UY)",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output BAM file with new tags for corrected UMI (UB) \
        [tagged.sorted.bam]",
        type=str,
        default="tagged.sorted.bam",
    )

    parser.add_argument(
        "-i",
        "--ref_interval",
        help="Size of genomic window (bp) to assign as gene name if no gene \
        assigned by featureCounts [1000]",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
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


def breadth_first_search(node, adj_list):
    searched = set()
    queue = set()
    queue.update((node,))
    searched.update((node,))

    while len(queue) > 0:
        node = queue.pop()
        for next_node in adj_list[node]:
            if next_node not in searched:
                queue.update((next_node,))
                searched.update((next_node,))

    return searched


def get_adj_list_directional(umis, counts, threshold=2):
    """
    identify all umis within the LEVENSHTEIN distance threshold
    and where the counts of the first umi is > (2 * second umi counts)-1
    """

    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)

    return adj_list


def get_connected_components_adjacency(umis, graph, counts):
    """
    find the connected UMIs within an adjacency dictionary
    """

    # TS: TO DO: Work out why recursive function doesn't lead to same
    # final output. Then uncomment below

    # if len(graph) < 10000:
    #    self.search = breadth_first_search_recursive
    # else:
    #    self.search = breadth_first_search

    found = set()
    components = list()

    for node in sorted(graph, key=lambda x: counts[x], reverse=True):
        if node not in found:
            # component = self.search(node, graph)
            component = breadth_first_search(node, graph)
            found.update(component)
            components.append(component)
    return components


def group_directional(clusters, adj_list, counts):
    """
    return groups for directional method
    """

    observed = set()
    groups = []
    for cluster in clusters:
        if len(cluster) == 1:
            groups.append(list(cluster))
            observed.update(cluster)
        else:
            cluster = sorted(cluster, key=lambda x: counts[x], reverse=True)
            # need to remove any node which has already been observed
            temp_cluster = []
            for node in cluster:
                if node not in observed:
                    temp_cluster.append(node)
                    observed.add(node)
            groups.append(temp_cluster)

    return groups


def cluster(counts_dict, threshold=3):
    adj_list = get_adj_list_directional(counts_dict.keys(), counts_dict, threshold)
    clusters = get_connected_components_adjacency(
        counts_dict.keys(), adj_list, counts_dict
    )
    final_umis = [list(x) for x in group_directional(clusters, adj_list, counts_dict)]
    return final_umis


def create_map_to_correct_umi(cluster_list):
    my_map = {y: x[0] for x in cluster_list for y in x}
    return my_map


def correct_umis(umis):
    # logging.info(f"clustering {len(umis)} UMIs")
    counts_dict = dict(umis.value_counts())
    umi_map = create_map_to_correct_umi(cluster(counts_dict))
    return umis.replace(umi_map)


def add_tags(chrom, umis, genes, args):
    """ """
    # Open temporary output BAM file for writing
    suff = f".{chrom}.bam"
    chrom_bam = tempfile.NamedTemporaryFile(
        prefix="tmp.align.", suffix=suff, dir=args.tempdir, delete=False
    )

    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_out = pysam.AlignmentFile(chrom_bam.name, "wb", template=bam)

    for align in bam.fetch(chrom):
        read_id = align.query_name

        # Corrected UMI = UB:Z
        align.set_tag("UB", umis[read_id], value_type="Z")

        # Annotated gene name = GN:Z
        align.set_tag("GN", genes[read_id], value_type="Z")

        bam_out.write(align)

    bam.close()
    bam_out.close()

    return chrom_bam.name


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
    n_aligns = int(sum([contig.mapped for contig in stats]))
    chroms = dict(
        [(contig.contig, contig.mapped) for contig in stats if contig.mapped > 0]
    )
    bam.close()
    return n_aligns, chroms


def create_region_name(align, args):
    """
    Create a 'gene name' based on the aligned chromosome and coordinates.
    The midpoint of the alignment determines which genomic interval to use
    for the 'gene name'.

    :param align: pysam BAM alignment
    :type align: class 'pysam.libcalignedsegment.AlignedSegment'
    :param args: object containing all supplied arguments
    :type args: class 'argparse.Namespace'
    :return: Newly created 'gene name' based on aligned chromosome and coords
    :rtype: str
    """
    chrom = align.reference_name
    start_pos = align.get_reference_positions()[0]
    end_pos = align.get_reference_positions()[-1]

    # Find the midpoint of the alignment
    midpoint = int((start_pos + end_pos) / 2)

    # Pick the genomic interval based on this alignment midpoint
    interval_start = np.floor(midpoint / args.ref_interval) * args.ref_interval
    interval_end = np.ceil(midpoint / args.ref_interval) * args.ref_interval

    # New 'gene name' will be <chr>_<interval_start>_<interval_end>
    gene = f"{chrom}_{int(interval_start)}_{int(interval_end)}"
    return gene


def read_bam(bam):
    """
    Read gene, cell barcode and UMI values from BAM entries
    """
    n_aligns, chroms = get_bam_info(bam)

    bam = pysam.AlignmentFile(bam, "rb")

    records = []
    for align in tqdm(bam.fetch(), total=n_aligns):
        read_id = align.query_name

        # Corrected cell barcode = CB:Z
        bc_corr = align.get_tag("CB")

        # Uncorrected UMI = UR:Z
        umi_uncorr = align.get_tag("UR")

        # Annotated gene name = GN:Z
        gene = align.get_tag("GN")

        if gene == "NA":
            # Group by region if no gene annotation
            gene = create_region_name(align, args)

        records.append((read_id, gene, bc_corr, umi_uncorr))

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "bc", "umi_uncorr"]
    )
    bam.close()

    return df


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


def launch_df_pool(df_grouped, func, threads, chunks):
    """ """
    p = multiprocessing.Pool(processes=threads)
    try:
        res = list(tqdm(p.imap(func, [df_ for idx, df_ in df_grouped]), total=chunks))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return pd.concat(res)


def func_group_apply(df):
    return df.groupby(["gene", "bc"])["umi_uncorr"].apply(correct_umis)


def process_bam_records(tup):
    """
    Read through all the alignments for specific chromosome to pull out
    the gene, barcode, and uncorrected UMI information. Use that to cluster
    UMIs and get the corrected UMI sequence. We'll then write out a temporary
    chromosome-specific BAM file with the corrected UMI tag (UB).

    :param tup: Tuple containing the input arguments
    :type tup: tup
    :return: Path to a temporary BAM file
    :rtype: str
    """
    input_bam = tup[0]
    chrom = tup[1]
    args = tup[2]

    # Open input BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    records = []
    for align in bam.fetch(contig=chrom):
        read_id = align.query_name

        # Make sure each alignment in this BAM has the expected tags
        for tag in ["UR", "CB", "GN"]:
            assert align.has_tag(tag), f"{tag} tag not found in f{read_id}"

        # Corrected cell barcode = CB:Z
        bc_corr = align.get_tag("CB")

        # Uncorrected UMI = UR:Z
        umi_uncorr = align.get_tag("UR")

        # Annotated gene name = GN:Z
        gene = align.get_tag("GN")

        # If no gene annotation exists
        if gene == "NA":
            # Group by region if no gene annotation
            gene = create_region_name(align, args)

        records.append((read_id, gene, bc_corr, umi_uncorr))

    # Done reading the chrom-specific alignments from the BAM
    bam.close()

    # Create a dataframe with chrom-specific data
    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "bc", "umi_uncorr"]
    )

    df["umi_corr"] = df.groupby(["gene", "bc"])["umi_uncorr"].apply(correct_umis)

    # Simplify to a read_id:umi_corr dictionary
    df = df.drop(["bc", "umi_uncorr"], axis=1).set_index("read_id")

    # Dict of corrected UMI for each read ID
    umis = df.to_dict()["umi_corr"]

    # Dict of gene names to add <chr>_<start>_<end> in place of NA
    genes = df.to_dict()["gene"]

    # Add corrected UMIs to each chrom-specific BAM entry via the UB:Z tag
    chrom_bam = add_tags(chrom, umis, genes, args)

    return chrom_bam


def main(args):
    init_logger(args)
    n_aligns, chroms = get_bam_info(args.bam)

    # Create temporary directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir)
    os.mkdir(args.tempdir)

    logger.info(f"Processing input BAM using {args.threads} threads")
    func_args = []
    chroms_sorted = dict(sorted(chroms.items(), key=lambda item: item[1]))
    for chrom in chroms_sorted.keys():
        func_args.append((args.bam, chrom, args))

    chrom_bam_fns = launch_pool(process_bam_records, func_args, args.threads)

    logger.info(f"Writing BAM with CR, CB, CY, UR, UB, and UY tags to {args.output}")
    tmp_bam = tempfile.NamedTemporaryFile(
        prefix="tmp.align.", suffix=".unsorted.bam", dir=args.tempdir, delete=False
    )
    merge_parameters = ["-f", tmp_bam.name] + list(chrom_bam_fns)
    pysam.merge(*merge_parameters)

    pysam.sort("-@", str(args.threads), "-o", args.output, tmp_bam.name)

    logger.info("Cleaning up temporary files")
    shutil.rmtree(args.tempdir)


if __name__ == "__main__":
    args = parse_args()

    main(args)
