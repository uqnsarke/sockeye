import argparse
import itertools
import logging
import subprocess

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


def add_tags(umis, genes, n_aligns, args):
    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_out = pysam.AlignmentFile(args.output, "wb", template=bam)

    for i, align in enumerate(tqdm(bam.fetch(), total=n_aligns)):
        read_id = align.query_name

        # Corrected UMI = UB:Z
        align.set_tag("UB", umis[read_id], value_type="Z")

        # Annotated gene name = GN:Z
        align.set_tag("GN", genes[read_id], value_type="Z")

        bam_out.write(align)

    bam.close()
    bam_out.close()


def count_bam_entries(bam):
    """
    Call `samtools idxstat <BAM> | cut -f3` to get number of alignments
    """
    stdout, stderr = run_subprocess(f"samtools idxstats {bam} | cut -f3")
    n_aligns = np.sum([int(x) for x in stdout.split("\n") if x != ""])
    return n_aligns


def read_bam(args):
    """
    Read gene, cell barcode and UMI values from BAM entries
    """
    logger.info(f"Counting alignments in {args.bam}")
    n_aligns = count_bam_entries(args.bam)

    bam = pysam.AlignmentFile(args.bam, "rb")

    logger.info(f"Reading read / barcode / UMI info from {args.bam}")
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

        records.append((read_id, gene, bc_corr, umi_uncorr))

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "bc", "umi_uncorr"]
    )
    bam.close()

    return df, n_aligns


def main(args):
    init_logger(args)

    df, n_aligns = read_bam(args)

    logger.info("Correcting UMIs: grouping by gene and corrected barcode")
    # Cluster UMIs that have been grouped by gene and cell barcode
    df["umi_corr"] = df.groupby(["gene", "bc"])["umi_uncorr"].apply(correct_umis)

    # Simplify to a read_id:umi_corr dictionary
    df = df.drop(["bc", "umi_uncorr"], axis=1).set_index("read_id")

    # Dict of corrected UMI for each read ID
    umis = df.to_dict()["umi_corr"]
    # Dict of gene names to add <chr>_<start>_<end> in place of NA
    genes = df.to_dict()["gene"]

    logger.info(f"Writing corrected UMI tags (UB) into {args.output}")
    # Add corrected UMIs to each BAM entry via the UB:Z tag
    add_tags(umis, genes, n_aligns, args)
    pysam.index(args.output)


if __name__ == "__main__":
    args = parse_args()

    main(args)
