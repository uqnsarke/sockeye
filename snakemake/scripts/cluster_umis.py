import argparse
import itertools
import logging
import subprocess

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
        help="BAM file of alignments with tags for gene, corrected barcode and \
        uncorrected UMI",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output BAM file with new tags for corrected UMI (UB:Z) [tagged.sorted.bam]",
        type=str,
        default="tagged.sorted.bam",
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


def add_umi_corr_tag(umis, n_aligns, args):
    bam = pysam.AlignmentFile(args.bam, "rb")
    bam_tagged = pysam.AlignmentFile(args.output, "wb", template=bam)

    logger.info(f"Writing corrected UMI tags (UB:Z) into {args.output}")
    for i, align in enumerate(tqdm(bam.fetch(), total=n_aligns)):
        read_id = align.query_name
        if umis.get(read_id):
            # Read has an annotated gene, barcode and corrected UMI
            align.set_tag("UB", umis[read_id], value_type="Z")
        else:
            # Read probably has no annotated gene (gene name = NA)
            align.set_tag("UB", "XXXXXXXXX", value_type="Z")

        bam_tagged.write(align)

    bam.close()
    bam_tagged.close()


def read_bam(args):
    """
    Read gene, cell barcode and UMI values from BAM entries
    """
    logger.info(f"Counting alignments in {args.bam}")
    stdout, stderr = run_subprocess(f"samtools view -c -F 260 {args.bam}")
    n_aligns = int(stdout.strip())

    bam = pysam.AlignmentFile(args.bam, "rb")

    logger.info(f"Reading read / barcode / UMI info from {args.bam}")
    records = []
    for i, align in enumerate(tqdm(bam.fetch(), total=n_aligns)):
        # Annotated gene name = GN:Z
        gene = align.get_tag("GN")
        if gene != "NA":
            read_id = align.query_name
            # Corrected cell barcode = CB:Z
            bc_corr = align.get_tag("CB")
            # Uncorrected UMI = UR:Z
            umi_uncorr = align.get_tag("UR")
            # UMI quality score = UY:Z
            # umi_qv = align.get_tag("UY")

            records.append((read_id, gene, bc_corr, umi_uncorr))

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "bc", "umi_uncorr"]
    )

    bam.close()

    return df, n_aligns


def main(args):
    init_logger(args)

    df, n_aligns = read_bam(args)

    # Cluster UMIs that have been grouped by gene and cell barcode
    df["umi_corr"] = df.groupby(["gene", "bc"])["umi_uncorr"].apply(correct_umis)

    # Simplify to a read_id:umi_corr dictionary
    df = df.drop(["gene", "bc", "umi_uncorr"], axis=1).set_index("read_id")
    umis = df.to_dict()["umi_corr"]

    # Add corrected UMIs to each BAM entry via the UB:Z tag
    add_umi_corr_tag(umis, n_aligns, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
