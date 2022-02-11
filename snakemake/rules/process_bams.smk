rule extract_barcodes:
    input:
        bam=BAM_SORT,
    output:
        bam=BAM_BC_UNCORR_TMP,
        tsv=BARCODE_COUNTS,
    params:
        barcodes=config["BC_SUPERLIST"],
        read1=config["READ_STRUCTURE"]["READ1"],
        read1_suff_length=config["BARCODE"]["READ1_SUFF_LENGTH"],
        barcode_length=config["READ_STRUCTURE"]["BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE"]["UMI_LENGTH"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/extract_barcode.py "
        "-t {threads} "
        "--read1_adapter {params.read1} "
        "--read1_suff_length {params.read1_suff_length} "
        "--barcode_length {params.barcode_length} "
        "--umi_length {params.umi_length} "
        "--output_bam {output.bam} "
        "--output_barcodes {output.tsv} "
        "{input.bam} {params.barcodes}"


rule cleanup_headers_1:
    input:
        BAM_BC_UNCORR_TMP,
    output:
        # bam=temp(BAM_BC_UNCORR),
        # bai=temp(BAM_BC_UNCORR_BAI),
        bam=BAM_BC_UNCORR,
        bai=BAM_BC_UNCORR_BAI,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools reheader --no-PG -c 'grep -v ^@PG' "
        "{input} > {output.bam}; "
        "samtools index {output.bam}"


rule generate_whitelist:
    input:
        BARCODE_COUNTS,
    output:
        whitelist=BARCODE_WHITELIST,
        kneeplot=BARCODE_KNEEPLOT,
    params:
        flags=config["BARCODE"]["KNEEPLOT_FLAGS"],
    conda:
        "../envs/kneeplot.yml"
    shell:
        "python {SCRIPT_DIR}/knee_plot.py "
        "{params.flags} "
        "--output_whitelist {output.whitelist} "
        "--output_plot {output.kneeplot} "
        "{input}"


checkpoint split_bam_by_chroms:
    input:
        bam=BAM_BC_UNCORR,
        bai=BAM_BC_UNCORR_BAI,
    output:
        split_dir=directory(SPLIT_DIR),
    conda:
        "../envs/barcodes.yml"
    threads: config["MAX_THREADS"]
    shell:
        "python {SCRIPT_DIR}/split_bam_by_chroms.py "
        "-t {threads} "
        "--output_dir {output.split_dir} "
        "{input.bam}"


rule assign_barcodes:
    input:
        bam=CHROM_BAM_BC_UNCORR,
        bai=CHROM_BAM_BC_UNCORR_BAI,
        whitelist=BARCODE_WHITELIST,
    output:
        bam=CHROM_BAM_BC_TMP,
        counts=CHROM_ASSIGNED_BARCODE_COUNTS,
    params:
        max_ed=config["BARCODE"]["MAX_ED"],
        min_ed_diff=config["BARCODE"]["MIN_ED_DIFF"],
        read1=config["READ_STRUCTURE"]["READ1"],
        read1_suff_length=config["BARCODE"]["READ1_SUFF_LENGTH"],
        barcode_length=config["READ_STRUCTURE"]["BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE"]["UMI_LENGTH"],
    threads: 1
    conda:
        "../envs/barcodes.yml"
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/assign_barcodes.py "
        "-t {threads} "
        "--output_bam {output.bam} "
        "--output_counts {output.counts} "
        "--max_ed {params.max_ed} "
        "--min_ed_diff {params.min_ed_diff} "
        "--read1_adapter {params.read1} "
        "--read1_suff_length {params.read1_suff_length} "
        "--barcode_length {params.barcode_length} "
        "--umi_length {params.umi_length} "
        "{input.bam} {input.whitelist} "


rule cleanup_headers_2:
    input:
        CHROM_BAM_BC_TMP,
    output:
        bam=CHROM_BAM_BC,
        bai=CHROM_BAM_BC_BAI,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools reheader --no-PG -c 'grep -v ^@PG' "
        "{input} > {output.bam}; "
        "samtools index {output.bam}"


# rule ont_featureCounts_genes:
#     input:
#         bam=CHROM_BAM_BC,
#         bai=CHROM_BAM_BC_BAI,
#     output:
#         featureCounts=temp(CHROM_FC_READ_ASSIGNS_TMP),
#         featureCounts_final=temp(CHROM_FC_READ_ASSIGNS),
#         gene_assigned=temp(CHROM_FC_GENES),
#         summary=temp(CHROM_FC_SUMMARY),
#     params:
#         rpath=str(SPLIT_DIR),
#         gtf=config["REF_GTF"],
#         minQV=config["FEATURECOUNTS"]["MINQV"],
#     threads: 1
#     conda:
#         "../envs/feature_counts.yml"
#     shell:
#         "featureCounts "
#         "-T {threads} "
#         "-a {params.gtf} "
#         "-g gene_name "
#         "-t gene "
#         "-o {output.gene_assigned} -L "
#         "-R CORE --primary "
#         "--Rpath {params.rpath} "
#         "-F GTF "
#         "-Q {params.minQV} "
#         "--largestOverlap "
#         "{input.bam} "
#         "&& cp {output.featureCounts} {output.featureCounts_final}"


rule bam_to_bed:
    input:
        bam=CHROM_BAM_BC,
        bai=CHROM_BAM_BC_BAI,
    output:
        bed=CHROM_BED_BC,
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools bamtobed -i {input.bam} > {output.bed}"


rule split_gtf_by_chroms:
    input:
        config["REF_GTF"],
    output:
        CHROM_GTF,
    params:
        chrom=lambda w: w.chrom,
    shell:
        "cat {input} | "
        "awk '$1==\"{params.chrom}\" {{print}}' > {output}"


rule assign_genes:
    input:
        bed=CHROM_BED_BC,
        chrom_gtf=CHROM_GTF,
    output:
        gene_assigned=CHROM_TSV_GENE_ASSIGNS,
    params:
        minQV=config["GENE_ASSIGNS"]["MINQV"],
    threads: 1
    conda:
        "../envs/assign_genes.yml"
    shell:
        "python {SCRIPT_DIR}/assign_genes.py "
        "--output {output.gene_assigned} "
        "{input.bed} {input.chrom_gtf}"


rule add_gene_tags_to_bam:
    input:
        bam=CHROM_BAM_BC,
        bai=CHROM_BAM_BC_BAI,
        genes=CHROM_TSV_GENE_ASSIGNS,
    output:
        bam=CHROM_BAM_BC_GENE_TMP,
    conda:
        "../envs/umis.yml"
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/add_gene_tags.py "
        "--output {output.bam} "
        "{input.bam} {input.genes}"


rule cleanup_headers_3:
    input:
        CHROM_BAM_BC_GENE_TMP,
    output:
        bam=CHROM_BAM_BC_GENE,
        bai=CHROM_BAM_BC_GENE_BAI,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools reheader --no-PG -c 'grep -v ^@PG' "
        "{input} > {output.bam}; "
        "samtools index {output.bam}"


rule cluster_umis:
    input:
        bam=CHROM_BAM_BC_GENE,
        bai=CHROM_BAM_BC_GENE_BAI,
    output:
        bam=CHROM_BAM_FULLY_TAGGED_TMP,
    params:
        interval=config["UMI"]["GENOMIC_INTERVAL"],
    conda:
        "../envs/umis.yml"
    threads: 1
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/cluster_umis.py "
        "--threads {threads} "
        "--ref_interval {params.interval} "
        "--output {output.bam} {input.bam}"


rule cleanup_headers_4:
    input:
        CHROM_BAM_FULLY_TAGGED_TMP,
    output:
        bam=CHROM_BAM_FULLY_TAGGED,
        bai=CHROM_BAM_FULLY_TAGGED_BAI,
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools reheader --no-PG -c 'grep -v ^@PG' "
        "{input} > {output.bam}; "
        "samtools index {output.bam}"


def gather_chrom_bams(wildcards):
    # throw an Exception if checkpoint is pending
    checkpoint_dir = checkpoints.split_bam_by_chroms.get(**wildcards).output[0]
    return expand(
        CHROM_BAM_FULLY_TAGGED,  # <-- SHOULD BE FINAL SPLIT BAM FILE
        run_id=wildcards.run_id,
        chrom=glob_wildcards(os.path.join(checkpoint_dir, "{chrom}.sorted.bam")).chrom,
    )


rule combine_chrom_bams:
    input:
        chroms=gather_chrom_bams,
    output:
        all=BAM_FULLY_TAGGED,
        bai=BAM_FULLY_TAGGED_BAI,
    params:
        split_dir=lambda w: str(SPLIT_DIR).format(run_id=w.run_id),
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools merge -o {output.all} {input}; "
        "samtools index {output.all}; "
        # "rm -rf {params.split_dir}"


# def gather_chrom_barcode_counts(wildcards):
#     # throw an Exception if checkpoint is pending
#     checkpoint_dir = checkpoints.split_bam_by_chroms.get(**wildcards).output[0]
#     return expand(
#         CHROM_ASSIGNED_BARCODE_COUNTS, # <-- SHOULD BE FINAL SPLIT BAM FILE
#         run_id=wildcards.run_id,
#         chrom=glob_wildcards(os.path.join(checkpoint_dir, "{chrom}.sorted.bam")).chrom,
#     )
#
#
# rule combine_chrom_assigned_barcode_counts:
#     input:
#         chroms=gather_chrom_barcode_counts,
#     output:
#         all=ASSIGNED_BARCODE_COUNTS,
#     run:
#         import pandas as pd
#
#         dfs = [pd.read_csv(fn, sep="\t", header=False, names=["barcode", "count"]) \
#             for fn in input.chroms
#         ]
#         df = pd.concat(dfs, axis=0)
#         df = df.groupby("barcode")["count"].sum()
#         df.to_csv(output.all, sep="\t", index=False)


rule umi_gene_saturation:
    input:
        bam=BAM_FULLY_TAGGED,
        bai=BAM_FULLY_TAGGED_BAI,
    output:
        plot=SAT_PLOT,
    conda:
        "../envs/plotting.yml"
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/saturation.py "
        "--output {output.plot} {input.bam}"
