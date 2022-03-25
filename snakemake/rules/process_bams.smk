rule extract_barcodes:
    input:
        bam=BAM_SORT,
    output:
        bam=temp(BAM_BC_UNCORR_TMP),
        tsv=BARCODE_COUNTS,
    params:
        barcodes=config["BC_SUPERLIST"],
        read1=config["READ_STRUCTURE_READ1"],
        read1_suff_length=config["BARCODE_READ1_SUFF_LENGTH"],
        barcode_length=config["READ_STRUCTURE_BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE_UMI_LENGTH"],
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
        bam=temp(BAM_BC_UNCORR),
        bai=temp(BAM_BC_UNCORR_BAI),
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
        flags=config["BARCODE_KNEEPLOT_FLAGS"],
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
        bam=temp(CHROM_BAM_BC_TMP),
        counts=CHROM_ASSIGNED_BARCODE_COUNTS,
    params:
        max_ed=config["BARCODE_MAX_ED"],
        min_ed_diff=config["BARCODE_MIN_ED_DIFF"],
        read1=config["READ_STRUCTURE_READ1"],
        read1_suff_length=config["BARCODE_READ1_SUFF_LENGTH"],
        barcode_length=config["READ_STRUCTURE_BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE_UMI_LENGTH"],
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
        bam=temp(CHROM_BAM_BC),
        bai=temp(CHROM_BAM_BC_BAI),
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools reheader --no-PG -c 'grep -v ^@PG' "
        "{input} > {output.bam}; "
        "samtools index {output.bam}"


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
        minQV=config["GENE_ASSIGNS_MINQV"],
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
        bam=temp(CHROM_BAM_BC_GENE_TMP),
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
        bam=temp(CHROM_BAM_BC_GENE),
        bai=temp(CHROM_BAM_BC_GENE_BAI),
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
        bam=temp(CHROM_BAM_FULLY_TAGGED_TMP),
    params:
        interval=config["UMI_GENOMIC_INTERVAL"],
        cell_gene_max_reads=config["UMI_CELL_GENE_MAX_READS"],
    conda:
        "../envs/umis.yml"
    threads: config["UMI_CLUSTER_MAX_THREADS"]
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/cluster_umis.py "
        "--threads {threads} "
        "--ref_interval {params.interval} "
        "--cell_gene_max_reads {params.cell_gene_max_reads} "
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
        chrom=glob_wildcards(
            # os.path.join(checkpoint_dir, "{chrom, r'\w+(\.\d+)?'}.sorted.bam")
            os.path.join(checkpoint_dir, "{chrom}.sorted.bam")
        ).chrom,
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
        "rm -rf {params.split_dir}"


rule construct_expression_matrix:
    input:
        bam=BAM_FULLY_TAGGED,
        bai=BAM_FULLY_TAGGED_BAI,
    output:
        tsv=MATRIX_COUNTS_TSV,
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/gene_expression.py "
        "--output {output.tsv} {input.bam}"


rule process_expression_matrix:
    input:
        tsv=MATRIX_COUNTS_TSV,
    output:
        tsv=MATRIX_PROCESSED_TSV,
    params:
        min_genes=config["MATRIX_MIN_GENES"],
        min_cells=config["MATRIX_MIN_CELLS"],
        max_mito=config["MATRIX_MAX_MITO"],
        norm_count=config["MATRIX_NORM_COUNT"],
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/process_matrix.py "
        "--min_genes {params.min_genes} "
        "--min_cells {params.min_cells} "
        "--max_mito {params.max_mito} "
        "--norm_count {params.norm_count} "
        "--output {output.tsv} {input.tsv}"


rule umap_reduce_expression_matrix:
    input:
        tsv=MATRIX_PROCESSED_TSV,
    output:
        tsv=MATRIX_UMAP_TSV,
    conda:
        "../envs/umap.yml"
    shell:
        "python {SCRIPT_DIR}/umap_reduce.py "
        "--output {output.tsv} {input.tsv}"


rule umap_plot_total_umis:
    input:
        umap=MATRIX_UMAP_TSV,
        matrix=MATRIX_PROCESSED_TSV,
    output:
        plot=MATRIX_UMAP_PLOT_TOTAL,
    conda:
        "../envs/plotting.yml"
    shell:
        "python {SCRIPT_DIR}/plot_umap.py "
        "--output {output.plot} "
        "{input.umap} {input.matrix}"


rule umap_plot_genes:
    input:
        umap=MATRIX_UMAP_TSV,
        matrix=MATRIX_PROCESSED_TSV,
    output:
        plot=MATRIX_UMAP_PLOT_GENE,
    params:
        gene=lambda wc: wc.plot_gene,
    conda:
        "../envs/plotting.yml"
    shell:
        "python {SCRIPT_DIR}/plot_umap.py "
        "--gene {params.gene} "
        "--output {output.plot} "
        "{input.umap} {input.matrix}"


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
