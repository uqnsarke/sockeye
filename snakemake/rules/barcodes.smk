rule extract_barcodes:
    input:
        STRANDED_FQ,
    output:
        all=BARCODES,
        hq=BARCODES_HQ,
    params:
        read1=config["READ_STRUCTURE"]["READ1"],
        read1_suff_length=config["BARCODE"]["READ1_SUFF_LENGTH"],
        batch_size=config["BARCODE"]["BATCH_SIZE"],
        barcode_length=config["READ_STRUCTURE"]["BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE"]["UMI_LENGTH"],
        bc_superlist=config["BC_SUPERLIST"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/extract_barcode.py "
        "-t {threads} "
        "--read1_adapter {params.read1} "
        "--read1_suff_length {params.read1_suff_length} "
        "--batch_size {params.batch_size} "
        "--barcode_length {params.barcode_length} "
        "--umi_length {params.umi_length} "
        "--bc_superlist {params.bc_superlist} "
        "--output_bc {output.all} "
        "--output_hq_bc {output.hq} "
        "{input}"


rule cluster_hq_barcodes:
    input:
        all=BARCODES,
        hq=BARCODES_HQ,
    output:
        BARCODES_CLUSTERS,
    params:
        bc_len=config["READ_STRUCTURE"]["BARCODE_LENGTH"],
        min_id=config["BARCODE"]["MIN_ID"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/call_vsearch_cluster.py "
        "--output {output} "
        "--bc_len {params.bc_len} "
        "--min_id {params.min_id} "
        "-t {threads} "
        "{input.all} {input.hq}"


rule generate_whitelist:
    input:
        BARCODES_CLUSTERS,
    output:
        whitelist=BARCODES_WHITELIST,
        kneeplot=BARCODES_KNEEPLOT,
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


rule assign_barcodes:
    input:
        barcodes=BARCODES,
        whitelist=BARCODES_WHITELIST,
    output:
        BARCODES_ASSIGNS,
    params:
        batch_size=config["BARCODE"]["BATCH_SIZE"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/assign_barcodes.py "
        "-t {threads} "
        "--output {output} "
        "--batch_size {params.batch_size} "
        "{input.barcodes} {input.whitelist}"


rule filter_barcode_assigns:
    input:
        BARCODES_ASSIGNS,
    output:
        bc=BARCODES_ASSIGNS_FILTERED,
        counts=BARCODES_ASSIGNS_FILTERED_COUNTS,
    params:
        max_ed=config["BARCODE"]["MAX_ED"],
        min_ed_diff=config["BARCODE"]["MIN_ED_DIFF"],
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/filter_barcode_assigns.py "
        "--max_ed {params.max_ed} "
        "--min_ed_diff {params.min_ed_diff} "
        "--output_bc {output.bc} "
        "--output_counts {output.counts} "
        "{input}"
