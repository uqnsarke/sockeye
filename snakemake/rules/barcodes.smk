rule extract_barcodes:
    input:
        STRANDED_FQ,
    output:
        fastq=BARCODES_UNCORR_READS,
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
        "--output_reads {output.fastq} "
        "--output_hq_bc {output.hq} "
        "{input}"


rule cluster_hq_barcodes:
    input:
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
        "{input.hq}"


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
        fastq=BARCODES_UNCORR_READS,
        whitelist=BARCODES_WHITELIST,
    output:
        all=BARCODES_CORR_READS,
        filtered=BARCODES_CORR_READS_FILTERED,
    params:
        batch_size=config["BARCODE"]["BATCH_SIZE"],
        max_ed=config["BARCODE"]["MAX_ED"],
        min_ed_diff=config["BARCODE"]["MIN_ED_DIFF"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/assign_barcodes.py "
        "-t {threads} "
        "--output_all {output.all} "
        "--output_filtered {output.filtered} "
        "--batch_size {params.batch_size} "
        "--max_ed {params.max_ed} "
        "--min_ed_diff {params.min_ed_diff} "
        "{input.fastq} {input.whitelist}"
