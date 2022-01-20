rule extract_barcodes:
    input:
        bam=BAM_SORT,
        barcodes=config["BC_SUPERLIST"],
    output:
        bam=BAM_BC_UNCORR,
        tsv=BARCODE_COUNTS,
    params:
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
        "{input.bam} {input.barcodes}; "
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


rule assign_barcodes:
    input:
        bam=BAM_BC_UNCORR,
        whitelist=BARCODE_WHITELIST,
    output:
        BAM_BC_CORR_UMI_UNCORR,
    params:
        max_ed=config["BARCODE"]["MAX_ED"],
        min_ed_diff=config["BARCODE"]["MIN_ED_DIFF"],
        read1=config["READ_STRUCTURE"]["READ1"],
        read1_suff_length=config["BARCODE"]["READ1_SUFF_LENGTH"],
        barcode_length=config["READ_STRUCTURE"]["BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE"]["UMI_LENGTH"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/assign_barcodes.py "
        "-t {threads} "
        "--output {output} "
        "--max_ed {params.max_ed} "
        "--min_ed_diff {params.min_ed_diff} "
        "--read1_adapter {params.read1} "
        "--read1_suff_length {params.read1_suff_length} "
        "--barcode_length {params.barcode_length} "
        "--umi_length {params.umi_length} "
        "{input.bam} {input.whitelist}; "
        "samtools index {output}"
