rule extract_umis:
    input:
        fastq=STRANDED_FQ,
        barcodes=BARCODES_ASSIGNS_FILTERED,
    output:
        UMI_EXTRACTED_READS,
    params:
        read1=config["READ_STRUCTURE"]["READ1"],
        read1_suff_length=config["BARCODE"]["READ1_SUFF_LENGTH"],
        batch_size=config["BARCODE"]["BATCH_SIZE"],
        barcode_length=config["READ_STRUCTURE"]["BARCODE_LENGTH"],
        umi_length=config["READ_STRUCTURE"]["UMI_LENGTH"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/barcodes.yml"
    shell:
        "python {SCRIPT_DIR}/extract_umi.py "
        "-t {threads} "
        "--read1_adapter {params.read1} "
        "--read1_suff_length {params.read1_suff_length} "
        "--batch_size {params.batch_size} "
        "--barcode_length {params.barcode_length} "
        "--umi_length {params.umi_length} "
        "--output {output} "
        "{input.fastq} {input.barcodes}"
