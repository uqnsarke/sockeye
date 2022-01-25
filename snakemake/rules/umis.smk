rule cluster_umis:
    input:
        bam=BAM_BC_CORR_UMI_UNCORR_GENE,
        bai=BAM_BC_CORR_UMI_UNCORR_GENE_BAI,
    output:
        bam=BAM_BC_CORR_UMI_CORR_GENE,
        bai=BAM_BC_CORR_UMI_CORR_GENE_BAI,
    conda:
        "../envs/umis.yml"
    threads: config["MAX_THREADS"]
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/cluster_umis.py "
        "--threads {threads} "
        "--output {output.bam} {input.bam}"


rule umi_gene_saturation:
    input:
        bam=BAM_BC_CORR_UMI_CORR_GENE,
        bai=BAM_BC_CORR_UMI_CORR_GENE_BAI,
    output:
        plot=SAT_PLOT,
    conda:
        "../envs/plotting.yml"
    shell:
        "touch {input.bai}; "
        "python {SCRIPT_DIR}/saturation.py "
        "--output {output.plot} {input.bam}"
