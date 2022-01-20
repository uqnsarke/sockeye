rule call_paftools:
    input:
        gtf=str(REF_GENES_GTF),
    output:
        bed=REF_GENES_BED,
    conda:
        "../envs/minimap2.yml"
    shell:
        "paftools.js gff2bed -j {input.gtf} > {output.bed}"


rule get_chrom_sizes:
    input:
        genome=config["REF_GENOME_FASTA"] + ".fai",
    output:
        chrsizes=REF_CHROM_SIZES,
    shell:
        "cut -f1,2 {input.genome} | sort -V > {output.chrsizes}"


rule align_to_ref:
    input:
        fastq=STRANDED_FQ,
        bed=REF_GENES_BED,
        chrom_sizes=REF_CHROM_SIZES,
    output:
        sam_tmp=temp(SAM_TMP),
        unsort_bam=temp(BAM_UNSORT_TMP),
        sort_bam=BAM_SORT,
    params:
        ref=config["REF_GENOME_FASTA"],
    threads: config["MAX_THREADS"]
    conda:
        "../envs/minimap2.yml"
    shell:
        "minimap2 -ax splice -uf --MD -t {threads} "
        "--junc-bed {input.bed} "
        "--secondary=no "
        "{params.ref} {input.fastq} > "
        "{output.sam_tmp} && "
        "samtools view {output.sam_tmp} "
        "-t {input.chrom_sizes} -o {output.unsort_bam}; "
        "samtools sort {output.unsort_bam} -o {output.sort_bam}; "
        "samtools index {output.sort_bam}"
