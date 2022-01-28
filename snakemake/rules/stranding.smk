rule cp_batch_fastqs:
    input:
        dir=lambda wildcards: sample_sheet.loc[wildcards.run_id],
    output:
        fofn=FOFN,
    params:
        dir=directory(COPIED_DIR),
    shell:
        "mkdir -p {params.dir}; "
        "ls -d {input.dir}/* > {output.fofn}"


checkpoint call_cat_fastq:
    input:
        FOFN,
    output:
        dir=directory(CHUNKED_FASTQ_DIR),
    threads: config["MAX_THREADS"]
    shell:
        "python {SCRIPT_DIR}/cat_fastqs.py "
        "--threads {threads} "
        "--output_dir {output.dir} "
        "{input}"


rule call_adapter_scan:
    input:
        CHUNKED_FASTQS,
    output:
        tsv=READ_CONFIG_CHUNKED,
        fastq=STRANDED_FQ_CHUNKED,
    threads: 1
    params:
        batch_size=config["READ_STRUCTURE"]["BATCH_SIZE"],
    conda:
        "../envs/stranding.yml"
    shell:
        "python {SCRIPT_DIR}/adapter_scan_vsearch.py "
        "-t {threads} "
        "--output_fastq {output.fastq} "
        "--output_tsv {output.tsv} "
        "--batch_size {params.batch_size} "
        "{input} "


def gather_tsv_files_from_run(wildcards):
    # throw and Exception if checkpoint is pending
    checkpoint_dir = checkpoints.call_cat_fastq.get(**wildcards).output[0]
    return expand(
        READ_CONFIG_CHUNKED,
        run_id=wildcards.run_id,
        batch_id=glob_wildcards(
            os.path.join(checkpoint_dir, "{batch_id}.fastq.gz")
        ).batch_id,
    )


rule combine_adapter_tables:
    input:
        tsv_files=gather_tsv_files_from_run,
    output:
        READ_CONFIG,
    run:
        import os
        import pandas as pd

        dfs = [pd.read_csv(fn, sep="\t") for fn in input.tsv_files]
        df = pd.concat(dfs, axis=0)
        df.to_csv(output[0], sep="\t", index=False)
        [os.remove(fn) for fn in input]


def gather_fastq_files_from_run(wildcards):
    checkpoint_dir = checkpoints.call_cat_fastq.get(**wildcards).output[0]
    return expand(
        STRANDED_FQ_CHUNKED,
        run_id=wildcards.run_id,
        batch_id=glob_wildcards(
            os.path.join(checkpoint_dir, "{batch_id}.fastq.gz")
        ).batch_id,
    )


rule combine_stranded_fastqs:
    input:
        fastq_files=gather_fastq_files_from_run,
    output:
        STRANDED_FQ,
    run:
        import shutil

        with open(output[0], "wb") as f_out:
            for chunked_fastq in input.fastq_files:
                with open(chunked_fastq, "rb") as f_:
                    shutil.copyfileobj(f_, f_out)
        [os.remove(fn) for fn in input]


rule summarize_adapter_table:
    input:
        READ_CONFIG,
    output:
        CONFIG_STATS,
    params:
        run_id="{run_id}",
    resources:
        mem=30,
    run:
        import pandas as pd
        import json

        df = pd.read_csv(input[0], sep="\t")
        stats = {}
        stats[params.run_id] = {}
        stats[params.run_id]["general"] = {}
        stats[params.run_id]["general"]["n_reads"] = df.shape[0]
        stats[params.run_id]["general"]["rl_mean"] = df["readlen"].mean()
        stats[params.run_id]["general"]["rl_std_dev"] = df["readlen"].std()
        stats[params.run_id]["general"]["n_fl"] = df[df["fl"] == True].shape[0]
        stats[params.run_id]["general"]["n_stranded"] = df[
            df["stranded"] == True
        ].shape[0]

        stats[params.run_id]["strand_counts"] = {}
        stats[params.run_id]["strand_counts"]["n_plus"] = df[
            df["orig_strand"] == "+"
        ].shape[0]
        stats[params.run_id]["strand_counts"]["n_minus"] = df[
            df["orig_strand"] == "-"
        ].shape[0]

        stats[params.run_id]["detailed_config"] = {}
        for category, n in df["orig_adapter_config"].value_counts().items():
            stats[params.run_id]["detailed_config"][category] = n

        stats[params.run_id]["summary_config"] = {}
        for label, n in df["lab"].value_counts().items():
            stats[params.run_id]["summary_config"][label] = n

        with open(output[0], "w") as f:
            json.dump(stats, f, indent=4)
