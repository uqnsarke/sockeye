gzip -cd refdata-gex-GRCh38-2020-A/fasta/genome.fa.gz > refdata-gex-GRCh38-2020-A/fasta/genome.fa
gzip -cd refdata-gex-GRCh38-2020-A/genes/genes.gtf.gz > refdata-gex-GRCh38-2020-A/genes/genes.gtf
rm -rf ./output
snakemake --snakefile ../Snakefile --use-conda --configfile config.test.yml -j 2 -pr all
