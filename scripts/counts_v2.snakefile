'''
Snakefile for running kallisto on RNAseq reads of pure Magnaporthe samples.

    Pipeline:
        1. count reads per Mo transcript per sample (kallisto)
        2. process count files into count matrix
        3. summarise number of reads processed

Run from do_count_matrix_v2.sh
'''

samples = []
fq1 = []
fq2 = []
with open("lib/samples_to_files.csv") as csv:
    for l in csv:
        l = l.rstrip("\n")
        if not l.startswith("name"):
            els = l.split(",")
            samples.append(els[3])
            fq1.append( els[6] )
            fq2.append( els[7] )

def sample_to_read(sample, samples, reads):
    return reads[samples.index(sample)]

configfile: "config.yaml"

# Removed percentage reads for taxid and read counts as no longer needed for pure Mo samples
rule counts:
    input: 
        counts="results/all_samples_tpm_matrix.txt",
        rawreadcount="results/read_counts.txt",
        metadata="results/run_metadata.txt"

# Following rules removed:
# count_reads_in_fq - ADDED AGAIN
# extract_percent_per_taxid
# count_raw_reads_in_fq - ADDED AGAIN
# run_kraken
# extract_reads
# count_extracted_reads_in_fq


rule count_raw_reads_in_fq:
    input:
        fq1=lambda wildcards: sample_to_read(wildcards.sample, samples, fq1)
    output:
        count=config['scratch'] + "{sample}/read_count.txt"
    threads: 1
    params:
        queue="tsl-short",
        mem="16G"
    shell: "zcat {input.fq1} | wc -l > {output.count}"

rule count_reads_in_fq:
    input:
        raw=expand(config['scratch'] + "{sample}/read_count.txt", sample=samples)
    threads: 1
    params:
        mem="16G",
        queue="tsl-short"
    output: "results/read_counts.txt"
    shell: "bash scripts/readcountsummary.sh {input.raw} > {output}"

rule kallisto_quant:
    input: 
        idx=config['scratch'] + "kallisto_indices/magnaporthe.idx",
        fq1=lambda wildcards: sample_to_read(wildcards.sample, samples, fq1),
        fq2=lambda wildcards: sample_to_read(wildcards.sample, samples, fq2)
    output: config['scratch'] + "{sample}/kallisto/abundance.tsv",
    params:
        folder=config['scratch'] + "{sample}/kallisto",
        mem="16G",
        queue="tsl-long",
        bootstraps=100
    threads: 1
    shell: "bash scripts/kallisto.sh {input.idx} {threads} {params.bootstraps} {input.fq1} {input.fq2} {params.folder}"

rule kallisto_index:
    input: config['reference_genome']
    output: config['scratch'] + "kallisto_indices/magnaporthe.idx"
    params:
        mem="32G",
        queue="tsl-short"
    shell: "bash scripts/idx.sh {input} {output}"

rule combine_counts:
    input: expand(config['scratch'] + "{sample}/kallisto/abundance.tsv", sample=samples)
    output: "results/all_samples_tpm_matrix.txt"
    params:
        mem="16G",
        queue="tsl-short"
    threads: 1
    shell: "python scripts/combine_tpm.py {input} > {output} "

rule sleuth_files:
    input: expand(config['scratch'] + "{sample}/kallisto/abundance.tsv", sample=samples)
    output:
        gz="results/kallisto_abundances.gz",
        meta="results/run_metadata.txt"
    params:
        mem="16G",
        queue="tsl-short",
        temp_dir=config['scratch'] + "/tmp/"
    threads: 1
    shell: "python scripts/make_metadata.py {params.temp_dir} {output.gz} {input} > {output.meta}"