# ***********
# * FIREFLY *
# ***********


# **** Imports ****

import glob
from snakemake.utils import R, report


# **** Variables ****

configfile: "config.yaml"

VERSION='0.2.0'


samples = config["samples"]

# raw file list
glob_pat_r1 = expand("{dir}/{sample}_*R1*.fastq*", dir=config["read_directory"], sample=samples)
glob_pat_r2 = expand("{dir}/{sample}_*R2*.fastq*", dir=config["read_directory"], sample=samples)
raw_r1 = [glob.glob(x) for x in glob_pat_r1]
raw_r2 = [glob.glob(x) for x in glob_pat_r2]


# **** Rules ****

rule all:
    input: 
        "results/taxonomy_out.rda", 
        "multiqc_data/multiqc_fastqc.txt", 
        "results/otus.tre"

rule qc:
    input: 
        "multiqc_data/multiqc_fastqc.txt",
        expand("clipped/{sample}_R1.cut", sample = samples)

rule clip_primers:
    input: 
        r1 = lambda wc: glob.glob("{dir}/{sample}_*R1*.fastq*".format(dir=config["read_directory"], sample=wc.sample)),
        r2 = lambda wc: glob.glob("{dir}/{sample}_*R2*.fastq*".format(dir=config["read_directory"], sample=wc.sample))
    output: 
        r1 = "clipped/{sample}_R1.cut", 
        r2 = "clipped/{sample}_R2.cut"
    log: 
        "logs/{sample}.log"
    conda:
        "envs/quality.yaml"
    shell: 
        """ 
	    cutadapt -g {config[fwd_primer]} -G {config[rev_primer]} \
            -a {config[rev_primer_rc]} -A {config[fwd_primer_rc]} \
            -m 50 -q {config[q_trim]} \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} >> {log}
        """

rule fastqc:
    input: 
        raw_r1 + raw_r2 
    output: 
        touch("fastqc.done")
    threads: 
        config["num_threads"] 
    conda:
        "envs/quality.yaml"
    shell: 
        "mkdir -p fastqc; fastqc -t {threads} -o fastqc {input}"

rule multiqc:
    input: 
        "fastqc.done"
    output: 
        "multiqc_data/multiqc_fastqc.txt"
    conda:
        "envs/quality.yaml"
    shell: 
        "multiqc -f fastqc"

rule filter_and_trim:
    input: 
        r1 = expand("clipped/{sample}_R1.cut", sample = samples),
        r2 = expand("clipped/{sample}_R2.cut", sample = samples)
    output: 
        r1 = expand("filtered/{sample}_R1.fq", sample = samples),
        r2 = expand("filtered/{sample}_R2.fq", sample = samples),
        filt_out = "results/filtered.rds"
    threads: 
        config["num_threads"]
    conda:
        "envs/dada2.yaml"
    script: 
        "scripts/filter_and_trim.R"

rule learn_errors:
    input: 
        r1 = rules.filter_and_trim.output.r1,
        r2 = rules.filter_and_trim.output.r2
    output: 
        "results/error_rates.rda" 
    threads: 
        config["num_threads"]
    conda:
        "envs/dada2.yaml"
    script: 
        "scripts/learn_errors.R"

rule infer_seqs:
    input: 
        rules.learn_errors.output
    output: 
        "results/merged.rds"
    threads: 
        config["num_threads"]
    conda:
        "envs/dada2.yaml"
    script: 
        "scripts/infer_seqs.R"

rule taxonomy:
    input: 
        rules.infer_seqs.output
    output: 
        taxonomy = "results/taxonomy_out.rda", 
        otus = "results/otus.fasta"
    params: 
        train_dir = "training"
    threads: 
        config["num_threads"]
    conda:
        "envs/dada2.yaml"
    script: 
        "scripts/assign_taxonomy.R"

rule ssu_align:
    input:
        rules.taxonomy.output.otus
    output: 
        "ssu_out/ssu_out.bacteria.stk"
    params: 
        dir = "ssu_out"
    log: 
        "logs/align.log"
    conda:
        "envs/align.yaml"
    shell: 
        "ssu-align -f {input} {params.dir} &>> {log}"

rule ssu_mask:
    input: 
        rules.ssu_align.output
    output: 
        "ssu_out/ssu_out.bacteria.mask.afa"
    params: 
        dir = "ssu_out"
    log: 
        "logs/align.log"
    conda:
        "envs/align.yaml"
    shell: 
        "ssu-mask --afa {params.dir} &>> {log}"

rule tree:
    input: 
        rules.ssu_mask.output
    output: 
        "results/otus.tre"
    log: 
        "logs/tree.log"
    conda:
        "envs/align.yaml"
    shell: 
        "FastTree -nt {input} >{output} 2> {log}"

