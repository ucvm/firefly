# *******************************
# * Vet med microbiome pipeline *
# *******************************

# **** Imports ****

import glob
from snakemake.utils import R, report


# **** Variables ****

configfile: "config.yaml"

VERSION='0.1.0'

# **** Conda ****

# set location of conda install - being explicit about this avoids problems with multiple conda installs
CONDA="/export/home/mworken/miniconda3/bin/"

# assumes env called 'firefly' already created, use the provided vmmp_env.yml file to create
ENV="source {CONDA}/activate firefly".format(CONDA = CONDA)

# prefix shell commands with the above source command
shell.prefix(ENV + '; ')

samples = config["samples"]

# raw file list
glob_pat_r1 = expand("{dir}/{sample}_*R1*.fastq*", dir=config["read_directory"], sample=samples)
glob_pat_r2 = expand("{dir}/{sample}_*R2*.fastq*", dir=config["read_directory"], sample=samples)
raw_r1 = [glob.glob(x) for x in glob_pat_r1]
raw_r2 = [glob.glob(x) for x in glob_pat_r2]


# **** Rules ****

rule all:
    input: "results/taxonomy_out.rda", "multiqc_data/multiqc_fastqc.txt", "results/otus.tre"

rule clip_primers:
    input: r1 = lambda wc: glob.glob("{dir}/{sample}_*R1*.fastq*".format(dir=config["read_directory"], sample=wc.sample)),
           r2 = lambda wc: glob.glob("{dir}/{sample}_*R2*.fastq*".format(dir=config["read_directory"], sample=wc.sample))
    output: r1="clipped/{sample}_R1.cut", r2="clipped/{sample}_R2.cut"
    log: "logs/{sample}.log"
    shell: """ 
	    cutadapt -e 0 -O 17 -g {config[fwd_primer]} -G {config[rev_primer]} \
            -a {config[rev_primer_rc]} -A {config[fwd_primer_rc]} \
            -m 50 -q {config[q_trim]} \
            -o {output.r1} -p {output.r2} {input.r1} {input.r2} >> {log}
            """
rule fastqc:
    input: raw_r1 + raw_r2 
    output: touch("fastqc.done")
    threads: config["num_threads"] 
    shell: "mkdir -p fastqc; fastqc -t {threads} -o fastqc {input}"

rule multiqc:
    input: "fastqc.done"
    output: "multiqc_data/multiqc_fastqc.txt"
    shell: "multiqc -f fastqc"

rule filter_and_trim:
    input: r1 = expand("clipped/{sample}_R1.cut", sample = samples),
           r2 = expand("clipped/{sample}_R2.cut", sample = samples)
    output: r1 = expand("filtered/{sample}_R1.fq", sample = samples),
            r2 = expand("filtered/{sample}_R2.fq", sample = samples),
            filt_out = "results/filtered.rds"
    threads: config["num_threads"]
    script: "filter_and_trim.R"

rule learn_errors:
    input: r1 = rules.filter_and_trim.output.r1,
           r2 = rules.filter_and_trim.output.r2
    output: "results/error_rates.rda" 
    threads: config["num_threads"]
    script: "learn_errors.R"

rule infer_seqs:
    input: rules.learn_errors.output
    output: "results/merged.rds"
    threads: config["num_threads"]
    script: "infer_seqs.R"

rule taxonomy:
    input: rules.infer_seqs.output
    output: taxonomy = "results/taxonomy_out.rda", otus = "results/otus.fasta"
    threads: config["num_threads"]
    script: "assign_taxonomy.R"

rule ssu_align:
    input: rules.taxonomy.output.otus
    output: "ssu_out/ssu_out.bacteria.stk"
    params: dir="ssu_out"
    log: "logs/align.log"
    shell: "~/bin/ssu-align -f {input} {params.dir} &>> {log}"

rule ssu_mask:
    input: rules.ssu_align.output
    output: "ssu_out/ssu_out.bacteria.mask.afa"
    params: dir="ssu_out"
    log: "logs/align.log"
    shell: "~/bin/ssu-mask --afa {params.dir} &>> {log}"

rule tree:
    input: rules.ssu_mask.output
    output: "results/otus.tre"
    log: "logs/tree.log"
    shell: "FastTree -nt {input} >{output} 2> {log}"

