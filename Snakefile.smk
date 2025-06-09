# Snakemake workflow for Petrobium genome assembly
"""
Author: Kit Nemeth & Susheel Bhanu BUSI
Affiliation: UKCEH
Date: [2025-06-04]
Run: snakemake -s Snakefile.smk --use-conda --cores 4 -rp
Latest modification:
Purpose: Adding in steps for genome assembly using Petrobium long-read (nanopore) data.
"""

########
# MODULES
import os, re
import sys
import glob
import pandas as pd

########
# PATHS
FASTQ_DIR = "/hdd0/susbus/petrobium"
DBS_DIR = "/hdd0/susbus/databases"   # BUSCO databases
ENV_DIR = "/home/krinem/UKCEH_Ecological_Genetics/YAML"  # Conda environments directory
RESULTS_DIR = "/ssd0/krinem/petrobium_results"     # Results directory
SRC_DIR = "/home/krinem/UKCEH_Ecological_Genetics/Scripts"  # Scripts directory

########
# INPUT
SAMPLES = ["petrobium"]  # List of sample names

rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "preprocessed/{sample}_filtered.fq"),sample=SAMPLES)

########
# RULES
rule concatenate:
    input:
        fastq=FASTQ_DIR
    output:
        concatenated_fq=os.path.join(RESULTS_DIR, "concatenated/{sample}_concatenated.fastq.gz")
    message:
        "Concatenating all fastq files"
    shell:
        "cd {input.fastq} && cat {input.fastq}/*.fastq.gz > {output.concatenated_fq}"

rule quality_filtering:
    input:
        rules.concatenate.output.concatenated_fq
    output:
        filtered_fq=os.path.join(RESULTS_DIR, "preprocessed/{sample}_filtered.fq")
    message:
        "Filtering raw reads based on Quality score 10"
    conda:
        os.path.join(ENV_DIR, "nanofilt.yaml")
    params:
        quality=10, 
        length=500
    shell:
        "NanoFilt -q {params.quality} -l {params.length} {input} > {output.filtered_fq}"

rule assembler_flye:
    input:
        rules.quality_filtering.output.filtered_fq
    output:
        assembly=os.path.join(RESULTS_DIR, "assembly/flye/{sample}_assembly.fasta")
    message:
        "Assembling genome using Flye"
    conda:
        os.path.join(ENV_DIR, "flye.yaml")
    threads:
        48
    #shell:
        "flye --nano-raw {input} --out-dir $dirname({output.assembly}) --threads {threads}"

rule polish_medaka:
    input:
       asm=rules.assembler_flye.output.assembly,
       filt=rules.quality_filtering.output.filtered_fq
    output:
        polished_assembly=os.path.join(RESULTS_DIR, "polished/{sample}_polished.fasta")
    message:
        "Polishing assembly using Medaka"
    conda:
        os.path.join(ENV_DIR, "medaka.yaml")
    threads:
        48
    shell:
        "medaka_consensus -i {input.filt} -d {input.asm} -o $dirname({output.polished_assembly}) -t {threads}"
    
rule quality_assessment:
    input:
        polished_assembly=rules.polish_medaka.output.polished_assembly
    output:
        quality_ass=os.path.join(RESULTS_DIR, "quality/{sample}_quality_ass.txt")
    message:
        "Assessing quality of the polished assembly"
    conda:
        os.path.join(ENV_DIR, "busco.yaml")
    threads:
        48
    params:
        mode="genome"
    shell:
        "busco -i {input.polished_assembly} -m {param.mode} -c {threads} -o $dirname({output.quality_ass})"
