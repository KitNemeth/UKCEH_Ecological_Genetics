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

#rule assembly:

#rule polish:

#rule quality_assessment:
