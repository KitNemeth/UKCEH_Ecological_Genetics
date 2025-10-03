# Snakemake workflow for genome assembly (species-agnostic, with logs)
"""
Author: Kit Nemeth & Susheel Bhanu BUSI
Affiliation: UKCEH
Date: [2025-06-04]
Run example:
    snakemake -s Snakefile.smk --use-conda --cores 48 -rp --config species=juniperus
Latest modification:
Purpose: Genome assembly using long-read (Nanopore) data, with species name set at runtime.
"""

########
# MODULES
import os

########
# CONFIG
species = config.get("species", "petrobium")   # default to petrobium if not provided

########
# PATHS
FASTQ_DIR = f"/ssd0/krinem/Sequence_data/{species}"
DBS_DIR     = "/hdd0/susbus/databases"   # BUSCO databases
ENV_DIR     = "/home/krinem/UKCEH_Ecological_Genetics/YAML"  # Conda environments directory
RESULTS_DIR = f"/ssd0/krinem/{species}_results"     # Results directory
SRC_DIR     = "/home/krinem/UKCEH_Ecological_Genetics/Scripts"  # Scripts directory

########
# INPUT
SAMPLES = [species]

rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "preprocessed/{sample}_filtered.fq"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "assembly/flye/{sample}_assembly.fasta"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "polished/{sample}_polished.fasta"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "quality/{sample}_quality_ass.txt"), sample=SAMPLES)

########
# RULES
import glob, os

rule concatenate:
    input:
        fastq=lambda wildcards: glob.glob(os.path.join(FASTQ_DIR, "*.fastq.gz"))
    output:
        concatenated_fq=os.path.join(RESULTS_DIR, "concatenated/{sample}_concatenated.fastq.gz"),
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}_concatenate.log")
    message:
        "Concatenating all fastq files"
    shell:
        """
        mkdir -p $(dirname {output.concatenated_fq}) $(dirname {log})
        cat {input.fastq} > "{output.concatenated_fq}" 2> "{log}"
        """

rule quality_filtering:
    input:
        rules.concatenate.output.concatenated_fq
    output:
        filtered_fq=os.path.join(RESULTS_DIR, "preprocessed/{sample}_filtered.fq")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}_quality_filtering.log")
    message:
        "Filtering raw reads based on Quality score 10"
    conda:
        "/ssd0/krinem/miniforge3/envs/nanofilt_env"
    params:
        quality=10, 
        length=500
    shell:
        """
        mkdir -p $(dirname {output.filtered_fq}) $(dirname {log})
        gunzip -c {input} | NanoFilt -q {params.quality} -l {params.length} > {output.filtered_fq} 2> {log}
        """

rule assembler_flye:
    input:
        rules.quality_filtering.output.filtered_fq
    output:
        assembly=os.path.join(RESULTS_DIR, "assembly/flye/{sample}_assembly.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}_flye.log")
    message:
        "Assembling genome using Flye"
    conda:
        os.path.join(ENV_DIR, "flye.yaml")
    threads:
        48
    shell:
        """
        outdir={output.assembly}.dir
        mkdir -p $outdir $(dirname {log})
        flye --nano-raw {input} --out-dir $outdir --threads {threads} > {log} 2>&1
        mv $outdir/assembly.fasta {output.assembly}
        """

rule polish_medaka:
    input:
        asm=rules.assembler_flye.output.assembly,
        filt=rules.quality_filtering.output.filtered_fq
    output:
        polished_assembly=os.path.join(RESULTS_DIR, "polished/{sample}_polished.fasta")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}_medaka.log")
    message:
        "Polishing assembly using Medaka"
    conda:
        os.path.join(ENV_DIR, "medaka.yaml")
    threads:
        48
    shell:
        """
        outdir=$(dirname {output.polished_assembly})/{wildcards.sample}_medaka_out

        mkdir -p $outdir $(dirname {log})

        medaka_consensus \
            -i {input.filt} \
            -d {input.asm} \
            -o $outdir \
            -t {threads} \
            > {log} 2>&1

        cp $outdir/consensus.fasta {output.polished_assembly}
        """

rule quality_assessment:
    input:
        polished_assembly=rules.polish_medaka.output.polished_assembly
    output:
        quality_ass=os.path.join(RESULTS_DIR, "quality/{sample}_quality_ass.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}_busco.log")
    message:
        "Assessing quality of the polished assembly"
    conda:
        os.path.join(ENV_DIR, "busco.yaml")
    threads:
        48
    params:
        mode="genome"
    shell:
        """
        mkdir -p $(dirname {output.quality_ass}) $(dirname {log})
        busco -i {input.polished_assembly} -m {params.mode} -c {threads} -o $dirname({output.quality_ass}) > {log} 2>&1
        """
