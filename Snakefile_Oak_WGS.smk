# ==========================================================
# WGS Data Preparation and Adapter Trimming Pipeline
# ==========================================================

import os
#run with snakemake s- Snakefile_Oak_WGS.smk -cores -j ...
########
# CONFIG
configfile: "YAML/config.yaml"

########
# PATHS
data = config["data"]
ENV_DIR = "/ssd0/krinem/UKCEH_Ecological_Genetics/YAML"  # Conda environments directory
clean = config["clean"]
samples_file = config["samples_file"]
analyses = config["analyses"]
map_dir = config["map"]
REF_FASTA = config["reference"]["fasta"]
REF_INDEX = config["reference"]["index_prefix"]
benchmark_dir = "benchmarks"

########
# Read samplesfile.txt
samples = {}
with open(samples_file) as f:
    for line in f:
        seq_id, short = line.strip().split()
        samples[short] = seq_id

# List of sample names
SAMPLES = sorted(samples.keys())

########
# Rule: all
rule all:
    input:
        [
            # AdapterRemoval outputs
            #os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"),
            #os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"),
            #os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"),

            # Bowtie2 mapped BAMs
            #os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam"),

            # Filtered BAMs (q30)
            #os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam"),

            # Final outputs
            expand(
                os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30_rmDup.bam"),
                sample=SAMPLES
            )
        ]

########
# AdapterRemoval
rule adapterremoval:
    input:
        r1=os.path.join(data, "{sample}_R1_001.fastq.gz"),
        r2=os.path.join(data, "{sample}_R2_001.fastq.gz")
    output:
        # Only keep outputs needed for mapping
        r1=temp(os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz")),
        r2=temp(os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz")),
        collapsed=temp(os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz")),
        # Optional: keep .settings for logging / QC
        settings=os.path.join(clean, "{sample}", "{sample}.settings")
    log:
        os.path.join(clean, "{sample}", "adapterremoval.log")
    threads: 16
    benchmark:
        os.path.join(benchmark_dir, "adapterremoval", "{sample}.txt")
    conda:
        os.path.join(ENV_DIR, "adapterremoval.yaml")
    retries: 3      # <--- retry up to 3 times if the job fails
    shell:
         """
        mkdir -p $(dirname {output.r1}) $(dirname {output.r2}) $(dirname {output.collapsed})
        export TMPDIR=/hdd0/krinem/tmp
        AdapterRemoval \
          --file1 {input.r1} \
          --file2 {input.r2} \
          --settings {output.settings} \
          --output1 {output.r1} \
          --output2 {output.r2} \
          --outputcollapsed {output.collapsed} \
          --trimns \
          --trimqualities \
          --minquality 20 \
          --minlength 25 \
          --collapse \
          > {log} 2>&1
        """

########
# Build Bowtie2 index
rule bowtie2_index:
    input:
        fasta=REF_FASTA
    output:
        expand(
            REF_INDEX + ".{n}.bt2",
            n=[1, 2, 3, 4, "rev.1", "rev.2"]
        )
    log:
        os.path.join(benchmark_dir, "bowtie2_index", "index.log")
    threads: 4
    benchmark:
        os.path.join(benchmark_dir, "bowtie2_index", "index.txt")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        bowtie2-build {input.fasta} {REF_INDEX} > {log} 2>&1
        """

########
# Bowtie2 mapping (per sample, parallelizable with -j)
rule map_reads:
    input:
        r1=os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"),
        r2=os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"),
        collapsed=os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"),
        index=expand(
            REF_INDEX + ".{n}.bt2",
            n=[1, 2, 3, 4, "rev.1", "rev.2"]
        )
    output:
        bam=temp(os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam"))
    log:
        os.path.join(map_dir, "{sample}", "map_reads.log")
    threads: 16
    benchmark:
        os.path.join(benchmark_dir, "map_reads", "{sample}.txt")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam})

        # Redirect all output (stdout + stderr) to log
        {{
            echo "=== Starting Bowtie2 mapping for {wildcards.sample} ==="
            date

            echo "--- bowtie2 ---"
            TMPDIR=/tmp bowtie2 -p {threads} -x {REF_INDEX} --no-unal \
                -1 {input.r1} -2 {input.r2} -U {input.collapsed}

            echo "--- samtools view ---"
            samtools view -bS - 

            echo "--- samtools sort ---"
            samtools sort -m 4G -@ 8 -T /tmp/sort_{wildcards.sample} -o {output.bam}

            echo "=== Finished mapping for {wildcards.sample} ==="
            date
        }} &> {log}
        """

########
# Filter BAMs for Q30
rule filter_q30:
    input:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam")  # temp: deleted after dedup
    output:
        bam=temp(os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam"))  # temp: deleted after dedup
    log:
        os.path.join(map_dir, "{sample}", "filter_q30.log")
    benchmark:
        os.path.join(benchmark_dir, "filter_q30", "{sample}.txt")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view -b -q 30 {input.bam} > {output.bam} 2> {log}
        """

########
# Remove duplicates using Picard
rule remove_duplicates:
    input:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam")
    output:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30_rmDup.bam"),  # final BAM: kept permanently
        metrics=os.path.join(map_dir, "{sample}", "{sample}_Qrob_Dup.txt")      # kept permanently
    log:
        os.path.join(map_dir, "{sample}", "remove_duplicates.log")
    benchmark:
        os.path.join(benchmark_dir, "remove_duplicates", "{sample}.txt")
    conda:
        os.path.join(ENV_DIR, "picard.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        picard MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            > {log} 2>&1
        """
