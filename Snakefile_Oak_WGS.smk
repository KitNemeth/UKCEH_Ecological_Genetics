
# ==========================================================
# WGS Data Preparation and Adapter Trimming Pipeline
# ==========================================================

import os

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
        expand(os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30_rmDup.bam"), sample=SAMPLES)

########
# AdapterRemoval
rule adapterremoval:
    input:
        r1=os.path.join(data, "{sample}_R1_001.fastq.gz"),
        r2=os.path.join(data, "{sample}_R2_001.fastq.gz")
    output:
        r1=os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"),
        r2=os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"),
        collapsed=os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"),
        settings=os.path.join(clean, "{sample}", "{sample}.settings")
    threads: 16
    resources:
        mem_mb=16000,
        io_slot=1   # limit number of concurrent AdapterRemoval jobs
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
          --collapse
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
    threads: 4
    benchmark:
        os.path.join(benchmark_dir, "bowtie2_index", "index.txt")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        bowtie2-build {input.fasta} {REF_INDEX}
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
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam")  # temp: deleted after q30 filter
    threads: 16
    resources:
        mem_mb=16000
    benchmark:
        os.path.join(benchmark_dir, "map_reads", "{sample}.txt")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        bowtie2 -p {threads} -x {REF_INDEX} --no-unal \
            -1 {input.r1} \
            -2 {input.r2} \
            -U {input.collapsed} | \
        samtools view -bS - | \
        samtools sort -o {output.bam}
        """

########
# Filter BAMs for Q30
rule filter_q30:
    input:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam")  # temp: deleted after dedup
    output:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam")  # temp: deleted after dedup
    benchmark:
        os.path.join(benchmark_dir, "filter_q30", "{sample}.txt")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view -b -q 30 {input.bam} > {output.bam}
        """

########
# Remove duplicates using Picard
rule remove_duplicates:
    input:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam")
    output:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30_rmDup.bam"),  # final BAM: kept permanently
        metrics=os.path.join(map_dir, "{sample}", "{sample}_Qrob_Dup.txt")      # kept permanently
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
            -M {output.metrics}
        """
