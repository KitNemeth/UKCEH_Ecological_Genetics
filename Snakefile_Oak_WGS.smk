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
REF = config["REF"]

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
# Batching configuration for Bowtie2
BATCH_SIZE = 15
BATCHES = [SAMPLES[i:i+BATCH_SIZE] for i in range(0, len(SAMPLES), BATCH_SIZE)]
BATCH_IDS = list(range(len(BATCHES)))

########
# Rule: all
rule all:
    input:
        # AdapterRemoval outputs
        expand(os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_singleton_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_collapsed_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_discarded.fastq.gz"), sample=SAMPLES),

        # Bowtie2 mapped BAMs
        expand(os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam"), sample=SAMPLES),

        # Filtered BAMs (q30)
        expand(os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam"), sample=SAMPLES),

        # Deduplicated BAMs
        expand(os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30_rmDup.bam"), sample=SAMPLES)

########
# RULES

# AdapterRemoval
rule adapterremoval:
    input:
        r1=lambda wc: os.path.join(data, f"{wc.sample}_R1_001.fastq.gz"),
        r2=lambda wc: os.path.join(data, f"{wc.sample}_R2_001.fastq.gz")
    output:
        r1=os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"),
        r2=os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"),
        singleton=os.path.join(clean, "{sample}", "{sample}_singleton_truncated.fastq.gz"),
        collapsed=os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"),
        collapsed_trunc=os.path.join(clean, "{sample}", "{sample}_collapsed_truncated.fastq.gz"),
        discarded=os.path.join(clean, "{sample}", "{sample}_discarded.fastq.gz"),
        settings=os.path.join(clean, "{sample}", "{sample}.settings")
    conda:
        os.path.join(ENV_DIR, "adapterremoval.yaml")
    shell:
        """
        mkdir -p $(dirname {output.r1})
        AdapterRemoval \
          --file1 "{input.r1}" \
          --file2 "{input.r2}" \
          --settings "{output.settings}" \
          --output1 "{output.r1}" \
          --output2 "{output.r2}" \
          --singleton "{output.singleton}" \
          --outputcollapsed "{output.collapsed}" \
          --outputcollapsedtruncated "{output.collapsed_trunc}" \
          --discarded "{output.discarded}" \
          --trimns \
          --trimqualities \
          --minquality 20 \
          --minlength 25 \
          --maxns 20 \
          --collapse
        """

# Bowtie2 mapping in batches
rule map_reads_batch:
    input:
        r1=lambda wc: [os.path.join(clean, s, f"{s}_R1_truncated.fastq.gz") for s in BATCHES[wc.batch]],
        r2=lambda wc: [os.path.join(clean, s, f"{s}_R2_truncated.fastq.gz") for s in BATCHES[wc.batch]],
        collapsed=lambda wc: [os.path.join(clean, s, f"{s}_collapsed.fastq.gz") for s in BATCHES[wc.batch]]
    output:
        bam=lambda wc: [os.path.join(map_dir, s, f"{s}_Qrob.bam") for s in BATCHES[wc.batch]]
    params:
        batch=lambda wc: BATCHES[wc.batch]
    conda: os.path.join(ENV_DIR, "bowtie2.yaml")
    threads: 4  # multi-threaded per sample
    shell:
        """
        for s in {params.batch}; do
            mkdir -p {map_dir}/$s
            bowtie2 -p {threads} -x {REF} --no-unal \
                -1 {clean}/$s/${s}_R1_truncated.fastq.gz \
                -2 {clean}/$s/${s}_R2_truncated.fastq.gz \
                -U {clean}/$s/${s}_collapsed.fastq.gz | \
            samtools view -bS - | \
            samtools sort -o {map_dir}/$s/${s}_Qrob.bam
        done
        """

# Filter BAMs for Q30
rule filter_q30:
    input:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob.bam")
    output:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam")
    conda:
        os.path.join(ENV_DIR, "bowtie2.yaml")
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view -b -q 30 {input.bam} > {output.bam}
        """

# Remove duplicates using Picard
rule remove_duplicates:
    input:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30.bam")
    output:
        bam=os.path.join(map_dir, "{sample}", "{sample}_Qrob_q30_rmDup.bam"),
        metrics=os.path.join(map_dir, "{sample}", "{sample}_Qrob_Dup.txt")
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
