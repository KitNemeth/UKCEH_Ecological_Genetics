# ==========================================================
# WGS Data Preparation and Adapter Trimming Pipeline
# ==========================================================

import os
import glob

########
# CONFIG
configfile: "YAML/config.yaml"

########
# PATHS
data = config["data"]
clean = config["clean"]
samples_file = config["samples_file"]

# Ensure clean directory exists
os.makedirs(clean, exist_ok=True)

# Read samplesfile.txt (provided manually)
samples = {}
with open(samples_file) as f:
    for line in f:
        seq_id, short = line.strip().split()
        samples[short] = seq_id

# List of sample names for expand()
SAMPLES = sorted(samples.keys())

########
# RULES

# Rule: create per-sample directories
rule make_sample_dirs:
    output:
        directory(os.path.join(clean, "{sample}"))

# AdapterRemoval rule
rule adapterremoval:
    input:
        r1=lambda wc: glob.glob(os.path.join(data, f"{samples[wc.sample]}*_R1*.fastq*"))[0],
        r2=lambda wc: glob.glob(os.path.join(data, f"{samples[wc.sample]}*_R2*.fastq*"))[0],
        dirs=rules.make_sample_dirs.output
    output:
        r1=os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"),
        r2=os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"),
        singleton=os.path.join(clean, "{sample}", "{sample}_singleton_truncated.fastq.gz"),
        collapsed=os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"),
        collapsed_trunc=os.path.join(clean, "{sample}", "{sample}_collapsed_truncated.fastq.gz"),
        discarded=os.path.join(clean, "{sample}", "{sample}_discarded.fastq.gz"),
        settings=os.path.join(clean, "{sample}", "{sample}.settings")
    shell:
        """
        module load tools
        module load gcc/15.1.0
        module load adapterremoval/2.3.1

        AdapterRemoval \
          --file1 {input.r1} \
          --file2 {input.r2} \
          --output1 {output.r1} \
          --output2 {output.r2} \
          --singleton {output.singleton} \
          --outputcollapsed {output.collapsed} \
          --outputcollapsedtruncated {output.collapsed_trunc} \
          --discarded {output.discarded} \
          --settings {output.settings} \
          --trimns \
          --trimqualities \
          --minquality 20 \
          --minlength 25 \
          --maxns 20 \
          --collapse \
          --gzip
        """

# Rule: all
rule all:
    input:
        expand(os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_singleton_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_collapsed_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_discarded.fastq.gz"), sample=SAMPLES)


