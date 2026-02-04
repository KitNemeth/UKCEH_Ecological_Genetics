# ==========================================================
# WGS Data Preparation and Adapter Trimming Pipeline
# ==========================================================
"""
Author      : Kit Nemeth
Affiliation : UKCEH
Date        : 2026-02-04

Purpose:
    Preprocessing of short-read Whole Genome Sequencing (WGS) data:
    - Generate per-sample directories
    - Trim adapters and low-quality bases using AdapterRemoval
    - Produce clean paired-end and collapsed reads ready for downstream analysis

Usage Example:
    snakemake -s WGS_prep.smk --use-conda --cores 48 -rp

Notes:
    - Supports variable FASTQ naming conventions, including differing lanes (_Lxxx), sample indexes (_Sxx), and optional R1/R2 suffixes
    - Automatically generates a sample table (samplesfile.txt) if missing
    - Compatible with HPC module systems; loads AdapterRemoval via module load commands
    - All intermediate outputs are stored in the 'clean' directory
"""

########
# MODULES
import os
import glob

########
# CONFIG
configfile: "YAML/config.yaml"

########
# PATHS
SAMP = config["SAMP"]
data = config["data"]
clean = config["clean"]
analyses = config["analyses"]
samples_file = config["samples_file"]

# read sample table
samples = {}
with open(config["samples_file"]) as f:
    for line in f:
        seq_id, sample = line.strip().split()
        samples[sample] = seq_id

# Ensure the main clean directory exists
os.makedirs(clean, exist_ok=True)
os.makedirs(analyses, exist_ok=True)

# Collect FASTQ files (ignore hidden files like .DS_Store or ._*)
fastq_files = [f for f in glob.glob(os.path.join(data, "*.fastq")) if not os.path.basename(f).startswith("._")]

# Build sample dictionary: short_name -> full seq_id
samples = {}

# Generate samples_file automatically if it doesn't exist
if not os.path.exists(samples_file):
    with open(samples_file, "w") as f:
        for fq in fastq_files:
            base = os.path.basename(fq)
            match = re.match(r"(.+)_R[12].*\.fastq", base)
            if match:
                seq_id = match.group(1)
                # Customize short name as needed; here using last underscore field
                short = seq_id.split("_")[-1]
                if short not in samples:
                    samples[short] = seq_id
                    f.write(f"{seq_id} {short}\n")
else:
    # Read existing samples_file
    with open(samples_file) as f:
        for line in f:
            seq_id, short = line.strip().split()
            samples[short] = seq_id

# Make a list of all sample short names
SAMPLES = list(samples.keys())

########
# RULES

# Create per-sample directories
rule make_sample_dirs:
    output:
        directory(os.path.join(clean, "{sample}"))

# AdapterRemoval
rule adapterremoval:
    input:
        r1=lambda wc: glob.glob(os.path.join(data, f"{samples[wc.sample]}*_R1*.fastq"))[0],
        r2=lambda wc: glob.glob(os.path.join(data, f"{samples[wc.sample]}*_R2*.fastq"))[0],
        outdir=lambda wc: os.path.join(clean, wc.sample)
    output:
        r1=os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq"),
        r2=os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq"),
        singleton=os.path.join(clean, "{sample}", "{sample}_singleton_truncated.fastq"),
        collapsed=os.path.join(clean, "{sample}", "{sample}_collapsed.fastq"),
        collapsed_trunc=os.path.join(clean, "{sample}", "{sample}_collapsed_truncated.fastq"),
        discarded=os.path.join(clean, "{sample}", "{sample}_discarded.fastq"),
        settings=os.path.join(clean, "{sample}", "{sample}.settings")
    shell:
        """
        module load tools
        module load gcc/15.1.0
        module load adapterremoval/2.3.1

        AdapterRemoval \
          --file1 {input.r1} \
          --file2 {input.r2} \
          --settings {output.settings} \
          --output1 {output.r1} \
          --output2 {output.r2} \
          --singleton {output.singleton} \
          --outputcollapsed {output.collapsed} \
          --outputcollapsedtruncated {output.collapsed_trunc} \
          --discarded {output.discarded} \
          --trimns \
          --trimqualities \
          --minquality 20 \
          --minlength 25 \
          --maxns 20 \
          --collapse
        """

# Final "all" rule
rule all:
    input:
        expand(os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq"), sample=SAMPLES)

