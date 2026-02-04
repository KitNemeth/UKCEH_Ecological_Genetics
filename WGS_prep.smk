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
# Mount
# --- Check that the CIFS share is mounted ---
if not os.path.ismount("/home/krinem/mount/cifs_07793_newLEAF"):
    raise RuntimeError("CIFS share is not mounted. Run 'mountp_newleaf' first.")

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

# Ensure directories exist 
os.makedirs(clean, exist_ok=True)

# Match both uncompressed and gzipped FASTQ
fastq_files = [f for f in glob.glob(os.path.join(data, "*.fastq*")) 
               if not os.path.basename(f).startswith("._")]


# Build sample dictionary: short_name -> full seq_id
samples = {}

# Generate samples_file if missing
if not os.path.exists(samples_file):
    with open(samples_file, "w") as f:
        for fq in fastq_files:
            base = os.path.basename(fq)
            match = re.match(r"(.+)_R[12].*\.fastq(?:\.gz)?$", base)
            if match:
                seq_id = match.group(1)
                # Take first N fields as short name (customize as needed)
                short = "_".join(seq_id.split("_")[:5])
                if short not in samples:
                    samples[short] = seq_id
                    f.write(f"{seq_id}\t{short}\n")
else:
    # Read existing samples_file
    with open(samples_file) as f:
        for line in f:
            seq_id, short = line.strip().split()
            samples[short] = seq_id

# List of sample names for expand()
SAMPLES = list(samples.keys())

########
# RULE ALL
rule all:
    input:
        # AdapterRemoval outputs
        expand(os.path.join(clean, "{sample}", "{sample}_R1_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_R2_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_singleton_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_collapsed.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_collapsed_truncated.fastq.gz"), sample=SAMPLES),
        expand(os.path.join(clean, "{sample}", "{sample}_discarded.fastq.gz"), sample=SAMPLES)

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









