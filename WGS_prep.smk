# ==========================================================
# WGS Data Preparation and Adapter Trimming Pipeline
# ==========================================================

import os
import glob
import re

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

########
# Collect FASTQ files (gzipped or not)
fastq_files = [f for f in glob.glob(os.path.join(data, "*.fastq*"))
               if not os.path.basename(f).startswith("._")]

########
# Build samples dictionary: short_name -> base seq_id
samples = {}

# Generate samples_file if missing
if not os.path.exists(samples_file):
    with open(samples_file, "w") as f:
        for fq in fastq_files:
            base = os.path.basename(fq)
            # Match and strip _R1/_R2 and _001 etc.
            match = re.match(r"(.+)_R[12]_.*\.fastq(?:\.gz)?$", base)
            if match:
                seq_id = match.group(1)
                if seq_id not in samples:
                    samples[seq_id] = seq_id
                    f.write(f"{seq_id}\t{seq_id}\n")
else:
    # Read existing samples_file
    with open(samples_file) as f:
        for line in f:
            seq_id, short = line.strip().split()
            samples[short] = seq_id

# List of sample names for expand()
SAMPLES = list(samples.keys())

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
