# ==========================================================
# Genotype Likelihood–Based Population Genomics Pipeline
# ==========================================================

import os

# Load config
configfile: "/home/krinem/UKCEH_Ecological_Genetics/YAML/angsdPopGen.yaml"

# Paths
ENV_DIR = "/ssd0/krinem/UKCEH_Ecological_Genetics/YAML"
benchmark_dir = "benchmarks"

bamlist = config["bamlist"]
REF_FASTA = config["reference"]["fasta"]

# Ensure log directory exists
os.makedirs("logs/angsd", exist_ok=True)
os.makedirs("logs/ngsadmix", exist_ok=True)
os.makedirs("logs/pcangsd", exist_ok=True)

# ----------------------------------------------------------
# Rule: all
# ----------------------------------------------------------

rule all:
    input:
        config["angsd"]["out"] + ".beagle.gz",
        expand("analyses/ngsadmix/admix_K{K}-1.qopt",
               K=config["ngsadmix"]["K"]),
        config["pcangsd"]["out_prefix"] + ".cov.npy"

# ----------------------------------------------------------
# 1. ANGSD
# ----------------------------------------------------------

rule angsd:
    input:
        bamlist=bamlist,
        ref=REF_FASTA
    output:
        geno=config["angsd"]["out"] + ".beagle.gz"
    params:
        out=config["angsd"]["out"],
        a=config["angsd"]
    threads: config["angsd"]["threads"]
    log:
        "logs/angsd/angsd.log"
    conda:
        os.path.join(ENV_DIR, "angsd.yaml")
    shell:
        """
        angsd \
            -b {input.bamlist} \
            -ref {input.ref} \
            -out {params.out} \
            -uniqueOnly {params.a[uniqueOnly]} \
            -remove_bads {params.a[remove_bads]} \
            -trim {params.a[trim]} \
            -C {params.a[C]} \
            -baq {params.a[baq]} \
            -minMapQ {params.a[minMapQ]} \
            -minQ {params.a[minQ]} \
            -docounts {params.a[docounts]} \
            -setMaxDepthInd {params.a[setMaxDepthInd]} \
            -gl {params.a[gl]} \
            -domajorminor {params.a[doMajorMinor]} \
            -domaf {params.a[doMaf]} \
            -doglf {params.a[doGlf]} \
            -dopost {params.a[doPost]} \
            -SNP_pval {params.a[SNP_pval]} \
            -dogeno {params.a[doGeno]} \
            --ignore-RG {params.a[ignore_RG]} \
            -geno_minDepth {params.a[geno_minDepth]} \
            -geno_maxDepth {params.a[geno_maxDepth]} \
            -postCutoff {params.a[postCutoff]} \
            -P {threads} \
            > {log} 2>&1
        """

# ----------------------------------------------------------
# 2. NGSadmix
# ----------------------------------------------------------

rule ngsadmix:
    input:
        geno=rules.angsd.output.geno
    output:
        "analyses/ngsadmix/admix_K{K}-1.qopt"
    params:
        threads=config["ngsadmix"]["threads"]
    threads: config["ngsadmix"]["threads"]
    log:
        "logs/ngsadmix/K{K}.log"
    conda:
        os.path.join(ENV_DIR, "angsd.yaml")
    shell:
        """
        NGSadmix \
            -likes {input.geno} \
            -K {wildcards.K} \
            -P {threads} \
            -outfiles analyses/ngsadmix/admix_K{wildcards.K}-1 \
            > {log} 2>&1
        """

# ----------------------------------------------------------
# 3. PCAngsd
# ----------------------------------------------------------

rule pcangsd:
    input:
        geno=rules.angsd.output.geno
    output:
        config["pcangsd"]["out_prefix"] + ".cov.npy"
    params:
        out=config["pcangsd"]["out_prefix"],
        threads=config["pcangsd"]["threads"]
    threads: config["pcangsd"]["threads"]
    log:
        "logs/pcangsd/pcangsd.log"
    conda:
        os.path.join(ENV_DIR, "pcangsd.yaml")
    shell:
        """
        pcangsd \
            -b {input.geno} \
            -o {params.out} \
            -t {threads} \
            > {log} 2>&1
        """