# ----------------------------------------------------------
# Load config
# ----------------------------------------------------------

import os

configfile: "YAML/bcftoolsconfig.yaml"

OUTPREFIX = config["outprefix"]
POPDIR = config["popdir"]
ENV_DIR = "/ssd0/krinem/UKCEH_Ecological_Genetics/YAML"  # Conda environments directory

# ----------------------------------------------------------
# Detect all population files
# ----------------------------------------------------------

POPLIST = glob_wildcards(f"{POPDIR}" + "/{pop}.txt").pop

import itertools
PAIRS = list(itertools.combinations(POPLIST, 2))

# ----------------------------------------------------------
# Final target: all pairwise FST outputs
# ----------------------------------------------------------

rule all:
    input:
        expand("fst/{pop1}_vs_{pop2}.fst",
               pop1=[p[0] for p in PAIRS],
               pop2=[p[1] for p in PAIRS])

# ----------------------------------------------------------
# 1. mpileup
# ----------------------------------------------------------

rule mpileup:
    input:
        bamlist = config["bamlist"],
        ref     = config["reference"]["fasta"],
        fai     = config["reference"]["fai"]
    output:
        OUTPREFIX + ".mpileup.bcf"
    log:
        "logs/mpileup.log"
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    shell:
        """
        bcftools mpileup -Ou -f {input.ref} -b {input.bamlist} \
            | bcftools view -Ob -o {output} \
            > {log} 2>&1
        """

# ----------------------------------------------------------
# 2. Variant calling
# ----------------------------------------------------------

rule call:
    input:
        OUTPREFIX + ".mpileup.bcf"
    output:
        OUTPREFIX + ".raw.bcf"
    log:
        "logs/call.log"
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    shell:
        """
        bcftools call -mv -Ob -o {output} {input} \
            > {log} 2>&1
        """

# ----------------------------------------------------------
# 3. Filtering
# ----------------------------------------------------------

rule filter:
    input:
        OUTPREFIX + ".raw.bcf"
    output:
        OUTPREFIX + ".filtered.vcf.gz"
    log:
        "logs/filter.log"
    conda:
        os.path.join(ENV_DIR, "bcftools.yaml")
    shell:
        """
        bcftools filter -Oz -o {output} {input} \
            > {log} 2>&1
        tabix -p vcf {output}
        """

# ----------------------------------------------------------
# 4. Pairwise FST (vcftools)
# ----------------------------------------------------------

rule fst:
    input:
        vcf = OUTPREFIX + ".filtered.vcf.gz",
        pop1 = POPDIR + "/{pop1}.txt",
        pop2 = POPDIR + "/{pop2}.txt"
    output:
        "fst/{pop1}_vs_{pop2}.fst"
    log:
        "logs/fst/{pop1}_vs_{pop2}.log"
    conda:
        os.path.join(ENV_DIR, "vcftools.yaml")
    shell:
        """
        vcftools --gzvcf {input.vcf} \
            --weir-fst-pop {input.pop1} \
            --weir-fst-pop {input.pop2} \
            --out fst/{wildcards.pop1}_vs_{wildcards.pop2} \
            > {log} 2>&1
        """