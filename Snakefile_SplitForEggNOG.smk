## Split cleaned and denoized reads (interleaved) into a chunks for DIAMOND and EggNog

## Author - Vladimir Mikryukov

## Software dependencies:
# - seqkit >= 0.15.0
# - BBMap >= 38.87
# - awk, sed, gzip


import glob
import re
from snakemake.utils import R, report


files = glob_wildcards("02_Denoized_Pool/{sample}.fq.gz")

rule all:
    input:
      expand("03_Chunks/{sample}/{sample}.part_{{part}}.fasta.gz", zip, sample = files)

# Convert reads to fasta and split them into chunks (for EggNog)
# use dynamic files
## TO DO - better to use Snakemake's data-dependent conditional execution mechanism
rule split_to_chunks:
    input:
        RR = "02_Denoized_Pool/{sample}.fq.gz"
    output:
        dynamic("03_Chunks/{sample}/{sample}.part_{part}.fasta.gz")
        # directory("03_Chunks/{sample}")
    log:
        "logs/03_Chunks/{sample}.{part}.log"
    threads: 3
    conda:
        "envs/split.yaml"
    params:
        NREADS=1500000
    shadow: "shallow"
    message:
        "Splitting sample {wildcards.sample}."
    shell: """

    ## Rename reads, convert FASTQ to FASTA
    rename.sh in={input.RR} int=t prefix={wildcards.sample} out=stdout.fq 2> {log} \
      | seqkit fq2fa -w 0 2>> {log} \
      | sed '/^>/ s/ /_/g' \
      | gzip -3 > {wildcards.sample}.fasta.gz \

    ## Split into chunks
    seqkit split -w 0 -s {params.NREADS} {wildcards.sample}.fasta.gz -O 03_Chunks/{wildcards.sample} 2>> {log}

    """
