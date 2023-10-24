## Extract best hits from DIAMOND matches (EggNog database) 
## and find unique target IDs that will be used for eggnog-mapper

## Author - Vladimir Mikryukov

## Software dependencies:
# - R >= 4.0.0
# - data.table >= 1.14.0 
# - stringi >= 1.5.3
# - optparse >= 1.6.6
# - R.utils >= 2.10.1


import glob
import re
import os
from snakemake.utils import R, report


files, = glob_wildcards("00_Diamond_AllHits/{sample}.m8.gz")

rule all:
    input:
        expand("01_Diamond_BestHits/{sample}.m8.gz", zip, sample = files),

## Extract best hits
rule get_best_hits:
    input:
        "00_Diamond_AllHits/{sample}.m8.gz"
    output:
        "01_Diamond_BestHits/{sample}.m8.gz"
    log:
        "logs/01_Diamond_BestHits/{sample}.log"
    threads: 1
    params:
        pident=50,
        evalue=1.0e-8
    message:
        "Extracting best hits from {wildcards.sample}."
    shell: """

    scripts/Extract_Diamond_hits.R \
        --input {input} \
        --pident {params.pident} \
        --evalue {params.evalue} \
        --output {output} \
        > {log}

    """
## ~ 50 seconds and ~3.6 GB RAM on 1 core for 15.3 mln records (~280 MB)

