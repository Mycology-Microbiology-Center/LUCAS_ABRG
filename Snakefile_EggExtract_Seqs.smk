## Extract sequences of the best hits from DIAMOND matches (EggNog database) 

## Author - Vladimir Mikryukov

## Software dependencies:
# - runiq >= 1.2.1
# - ripgrep >= 13.0.0
# - sed
# - awk


import glob
import re
import os
from snakemake.utils import R, report


files, = glob_wildcards("05_EggBest/{sample}.m8.gz")

rule all:
    input:
        expand("05_EggBest_Seqs/{sample}.fasta.gz", zip, sample = files),


## Extract best hit IDs
rule get_best_hit_ids:
    input:
        "05_EggBest/{sample}.m8.gz"
    output:
        temp("tmp_IDS/{sample}.m8")
    threads: 1
    message:
        "Extracting best hit IDs from {wildcards.sample}."
    shell: """

    zcat {input} \
      | awk '{{ print $1 }}' \
      | sed 's/_12$/_/g; s/_1$/_/g; s/_2$/_/g' \
      | runiq -f digest - \
      > {output}

    """


## Extract sequences by IDs
# (sequences should be without line breaks)
rule get_seqs:
    input:
        IDS = rules.get_best_hit_ids.output,
        SQS = "03_Chunks_Pool/{sample}.fasta.gz"
    output:
        "05_EggBest_Seqs/{sample}.fasta.gz"
    threads: 1
    params:
        gzip_flags = "-7"
    message:
        "Extracting sequences from {wildcards.sample}."
    shell: """

    rg -z -A 1 \
      -f {input.IDS} \
      --context-separator "" \
      --threads {threads} \
      {input.SQS} \
      | sed '/^$/d' \
      | gzip {params.gzip_flags} > {output}
    
    """

