## Taxonomy profiling with Kraken2 and abundance estimation with Bracken

## Author - Vladimir Mikryukov

## Software dependencies:
# - kraken2 >= 2.1.2
# - bracken >= 2.6.1
# - BBMap

## Input = interleaved filtered reads (the same as used for Diamond)

##### Imports

import glob
import re
from snakemake.utils import R, report


## Input = interleaved paired-end reads
directories, files = glob_wildcards("03_Chunks/{dir}/{sample}.fasta.gz")

## Dummy rule to collect all targets (this will only be executed if no other target is defined)
rule all:
    input:
        expand("07_Kraken/{dir}/{dir}__KrakRep.log", zip, dir = directories),
        expand("07_Kraken/{dir}/{dir}__KrakOut.txt.gz", zip, dir = directories),
        expand("07_Kraken/{dir}/{dir}__Bracken_species.txt", zip, dir = directories),
        expand("07_Kraken/{dir}/{dir}__Bracken_family.txt", zip, dir = directories),
        expand("07_Kraken/{dir}/{dir}__Bracken_phylum.txt", zip, dir = directories),

## Merge chunks
rule merge_chunks:
    output:
        temp("99_cat/{dir}/{dir}.fasta.gz")
    threads: 1
    message:
        "Concatenating sample - {wildcards.dir}."
    shell: """
        mkdir -p 99_cat/{wildcards.dir}
        cat 03_Chunks/{wildcards.dir}/*.fasta.gz > {output}
    """

## Deinterleave reads
rule deinterleave:
    input:
        rules.merge_chunks.output,
    output:
        R1 = temp("99_Deint/{dir}/{dir}__1.fa.gz"),
        R2 = temp("99_Deint/{dir}/{dir}__2.fa.gz")
    threads: 5
    conda:
        "envs/bbmap.yaml"
    log:
        "logs/{dir}/1_deint.log"
    params:
        compress=4
    # shadow: "shallow"
    message:
        "Deinterleaving {wildcards.dir}."
    shell: """

    reformat.sh \
      in={input} \
      int=t \
      out1={output.R1} \
      out2={output.R2} \
      zl={params.compress} \
      t={threads} \
      &> {log}

    """


## Generate profiles of microbial clades and their abundances with Kraken2
rule kraken:
    input:
        R1 = rules.deinterleave.output.R1,
        R2 = rules.deinterleave.output.R2
    output:
        KReport = "07_Kraken/{dir}/{dir}__KrakRep.log",
        KRakOut = "07_Kraken/{dir}/{dir}__KrakOut.txt.gz"
    shadow: "shallow"
    conda:
        "envs/kraken2.yaml"
    threads: 5
    log:
        "logs/{dir}/2_Kraken.log"
    params:
        DB="/DB/k2_pluspf_20230605/"
    message:
        "Kraken - sample {wildcards.dir}."
    shell: """

    kraken2 \
      --db {params.DB} \
      --use-names \
      --threads {threads} \
      --report {output.KReport} \
      --gzip-compressed \
      --paired {input.R1} {input.R2} \
      2> {log} \
      | gzip > {output.KRakOut}

    """

## Re-estimate reads assigned to species with Bracken
rule bracken_species:
    input:
        KReport = rules.kraken.output.KReport
    output:
        Brak = "07_Kraken/{dir}/{dir}__Bracken_species.txt"
    shadow: "shallow"
    conda:
        "envs/kraken2.yaml"
    threads: 1
    params:
        DB="/DB/k2_pluspf_20230605/",
        Threshold=5,
        ReadLen=150,
        Level="S"
    log:
        "logs/{dir}/3_Bracken_species.log"
    message:
        "Bracken - species - sample {wildcards.dir}."
    shell: """

    bracken \
      -d {params.DB} \
      -i {input.KReport} \
      -l {params.Level} \
      -t {params.Threshold} \
      -r {params.ReadLen} \
      -o {output.Brak} \
      >> {log}

    """

## Re-estimate reads assigned to family with Bracken
rule bracken_family:
    input:
        KReport = rules.kraken.output.KReport
    output:
        Brak = "07_Kraken/{dir}/{dir}__Bracken_family.txt"
    shadow: "shallow"
    conda:
        "envs/kraken2.yaml"
    threads: 1
    params:
        DB="/DB/k2_pluspf_20230605/",
        Threshold=5,
        ReadLen=150,
        Level="F"
    log:
        "logs/{dir}/4_Bracken_family.log"
    message:
        "Bracken - family - sample {wildcards.dir}."
    shell: """

    bracken \
      -d {params.DB} \
      -i {input.KReport} \
      -l {params.Level} \
      -t {params.Threshold} \
      -r {params.ReadLen} \
      -o {output.Brak} \
      >> {log}

    """

## Re-estimate reads assigned to species with Bracken
rule bracken_phylum:
    input:
        KReport = rules.kraken.output.KReport
    output:
        Brak = "07_Kraken/{dir}/{dir}__Bracken_phylum.txt"
    shadow: "shallow"
    conda:
        "envs/kraken2.yaml"
    threads: 1
    params:
        DB="/DB/k2_pluspf_20230605/",
        Threshold=5,
        ReadLen=150,
        Level="P"
    log:
        "logs/{dir}/4_Bracken_phylum.log"
    message:
        "Bracken - phylum - sample {wildcards.dir}."
    shell: """

    bracken \
      -d {params.DB} \
      -i {input.KReport} \
      -l {params.Level} \
      -t {params.Threshold} \
      -r {params.ReadLen} \
      -o {output.Brak} \
      >> {log}

    """
