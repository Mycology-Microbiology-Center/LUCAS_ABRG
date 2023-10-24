## Deep learning-based prediction of antibiotic resistance genes
## with DeepARG (short_reads_pipeline)

## Author - Vladimir Mikryukov

## NB. DeepARG does not restrict number of threads for Diamond,
##     `deeparg/entry.py` was modified and Diamond threads were hard-coded to 16 !

## Software dependencies: (see `envs/deeparg.yaml` and `rclone.yaml`)
# - deeparg == 1.0.2  (via pip)
# - diamond==0.9.24
# - python=2.7.18
# - trimmomatic
# - vsearch>=2.17.1
# - bedtools>=2.30.0
# - bowtie2>=2.3.5
# - tbb=2020.2
# - samtools>=1.12

import glob
import re
import os
from snakemake.utils import R, report


directories, files, = glob_wildcards("99_Links/{dir}/{sample}__1.txt")

rule all:
    input:
        expand("05_DeepARG/{dir}/{sample}.clean.fa.gz", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}.sam.gz", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.align.daa", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.align.daa.tsv.gz", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged.quant", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged.quant.subtype", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged.quant.type", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_deeparg.mapping.potential.ARG", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_sorted.bam", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_sorted.bam.merged", zip, dir = directories, sample = files),
        expand("05_DeepARG/{dir}/{sample}_sorted.bam.merged.quant", zip, dir = directories, sample = files)

## Download data
rule download:
    input:
        R1 = "99_Links/{dir}/{sample}__1.txt"
    output:
        R1 = temp("01_Deint/{dir}/{sample}__1.fq.gz"),
        R2 = temp("01_Deint/{dir}/{sample}__2.fq.gz")
    log:
        R1 = "logs/00_Download/{dir}/{sample}_R1.log",
        R2 = "logs/00_Download/{dir}/{sample}_R2.log"
    threads: 2
    resources:
        bigfile = 1,
        mem_mb = 400,
        runtime = 30
    conda:
        "envs/rclone.yaml"
    shadow: "shallow"
    message:
        "Downloading {wildcards.sample}."
    shell: """

    ## Download R1
    rclone copy \
      --multi-thread-streams {threads} \
      UT_OneDrive:/01_Deint/{wildcards.dir}/{wildcards.sample}__1.fq.gz \
      ./01_Deint/{wildcards.dir}/ \
      --log-file {log.R1}

    if [ $? -eq 0 ]
    then
      touch {output.R1}
      echo "R1 file downloaded" >> {log.R1}
    else
      echo "Could not download R1 file" >> {log.R1}
    fi

    ## Download R2
    rclone copy \
      --multi-thread-streams {threads} \
      UT_OneDrive:/01_Deint/{wildcards.dir}/{wildcards.sample}__2.fq.gz \
      ./01_Deint/{wildcards.dir}/ \
      --log-file {log.R2}

    if [ $? -eq 0 ]
    then
      touch {output.R2}
      echo "R2 file downloaded" >> {log.R2}
    else
      echo "Could not download R2 file" >> {log.R2}
    fi

    """


## ARG prediction
rule deeparg:
    input:
        R1 = rules.download.output.R1,
        R2 = rules.download.output.R2
    output:
        reads_clean = "05_DeepARG/{dir}/{sample}.clean.fa.gz",
        sam = "05_DeepARG/{dir}/{sample}.sam.gz",
        daa = "05_DeepARG/{dir}/{sample}_deeparg.align.daa",
        daa_tsv = "05_DeepARG/{dir}/{sample}_deeparg.align.daa.tsv.gz",
        ARG = "05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG",
        ARG_merg = "05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged",
        ARG_merg_quant = "05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged.quant",
        ARG_merg_quant_subtype = "05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged.quant.subtype",
        ARG_merg_quant_type = "05_DeepARG/{dir}/{sample}_deeparg.mapping.ARG.merged.quant.type",
        ARG_potent = "05_DeepARG/{dir}/{sample}_deeparg.mapping.potential.ARG",
        sort_bam = "05_DeepARG/{dir}/{sample}_sorted.bam",
        sort_bam_mrg = "05_DeepARG/{dir}/{sample}_sorted.bam.merged",
        sort_bam_mrg_quant = "05_DeepARG/{dir}/{sample}_sorted.bam.merged.quant"
    log:
        "logs/05_DeepARG/{dir}/{sample}.log",
    params:
        DBPATH = "/gpfs/space/home/amiri/DeepARG/",
        identity = 80,
        probability = 0.8,
        evalue = 1e-10,
        ident16s = 0.75
    threads: 16
    resources:
        mem_mb = 13000,
        runtime = 180
    shadow: "shallow"
    conda:
        "envs/deeparg.yaml"
    message:
        "DeepARG - {wildcards.sample}."
    shell: """

    mkdir -p tmp_{wildcards.dir}

    echo "Starting deeparg" >> {log}
    time \
    deeparg short_reads_pipeline \
        --forward_pe_file {input.R1} \
        --reverse_pe_file {input.R2} \
        --output_file tmp_{wildcards.dir}/reads \
        --deeparg_data_path {params.DBPATH} \
        --deeparg_identity {params.identity} \
        --deeparg_probability {params.probability} \
        --deeparg_evalue {params.evalue} \
        --bowtie_16s_identity {params.ident16s} \
        &>> {log}

    echo "Compressing output" >> {log}
    gzip -5 tmp_{wildcards.dir}/reads.clean &
    gzip -5 tmp_{wildcards.dir}/reads.clean.sam &
    gzip -5 tmp_{wildcards.dir}/reads.clean.deeparg.align.daa.tsv &
    wait

    echo "Moving output to the results folder" >> {log}
    mv tmp_{wildcards.dir}/reads.clean.gz {output.reads_clean}
    mv tmp_{wildcards.dir}/reads.clean.sam.gz {output.sam}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.align.daa.tsv.gz {output.daa_tsv}

    mv tmp_{wildcards.dir}/reads.clean.deeparg.align.daa {output.daa}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.mapping.ARG {output.ARG}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.mapping.ARG.merged {output.ARG_merg}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.mapping.ARG.merged.quant {output.ARG_merg_quant}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.mapping.ARG.merged.quant.subtype {output.ARG_merg_quant_subtype}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.mapping.ARG.merged.quant.type {output.ARG_merg_quant_type}
    mv tmp_{wildcards.dir}/reads.clean.deeparg.mapping.potential.ARG {output.ARG_potent}
    mv tmp_{wildcards.dir}/reads.clean.sorted.bam {output.sort_bam}
    mv tmp_{wildcards.dir}/reads.clean.sorted.bam.merged {output.sort_bam_mrg}
    mv tmp_{wildcards.dir}/reads.clean.sorted.bam.merged.quant {output.sort_bam_mrg_quant}

    ## Clean up
    echo "Cleaning" >> {log}
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__1.fq.gz.paired
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__1.fq.gz.paired.merged
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__1.fq.gz.paired.unmerged
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__1.fq.gz.unpaired
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__2.fq.gz.paired
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__2.fq.gz.paired.unmerged
    rm 01_Deint/{wildcards.dir}/{wildcards.sample}__2.fq.gz.unpaired

    echo "Sample" {wildcards.dir} "done" >> {log}
    """
