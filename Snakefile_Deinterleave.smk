# Deinterleave FASTQ files

## Software dependencies:
# - BBMap

import glob
import re
import os
from snakemake.utils import R, report


files, = glob_wildcards("00_Denoized/{sample}.fq.gz")

rule all:
    input:
        expand("01_Deint/{sample}__1.fq.gz", zip, sample = files),
        expand("01_Deint/{sample}__2.fq.gz", zip, sample = files)

rule deint:
    input:
        "00_Denoized/{sample}.fq.gz"
    output:
        R1 = "01_Deint/{sample}__1.fq.gz",
        R2 = "01_Deint/{sample}__2.fq.gz"
    log:
        "logs/{sample}.log"
    params:
        compress=7
    threads: 4
    conda:
        "envs/bbmap.yaml"
    shadow: "shallow"
    message:
        "Deinterleaving {wildcards.sample}."
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
