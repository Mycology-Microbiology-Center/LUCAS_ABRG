## Data-filtering pipeline

## Author - Vladimir Mikryukov

## Software dependencies:
# - cutadapt >= 3.1
# - fastp >= 0.20.1
# - vsearch >= 2.15.0
# - BBMap >= 38.87
# - awk
# - pigz

##### Variables

# configfile: "config.yaml"

PIPELINEVERSION='1.5'

# Adapters - Illumina adapters used for library preparation
ADAPTSEQ = "CTGTCTCTTATA"

# Minimum 5' end quality - trim that bases
# there are reads with the first nucleotide = N with low quality
MIN_5_Q=12

# Discard reads shorter than MINLEN (at the primer trimming stage)
MINLEN = 35

# Filtering based on percentage of unqualified bases - fastp
# how many percents of bases are allowed to be unqualified (Q < 24)
PHRED_MIN=24     # --qualified_quality_phred
PHRED_PERC=20    # --unqualified_percent_limit

# Read correction settings (based on PE overlap) - for fastp
OVER_MINLEN=40   # --overlap_len_require
OVER_DIFF=5      # --overlap_diff_limit 
OVER_DIFPERC=15  # --overlap_diff_percent_limit

# Poly-G tails trimming - fastp
POLY_G_MINLEN=4

# Max number for expected error for read filtering. Increase to reduce stringency
MAXEE = 3

# Max number of ambiguous nucleotides
MAXN = 3


##### Imports

import glob
import re
from snakemake.utils import R, report


directories, files = glob_wildcards("00_RawData/{dir}/{sample}_1.fq.gz")


## Dummy rule to collect all targets (this will only be executed if no other target is defined)
rule all:
    input:
        expand("02_Denoized/{dir}/{sample}__denoized2.fq.gz", zip, dir = directories, sample = files)


## Trim sequencing adapters from both paired-end reads with cutadapt
## Remove Novaseq poly-G tails
## error correction based on PE overlap
## quality filtering
## Use named pipes and process substitution to reduce IO
rule data_qc:
    input:
        R1 = "00_RawData/{dir}/{sample}_1.fq.gz",
        R2 = "00_RawData/{dir}/{sample}_2.fq.gz"
    output:
        R1 = temp("01_QC/{dir}/{sample}__R1.clean.fq.gz"),
        R2 = temp("01_QC/{dir}/{sample}__R2.clean.fq.gz")
    log:
        log1 = "logs/01_QC/{dir}/1_{sample}_cutadapt.log",
        log2 = "logs/01_QC/{dir}/2_{sample}_fastp.log",
        log2h= "logs/01_QC/{dir}/2_{sample}_fastp.html",
        log2j= "logs/01_QC/{dir}/2_{sample}_fastp.json"
        # , log3 = "logs/01_QC/{dir}/3_{sample}_vsearch.log"
    shadow: "shallow"
    threads: 5
    conda:
        "envs/datafiltering.yaml"
    message:
        "Quality filtering - sample {wildcards.dir}."
    shell: """

    ## Create named pipes
    mkfifo tmp_r1_cutadapt.fq
    mkfifo tmp_r2_cutadapt.fq

    mkfifo tmp_r1_fastp.fq
    mkfifo tmp_r2_fastp.fq


    ## Nextera trimming with cutadapt
    # + remove N from 5' end (in has very low quality)
    cutadapt \
      -a {ADAPTSEQ} -A {ADAPTSEQ} \
      --minimum-length {MINLEN} \
      -q {MIN_5_Q},0 \
      -o tmp_r1_cutadapt.fq -p tmp_r2_cutadapt.fq \
      {input.R1} {input.R2} \
      --cores 1 \
      > {log.log1} &

    ## Remove poly-G tails and correct reads by overlap
    ## draft quality filtering (no more than 20% of nucleotides with Phred < 24)
    fastp \
      --out1=tmp_r1_fastp.fq --out2=tmp_r2_fastp.fq \
      --disable_adapter_trimming \
      --qualified_quality_phred={PHRED_MIN} --unqualified_percent_limit={PHRED_PERC} \
      --length_required={MINLEN} \
      --correction --overlap_len_require={OVER_MINLEN} --overlap_diff_limit={OVER_DIFF} --overlap_diff_percent_limit={OVER_DIFPERC} \
      --trim_poly_g --poly_g_min_len={POLY_G_MINLEN} \
      --thread=1 \
      --html={log.log2h} --json={log.log2j} \
      --stdin --interleaved_in \
      --in1=<(paste tmp_r1_cutadapt.fq tmp_r2_cutadapt.fq | paste - - - - | awk -v OFS="\n" -v FS="\t" '{{print($1,$3,$5,$7,$2,$4,$6,$8)}}') \
      2> {log.log2} &
  
    cat tmp_r1_fastp.fq | gzip -3 > {output.R1} &
    cat tmp_r2_fastp.fq | gzip -3 > {output.R2}

    ## Remove temp files
    rm tmp_r1_cutadapt.fq
    rm tmp_r2_cutadapt.fq
    rm tmp_r1_fastp.fq
    rm tmp_r2_fastp.fq
    """
# first part with named pipes requires ~ 4 minutes per 1GB of gzip-compressed data and ~7 threads (with vsearch)

## Decontamination by mapping
## Remove synthetic artifacts and spike-ins
rule decontamination:
    input:
        R1 = rules.data_qc.output.R1,
        R2 = rules.data_qc.output.R2
    output:
        RR = temp("01_QC/{dir}/{sample}__RR.clean.fq.gz"),
        CC = temp("01_QC/{dir}/{sample}_contaminants.fq.gz")
    log:
        log4 = "logs/01_QC/{dir}/4_{sample}_bbmap.log",
        log5 = "logs/01_QC/{dir}/5_{sample}_bbduk.log",
        log6 = "logs/01_QC/{dir}/6_{sample}_artifacts.log",
        log7 = "logs/01_QC/{dir}/7_{sample}_ihist_overlap.log",
        log8 = "logs/01_QC/{dir}/7_{sample}_overlap_correct.log"
    shadow: "shallow"
    threads: 2
    conda:
        "envs/datafiltering.yaml"
    params:
        Contaminants="/mnt/Dat2/DB/Genomes/BBmap/Contaminants.fa.gz",
        Index="/mnt/Dat2/DB/Genomes/BBmap/",
        JavaMem="-Xmx30g"
    message:
        "Decontamination - sample {wildcards.dir}."
    shell: """

    ## 0. Remove contaminant reads (with high precision and lower sensitivity)
    ## 1. Additionally check for remaining adapters
    ## 2. Remove synthetic artifacts and spike-ins by kmer-matching, allowing 1 mismatch
    ## 3. Correct by overlap - bbmerge (ecco mix vstrict)
    bbmap.sh {params.JavaMem} \
        minratio=0.9 minid=0.94 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
        qtrim=lr trimq=10 untrim \
        idtag printunmappedcount kfilter=25 maxsites=1 k=14 \
        threads=1 \
        ref={params.Contaminants} \
        path={params.Index} \
        in1={input.R1} \
        in2={input.R2} \
        outu=stdout.fq \
        outm={output.CC} \
        2> {log.log4} \
      | bbduk.sh \
        in=stdin.fq int=t \
        outu=stdout.fq \
        threads=1 \
        ktrim=r k=23 mink=11 hdist=1 tbo tpe ref=adapters ftm=5 ordered \
        minlen={MINLEN} \
        maxns={MAXN} \
        2> {log.log5} \
      | bbduk.sh \
        in=stdin.fq int=t \
        out=stdout.fq \
        threads=1 \
        k=31 hdist=1 ref=artifacts,phix ordered cardinality \
        2> {log.log6} \
      | bbmerge.sh \
        in=stdin.fq int=t \
        out={output.RR} \
        threads=1 \
        ecco mix vstrict ordered \
        ihist={log.log7} \
        2> {log.log8}

    """


## Error correction - Phase 1
## Group overlapping reads into clumps (clusters share kmers)
rule error_correct_p1:
    input:
        RR = rules.decontamination.output.RR
    output:
        RR = temp("02_Denoized/{dir}/{sample}__denoized1.fq.gz")
    log:
        log1 = "logs/02_Denoized/{dir}/1_{sample}_clumpify.log"
    threads: 10
    conda:
        "envs/datafiltering.yaml"
    params:
        threads_half=5,
        JavaMem="-Xmx100g"
    shadow: "shallow"
    message:
        "Error correction - Phase 1 - sample {wildcards.dir}."
    shell: """

    ## Error-correct phase 2 - clumpify (ecc passes=4 reorder)
    clumpify.sh {params.JavaMem} \
      in={input.RR} \
      out={output.RR} \
      ecc passes=4 reorder \
      threads={params.threads_half} \
      2> {log.log1}

    """


## Error correction - Phase 2 
## Correction is handled by two algorithms, “pincer” and “tail”.
## Pincer corrects errors bidirectionally, using kmers on the left and right;
## therefore, it can only work on bases in the middle of the read, at least K away from either end.
## Tail is not as robust, but is able to work on the ends of the read. 
rule error_correct_p2:
    input:
        RR = rules.error_correct_p1.output.RR
    output:
        RR = "02_Denoized/{dir}/{sample}__denoized2.fq.gz"
        # R1 = temp("02_Denoized/{dir}/{sample}__R1.fq.gz"),
        # R2 = temp("02_Denoized/{dir}/{sample}__R1.fq.gz")
    log:
        log2 = "logs/02_Denoized/{dir}/2_{sample}_tadpole.log"
    threads: 10
    conda:
        "envs/datafiltering.yaml"
    params:
        JavaMem="-Xmx110g",
        threads_half=7,
        ignore_kmer_depth=1
    shadow: "shallow"
    message:
        "Error correction - Phase 2 - sample {wildcards.dir}."
    shell: """

    ## Error-correct phase 3 - tadpole (ecc k=60 ordered)
    # Low-depth reads can be discarded here with the "tossjunk", "tossdepth", or "tossuncorrectable" flags.
    tadpole.sh {params.JavaMem} \
      in={input.RR} \
      out={output.RR} \
      mode=correct \
      ecc k=60 ordered \
      prefilter={params.ignore_kmer_depth} \
      prealloc \
      threads={params.threads_half} \
      2> {log.log2}

    ## if not enough RAM - use bbcms.sh, less accurate, but scales well
    ## it never runs out of memory, since it uses a lossy data structure

    """
