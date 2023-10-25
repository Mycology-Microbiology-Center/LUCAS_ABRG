## Identification and taxonomic classification of small and large subunit rRNA in metagenomic data
## with metaxa2

## Author - Vladimir Mikryukov

## Software dependencies: (see `envs/metaxa.yaml` and `envs/rclone.yaml`)
# - metaxa>=2.2
# - blast-legacy
# - perl
# - bbmap
# - rclone

import glob
import re
import os
from snakemake.utils import R, report


files, = glob_wildcards("99_Links/{sample}.txt")

rule all:
  input:
    expand("02_Metaxa/{sample}/{sample}_SSU_eukaryota.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_archaea.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_bacteria.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_uncertain.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_hmmmer.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_blast.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_extraction.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_taxonomy.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_SSU_taxonomy_reliability.txt.gz", zip, sample = files),

    expand("02_Metaxa/{sample}/{sample}_LSU_eukaryota.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_archaea.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_bacteria.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_uncertain.fa.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_hmmmer.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_blast.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_extraction.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_taxonomy.txt.gz", zip, sample = files),
    expand("02_Metaxa/{sample}/{sample}_LSU_taxonomy_reliability.txt.gz", zip, sample = files)



## Download data
rule download:
  input:
      RR = "99_Links/{sample}.txt"
  output:
      R1 = temp("00_Download/{sample}/{sample}_R1.fasta"),
      R2 = temp("00_Download/{sample}/{sample}_R2.fasta")
  log:
    "logs/00_Download/{sample}.log"
  threads: 2
  resources:
    bigfile = 1,
    mem_mb = 1000,
    runtime = 20
  conda:
    "envs/rclone.yaml"
  shadow: "minimal"
  message:
    "Downloading {wildcards.sample}."
  shell: """

  ## Download data
  echo "Downloading: " {wildcards.sample} > {log}

  rclone copy \
    --multi-thread-streams {threads} \
    UT_Nilou:/Data/{wildcards.sample}.fasta.gz \
    ./00_Download/{wildcards.sample}/ \
    --log-file {log}

  if [ $? -eq 0 ]
  then
    echo "File downloaded" >> {log}
  else
    echo "Could not download file" >> {log}
  fi

  echo "Deinterleaving reads" >> {log}

  ## Deinterleave reads
  ## Split data into a paired reads
  ## (do not compress FASTA, as metaxa does not accept compressed reads)
  reformat.sh \
    int=t fastawrap=999 \
    in="00_Download/{wildcards.sample}/{wildcards.sample}.fasta.gz" \
    out1={output.R1} out2={output.R2} \
    2>> {log}

  ## Clean up
  rm "00_Download/{wildcards.sample}/{wildcards.sample}.fasta.gz"

  """



## SSU identification
rule metaxa_ssu:
  input:
      R1 = rules.download.output.R1,
      R2 = rules.download.output.R2
  output:
      ssu_eukaryota = "02_Metaxa/{sample}/{sample}_SSU_eukaryota.fa.gz",
      ssu_archaea = "02_Metaxa/{sample}/{sample}_SSU_archaea.fa.gz",
      ssu_bacteria = "02_Metaxa/{sample}/{sample}_SSU_bacteria.fa.gz",
      ssu_uncertain = "02_Metaxa/{sample}/{sample}_SSU_uncertain.fa.gz",
      ssu_hmmmer = "02_Metaxa/{sample}/{sample}_SSU_hmmmer.txt.gz",
      ssu_blast = "02_Metaxa/{sample}/{sample}_SSU_blast.txt.gz",
      ssu_extract = "02_Metaxa/{sample}/{sample}_SSU_extraction.txt.gz",
      ssu_tax = "02_Metaxa/{sample}/{sample}_SSU_taxonomy.txt.gz",
      ssu_taxrel = "02_Metaxa/{sample}/{sample}_SSU_taxonomy_reliability.txt.gz"
  log:
      log1 = "logs/02_Metaxa/{sample}_SSU.log",
      log2 = "logs/02_Metaxa/{sample}_SSU_summary.log"
  threads: 1
  resources:
      mem_mb = 1000,
      runtime = 30
  shadow: "minimal"
  conda:
    "envs/metaxa.yaml"
  message:
    "Metaxa SSU - {wildcards.sample}."
  shell: """

  #### Run Metaxa (+ use BLAST+ instead of legacy BLAST)
  
  echo "Starting Metaxa" >> {log.log1}
  time \
  metaxa2 \
    -1 {input.R1} -2 {input.R2} \
    --mode metagenome -t all \
    -E 1 -S 12 -N 2 -M 5 -H5 -R 80 \
    -T 0,60,70,75,85,90,97 \
    --selection_priority score \
    --blast_eval 1e-10 \
    --blast_wordsize 14 \
    --allow_single_domain 1e-9,0 \
    --complement T \
    --plus T --megablast F \
    --heuristics T \
    --summary T --graphical F --fasta T --split_pairs F \
    --table T --taxonomy T --reltax T \
    --not_found F \
    --align none --truncate T \
    -g ssu \
    --cpu {threads} --multi_thread F \
    -o SSU \
    2>> {log.log1}

    # --cpu {threads} --multi_thread T

  echo "Moving output to the results folder" >> {log.log1}
  cat SSU.summary.txt >> {log.log2}

  gzip -c SSU.eukaryota.fasta > {output.ssu_eukaryota}
  gzip -c SSU.archaea.fasta   > {output.ssu_archaea}
  gzip -c SSU.bacteria.fasta  > {output.ssu_bacteria}
  gzip -c SSU.uncertain.fasta > {output.ssu_uncertain}

  gzip -c SSU.hmmer.table        > {output.ssu_hmmmer}
  gzip -c SSU.blast.table        > {output.ssu_blast}
  gzip -c SSU.extraction.results > {output.ssu_extract}

  gzip -c SSU.taxonomy.txt > {output.ssu_tax}
  gzip -c SSU.taxonomy-reliability.txt > {output.ssu_taxrel}
  

  ## Redundant output
  # SSU.extraction.fasta     # all sequences
  # SSU.mitochondria.fasta
  # SSU.chloroplast.fasta

  echo "Metaxa SSU for sample" {wildcards.sample} "finished" >> {log.log1}
  """



## LSU identification
rule metaxa_lsu:
  input:
      R1 = rules.download.output.R1,
      R2 = rules.download.output.R2
  output:
      lsu_eukaryota = "02_Metaxa/{sample}/{sample}_LSU_eukaryota.fa.gz",
      lsu_archaea = "02_Metaxa/{sample}/{sample}_LSU_archaea.fa.gz",
      lsu_bacteria = "02_Metaxa/{sample}/{sample}_LSU_bacteria.fa.gz",
      lsu_uncertain = "02_Metaxa/{sample}/{sample}_LSU_uncertain.fa.gz",
      lsu_hmmmer = "02_Metaxa/{sample}/{sample}_LSU_hmmmer.txt.gz",
      lsu_blast = "02_Metaxa/{sample}/{sample}_LSU_blast.txt.gz",
      lsu_extract = "02_Metaxa/{sample}/{sample}_LSU_extraction.txt.gz",
      lsu_tax = "02_Metaxa/{sample}/{sample}_LSU_taxonomy.txt.gz",
      lsu_taxrel = "02_Metaxa/{sample}/{sample}_LSU_taxonomy_reliability.txt.gz"
  log:
      log1 = "logs/02_Metaxa/{sample}_LSU.log",
      log2 = "logs/02_Metaxa/{sample}_LSU_summary.log"
  threads: 1
  resources:
      mem_mb = 1000,
      runtime = 30
  shadow: "minimal"
  conda:
    "envs/metaxa.yaml"
  message:
    "Metaxa LSU - {wildcards.sample}."
  shell: """

  #### Run Metaxa (+ use BLAST+ instead of legacy BLAST)
  
  echo "Starting Metaxa" >> {log.log1}
  time \
  metaxa2 \
    -1 {input.R1} -2 {input.R2} \
    --mode metagenome -t all \
    -E 1 -S 12 -N 2 -M 5 -H5 -R 80 \
    -T 0,60,70,75,85,90,97 \
    --selection_priority score \
    --blast_eval 1e-10 \
    --blast_wordsize 14 \
    --allow_single_domain 1e-9,0 \
    --complement T \
    --plus T --megablast F \
    --heuristics T \
    --summary T --graphical F --fasta T --split_pairs F \
    --table T --taxonomy T --reltax T \
    --not_found F \
    --align none --truncate T \
    -g lsu \
    --cpu {threads} --multi_thread F \
    -o LSU \
    2>> {log.log1}

    # --cpu {threads} --multi_thread T

  echo "Moving output to the results folder" >> {log.log1}
  cat LSU.summary.txt >> {log.log2}

  gzip -c LSU.eukaryota.fasta > {output.lsu_eukaryota}
  gzip -c LSU.archaea.fasta   > {output.lsu_archaea}
  gzip -c LSU.bacteria.fasta  > {output.lsu_bacteria}
  gzip -c LSU.uncertain.fasta > {output.lsu_uncertain}

  gzip -c LSU.hmmer.table        > {output.lsu_hmmmer}
  gzip -c LSU.blast.table        > {output.lsu_blast}
  gzip -c LSU.extraction.results > {output.lsu_extract}

  gzip -c LSU.taxonomy.txt > {output.lsu_tax}
  gzip -c LSU.taxonomy-reliability.txt > {output.lsu_taxrel}
  

  ## Redundant output
  # LSU.extraction.fasta     # all sequences
  # LSU.mitochondria.fasta
  # LSU.chloroplast.fasta

  echo "Metaxa LSU for sample" {wildcards.sample} "finished" >> {log.log1}
  """

