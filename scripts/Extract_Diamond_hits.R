#!/usr/bin/env Rscript

## Author - Vladimir Mikryukov

# Script for the extraction of DIAMOND best hits (for interleaved pair-end data)
# Input file is a result of DIAMOND homology search in BLAST output format 6 (*.m8.gz)
#    1      qseqid = Query Seq-id
#    2      sseqid = Subject Seq-id
#    3      pident = Percentage of identical matches
#    4      length = Alignment length
#    5    mismatch = Number of mismatches
#    6     gapopen = Number of gap openings
#    7      qstart = Start of alignment in query
#    8        qend = End of alignment in query
#    9      sstart = Start of alignment in subject
#    10       send = End of alignment in subject
#    11     evalue = Expect value
#    12   bitscore = Bit score

## Query ID should contain read pair info (e.g., "QUERYNAME_1:" = R1; "QUERYNAME_2:" = R2)
## read pair info was added with BBMap

## Usage:
# ./Extract_Diamond_hits.R --input "sbs.m8.gz" --pident 50 --evalue 1.0e-8 --output "sbs_best.m8.gz"

## Output:
# Subset with best hits and 4 columns: qseqid, sseqid, evalue, bitscore
# Read pair info is added to the query name:
#   _12 = R1+R2 (the same target for both reads),
#   _1 = R1, _2 = R2 only


###############################
############################### Parse input parameters
###############################

## Script variables that we need to parse
#   INPUTFILE - input file name
#   OUTFILE   - name of the output file with results

cat("Parsing input options and arguments...\n")

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"

## Parse arguments
option_list <- list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character', help="Input file name"),
  make_option(c("-p", "--pident"), action="store", default=50, type='double', help="Minimum percentage of identity"),
  make_option(c("-e", "--evalue"), action="store", default=0.00000001, type='double', help="Maximum E-value"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character', help="Output file")
)
opt <- parse_args(OptionParser(option_list=option_list))

if(is.na(opt$input)){
  cat("Input file is not specified.\n", file=stderr())
  stop()
}
if(is.na(opt$output)){
  cat("Output file is not specified.\n", file=stderr())
  stop()
}

## Assign variables
INPUTFILE <- opt$input
OUTFILE <- opt$output
MINIDENT <- as.numeric(opt$pident)
MAXEVALUE <- as.numeric(opt$evalue)

# MINIDENT <- 50           # Minimum percentage identity of hits
# MAXEVALUE <- 1.0e-8      # Maximum E-value of hits

NCORES_DT <- 1             # Number of CPU cores for data.table

## Log assigned variables
cat(paste("Input file: ", INPUTFILE, "\n", sep=""))
cat(paste("Minimum percentage identity: ", MINIDENT, "\n", sep=""))
cat(paste("Maximum E-value: ", MAXEVALUE, "\n", sep=""))
cat(paste("Number of threads for data.table: ", NCORES_DT, "\n", sep=""))
cat(paste("Output file: ", OUTFILE, "\n", sep=""))
cat("\n")


###############################
############################### Load packages
###############################

cat("Loading main packages...\n")
suppressPackageStartupMessages( library(data.table) ); cat(paste("data.table", packageVersion("data.table"), "\n"))
suppressPackageStartupMessages( library(stringi) ); cat(paste("stringi", packageVersion("stringi"), "\n"))


## Set the number of threads that data.table should use
setDTthreads(threads = NCORES_DT)

set.seed(14789)


###############################
############################### Main workflow
###############################


cat("\nLoading the data...\n")

## Load results of homology search (performed with Diamond)
sb <- fread(INPUTFILE, header = FALSE)

cat(paste("..Number of rows: ", nrow(sb), "\n", sep = ""))


cat("Prepare data...\n")

## Assign column names
colnames(sb) <- c("qseqid", "sseqid", "pident", "length",
	"mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

## Extract read pair info to a separate column
sb[, R := stri_sub(qseqid, -2, -2)]

## Remove read pair info from a query name
sb[, qseqid := stri_sub(qseqid, 1, -4)]

## Subset columns
sb <- sb[, .(qseqid, sseqid, pident, length, evalue, bitscore, R)]


cat("Extracting best hits...\n")

## Extract SeqIDs with two hits
th <- sb[, .(.N), by = .(qseqid, sseqid)]
th <- th[N > 1]

sb2 <- merge(sb, th,
	by.x = c('qseqid', 'sseqid'), by.y = c('qseqid', 'sseqid'),
	all.x = FALSE, all.y = FALSE)

## Calc min e-value and max perc ID per match
sbm <- sb2[, .(pident = max(pident), length = mean(length), evalue = min(evalue), bitscore = max(bitscore)), by = .(qseqid, sseqid)]

## Filter by e-value and identity
sbf <- sbm[pident >= MINIDENT & evalue <= MAXEVALUE ]

## Sort by bitscore, length and identity 
sbf <- sbf[order(qseqid, -bitscore, -length, -pident)]

## Select best hit for R1+R2
sbf <- sbf[, .SD[1], qseqid]

## Clean up
rm(sbm, sb2, th)

## Extract hits that are not captured by two reads, filter and sort
tt <- sb[! qseqid %in% sbf$qseqid ]
tt <- tt[pident >= MINIDENT & evalue <= MAXEVALUE ]
tt <- tt[order(qseqid, R, -bitscore, -length, -pident)]
tt <- tt[, .SD[1], .(qseqid, R)]

## Add read pair info to the header
sbf[, qseqid := paste(qseqid, 12, sep = "_")]
tt[, qseqid := paste(qseqid, R, sep = "_")]

## Prepare data for orthology and functional annotation by EggNog-mapper
## (query - hit - evalue - score)
colz <- c("qseqid", "sseqid", "evalue", "bitscore")

res <- rbind(
  sbf[, ..colz],   # two hits
  tt[, ..colz]     # single hits
)

cat(paste("..Pair-end hits: ", nrow(sbf), "\n", sep = ""))
cat(paste("..Single-end hits: ", nrow(tt), "\n", sep = ""))
cat(paste("..Total number of hits extracted: ", nrow(res), "\n", sep = ""))


cat("Exporting results...\n")

## Export results
fwrite(x = res,
  file = OUTFILE,
  sep = "\t", quote = F, row.names = F, col.names = F)


cat("All done.\n")
cat("\n")
cat("Session info:\n")
sessionInfo()
cat("\n")
