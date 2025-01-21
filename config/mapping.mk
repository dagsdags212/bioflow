# Tool to perform read mapping.
MAPPER ?= bwa

# Path to reference file (FASTA).
REF ?= ref.fa

# Path to single-end FASTQ file.
SE ?= 

# Path to pair-end #1 FASTQ file.
PE1 ?= reads1.fq

# Path to pair-end #2 FASTQ file.
PE2 ?= reads2.fq

# Path to directory for storing output.
OUTDIR ?=

# Filename for output alignment.
PREFIX ?=

# Number of CPU cores.
THREADS ?= 4

# Score threshold for filtering mapped reads.
MINQUAL ?= 20
