# ===== GLOBAL options =====
# Program to perform assembly.
ASSEMBLER ?= megahit

# Path to single-end FASTQ file.
SE ?=

# Path to pair-end #1 FASTQ file.
PE1 ?=

# Path to pair-end #2 FASTQ file.
PE2 ?=

# Path to directory for storing output.
OUTDIR ?=

# Number of CPU cores.
THREADS ?=

# K-mer length.
KMER_SIZE ?= 31

# ===== MEGAHIT options =====


# ===== SPADES options =====
# Minimum contig length to filter.
MIN_CONTIG_LEN ?=

# Sequencing platform.
PLATFORM ?=

# ===== MINIA options =====
#  Prefix to include in output files.
PREFIX ?=



