# Path to FASTA file containing sequences to align
FA ?=

# Directory containing FASTA files to align
DIR ?=

# Aligner of choice
ALIGNER ?= mafft

# Scoring metrics for alignment
GOP ?= 1.53
GEP ?= 0.0

# Number of cycles for interative refinement
ITER ?= 1

# Number of permutations
PERMS ?=

# Seed number
SEED ?=

# Set default output to FASTA format
OUTFORMAT = fasta

# Filename for output
OUTFILE ?= aln

# For pairwise alignment
S1 ?=
S2 ?=
