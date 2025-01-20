# === GLOBAL options ===
# Sequence to annotate
FA ?=

# === PRODIGAL options ===
# Specify output format for Prodigal
PRODIGAL_OUTDIR ?= output/prodigal

# Specify output directory for Prodigal
PRODIGAL_OUTFMT ?= gff

# === BUSCO options ===
# BUSCO mode
MODE ?= genome

# Specify domain of source FASTA file
DOMAIN ?= prokaryote
