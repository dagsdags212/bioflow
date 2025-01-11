# 
# Perform pairwise and multiple sequence alignment.
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
.PHONY: help init params clean

# Formatting variables
dot := .
comma := ,
empty := 
space := $(empty) $(empty)

# Number of cores
THREADS ?= 8

# Environment manager
ENV_MANAGER ?= micromamba

# Conda environment
ENV := bwf-alignment

# Check if dependencies are installed
dependencies := mafft muscle clustalo jalview

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

# Output format of resulting alignment
OUTFMT ?= fasta

# Filename for output
OUTNAME ?= aln

# MAFFT flags
mafft_opts := --op $(GOP) --ep $(GEP) --thread $(THREADS) --maxiterate $(ITER) --auto --fft

# MUSCLE flags
muscle_opts := -html output/muscle/$(OUTNAME).$(OUTFMT).html

# ClustalO flags
clustal_opts := --threads $(THREADS) --iterations $(ITER) -v -v

# JalView flags
jalview_opts := --colour clustal --type html --headless --close --quit

# Display help message
help:
	@echo
	@echo "alignment.mk: perform sequence alignment"
	@echo
	@echo "Usage:"
	@echo "  make -f src/alignment.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  align      - perform iterative alignment of multiple sequences"
	@echo "  list       - display list of supported sequence aligners"
	@echo "  view       - generate an HTML file for the alignment using JalView"

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 8)"
	@echo "Alignment settings"
	@echo "  FA                path to FASTA file to align"
	@echo "  DIR               path to directory containing FASTA files to align"
	@echo "  ALIGNER           tool of choice for conducting sequence alignment (default: mafft)"
	@echo "  OUTNAME           filename of output alignment (default: aln)"
	@echo "  GOP            	 gap opening penalty (default: 1.53)"
	@echo "  GEP               gap extension penalty (default: 0.0)"
	@echo "  ITER              number of iterations for iterative refinemane (default: 3)"
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-alignment)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"

# Display supported aligners
list:
	@echo "Supported aligners:"
	@echo "  - ClustalO"
	@echo "  - MAFFT"
	@echo "  - MUSCLE"

# Create temporary files to be used for alignment
/tmp/seqs.fa:
ifdef FA
	@cat $(FA) > /tmp/seqs.fa
else ifdef DIR
	@cat $(shell find $(DIR) -type f -name "*.fa*" -o -name "*.fna*") > /tmp/seqs.fa
endif

# Create index for reference genome
align: /tmp/seqs.fa
ifeq ($(ALIGNER),mafft)
	@mkdir -p output/mafft/
	@echo "Aligning sequences with MAFFT"
	mafft $(mafft_opts) $< > output/mafft/$(OUTNAME).$(OUTFMT)
else ifeq ($(ALIGNER),muscle)
	@mkdir -p output/muscle/
	@echo "Aligning sequences with MUSCLE"
	muscle -super5 $< $(muscle_opts) -output output/muscle/$(OUTNAME).$(OUTFMT)
else ifeq ($(ALIGNER),clustalo)
	@mkdir -p output/clustalo
	@echo "Aligning sequences with ClustalO"
	clustalo -i $< $(clustal_opts) -o output/clustalo/$(OUTNAME).$(OUTFMT)
endif

view: $(shell find output/ -maxdepth 2 -name "*$(OUTNAME).$(OUTFMT)")
	@for aln in $^; do \
		echo "Rendering $$aln in HTML format"
		jalview --open $${aln} $(jalview_opts) --image $${aln%.*}.html; \
	done;

clean:
	rm -rf output/mafft/ output/muscle/ output/clustlo/
