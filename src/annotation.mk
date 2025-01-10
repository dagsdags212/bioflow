# 
# Perform gene prediction to annotate genomes.
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
.PHONY: help init params clean datasets

# Formatting variables
dot := .
comma := ,

# Number of cores
THREADS ?= 8

# Environment manager
ENV_MANAGER := micromamba

# Conda environment
ENV := bwf-annotation

# Check if dependencies are installed
dependencies := busco

# Sequence to annotate
FA ?=

# Output format for Prodigal
PRODIGAL_OUTFMT ?= gff

# Prodigal flags
prodigal_opts := -o output/prodigal/predicted_genes.$(PRODIGAL_OUTFMT) 
prodigal_opts += -a output/prodigal/predicted_proteins.faa -f $(PRODIGAL_OUTFMT)

# BUSCO flags
busco_opts := -c $(THREADS) -f

# BUSCO mode
MODE ?= genome

# Only specify BUSCO mode if parameters is valid (genome, transcriptome)
ifeq ($(MODE),genome)
busco_opts += -m genome
endif

ifeq ($(MODE),transcriptome)
busco_opts += -m transcriptome
endif

# Specify domain of source FASTA file
DOMAIN ?=

# Determine domain of input sequence
ifneq (,$(findstring euk,$(DOMAIN)))
busco_opts += --auto-lineage-euk
endif

ifneq (,$(findstring prok,$(DOMAIN)))
busco_opts += --auto-lineage-prok
endif

# Search all lineages if domain is not defined
ifndef DOMAIN
busco_opts += --auto-lineage
endif

# Display help message
help:
	@echo
	@echo "annotation.mk: predict genes and provide annotation for genomes"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/annotation.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  annotate   - run annotation pipeline on sequence file"
	@echo "  datasets   - list available BUSCO datasets"
	@echo "  predict    - run HMM-based gene prediction using Prodigal"

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo "Annotation settings"
	@echo "  FA                path to target FASTA file"
	@echo "  MODE              BUSCO analysis mode to run [genome|protein|transcriptome] (default: genome)"
	@echo "  DOMAIN            specify the domain where the input sequence belongs to [prokaryote|eukaryote]"
	@echo "  PRODIGAL_OUTFMT   specify output format for Prodigal [gff|gbk|sco] (default: gff)"
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-mapping)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"

# List available BUSCO datasets
datasets:
	@busco --list-datasets

# Run ab initio gene prediction using prodigal
predict:
	mkdir -p output/prodigal/
	prodigal -i $(FA) $(prodigal_opts)

# Run BUSCO annotation pipeline
annotate:
	mkdir -p output/busco/
	busco -i $(FA) -o output/busco/ $(busco_opts)

clean:
	rm -rf output/busco/ output/prodigal/
	rm -rf busco_downloads/
