# 
# Perform gene prediction to annotate genomes.
#

# Absoluate path for the `bioflow/src` directory
SRC_ROOT := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

# Absolute path for `bioflow` directory
PROJECT_ROOT := $(shell dirname $(SRC_ROOT))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := $(PROJECT_ROOT)/config

# import config variables
include $(CONFIG_ROOT)/_config.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

.PHONY: help params init clean

# Project root
ROOT_DIR = $(shell dirname $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST)))))

# Conda environment
ENV := bf-annotation
# Path to conda environment
ENV_DIR = $(shell $(ENV_MANAGER) info | grep "envs directories" | cut -d ":" -f 2 | xargs)/$(ENV)

# Check if dependencies are installed
dependencies := prokka busco artemis

# Sequence to annotate
FA ?=

# Output format for Prodigal
PRODIGAL_OUTFMT ?= gff

# Prodigal flags
prodigal_opts := -o output/prodigal/predicted.$(PRODIGAL_OUTFMT) 
prodigal_opts += -a output/prodigal/predicted.faa -f $(PRODIGAL_OUTFMT)
prodigal_opts += -s output/prodigal/predicted.txt

# BUSCO flags
busco_opts := -c $(THREADS) -f

# BUSCO mode
MODE ?= genome

# Only specify BUSCO mode if parameters is valid (genome, transcriptome)
ifeq ($(MODE),genome)
busco_opts += -m genome
else ifeq ($(MODE),transcriptome)
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

# Prokka flags
prokka_opts := --gffver 3 --force

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
	@echo "  view       - visualize the annotated genome using Artemis"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo
	@echo "Annotation settings"
	@echo "  FA                path to target FASTA file"
	@echo "  MODE              BUSCO analysis mode to run [genome|protein|transcriptome] (default: genome)"
	@echo "  DOMAIN            specify the domain where the input sequence belongs to [prokaryote|eukaryote]"
	@echo "  PRODIGAL_OUTFMT   specify output format for Prodigal [gff|gbk|sco] (default: gff)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-annotation)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# List available BUSCO datasets
datasets:
	@busco --list-datasets

# Run ab initio gene prediction using prodigal
predict:
	@mkdir -p output/prodigal/
	prodigal -i $(FA) $(prodigal_opts)

# Run BUSCO annotation pipeline
annotate:
	@# annotate using BUSCO
	@mkdir -p output/busco/
	@echo "Annotating $(FA) with BUSCO"
	busco -i $(FA) -o output/busco/ $(busco_opts)

	@# annotate using Prokka
	@mkdir -p output/prokka/
	@echo "Annotating $(FA) with Prokka"
	prokka $(prokka_opts) --outdir output/prokka/ $(FA)

view: $(wildcard output/prokka/*.gff) $(FA)
	art $(FA) +$<

clean:
	rm -rf output/busco/ output/prodigal/ output/prokka/
	rm -rf busco_downloads/
