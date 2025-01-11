# 
# Conduct phylogenetic tree inference from an alignment 
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
ENV := bwf-phylo

# Check if dependencies are installed
dependencies := raxml iqtree blast bedtools

# Path to alignment file
ALN ?=

# Prefered tools for running phylogenetic inferences
TOOL ?= raxml

# Output filename
OUTNAME ?=

# Model of evolution
MODEL ?= GTRGAMMA

# Specify seed
SEED ?=

# Number of ML searches to perform
N ?=

# Perform rapid bootstrappint
RAPID_BOOTSTRAP ?= true

# Global RAxML flags
raxml_global_opts := $(if $(OUTNAME),-n output/raxml/$(OUTNAME),-n output/raxml/result_)
raxml_global_opts += $(if $(MODEL),-m $(MODEL))
raxml_global_opts += $(if $(SEED),-p $(SEED),12345)
raxml_global_opts += $(if $(N),-# $(N),20)

# RAxML flags for ML search
raxml_ml_opts := $(raxml_global_opts)

# RAxML flags for bootstraping
raxml_bootstrap_opts := $(raxml_global_opts)

ifeq ($(RAPID_BOOTSTRAP),true)
	raxml_bootstrap_opts += -x $(SEED)
else
	raxml_bootstrap_opts += -b $(SEED)
endif

# IQTree flags
iqtree_opts :=

# Display help message
help:
	@echo
	@echo "phylo.mk: infer and generate phylogenies from sequence alignments"
	@echo
	@echo "Usage:"
	@echo "  make -f src/phylo.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  bootstrap - perform a simple bootstrap analysis"
	@echo "  ML        - perform a maximum likelihood search"
	@echo "  models    - list available evolutionary models"
	@echo "  init      - download all dependencies"
	@echo "  params    - display complete list of parameters"

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 8)"
	@echo "Alignment settings"
	@echo "  ALN               path to input alignment file for generating trees"
	@echo "  MODEL             specify evolutionary mode to use (default: GTRGAMMA)"
	@echo "  N                 number of bootstrap searches to perform"
	@echo "  OUTNAME           filename for output tree"
	@echo "  RAPID_BOOTSTRAP   perform rapid bootstrapping [true|false](default: false)"
	@echo "  SEED              provide a seed for ML/bootstrapping (default: 12345)"
	@echo "  TOOL              specify tools for performing tree inference (default: raxml)"
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-phylo)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"

# Display evolutionaryu models
models:
	@echo "Evolutionary models:"
	@echo "  - CAT"
	@echo "  - GAMMA"

# Perform phylogenetic inference using tool of preference
ML: $(ALN)
ifeq ($(TOOL),raxml)
	@mkdir -p output/raxml
	raxmlHPC -s $< $(raxml_ml_opts)
endif

# Peform bootstrapping
bootstrap: $(ALN)
ifeq ($(TOOL),raxml)
	@mkdir -p output/raxml
	raxmlHPC -s $< $(raxml_bootstrap_opts)
endif
