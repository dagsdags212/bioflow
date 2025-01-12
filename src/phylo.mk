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
MSA ?=

# Prefered tools for running phylogenetic inferences
TOOL ?= raxml

# Output filename
PREFIX ?=

# Model of evolution
MODEL ?= GTR+G

# Specify seed
SEED ?= 12345

# Number of ML searches to perform
N ?=

# Perform rapid bootstrappint
RAPID_BOOTSTRAP ?= true

# Global RAxML flags
raxml_global_opts := --threads $(THREADS) --model $(MODEL) --seed $(SEED)
raxml_global_opts += output/raxml/$(if $(PREFIX),$(PREFIX),result)

# RAxML flags for ML search
raxml_ml_opts := $(raxml_global_opts)
raxml_ml_opts += $(if $(N),--tree pars{$(N)})

# RAxML flags for bootstraping
raxml_bootstrap_opts := $(raxml_global_opts) --bootstrap
raxml_bootstrap_opts += $(if $(N),--bs-trees $(N))

# TODO: add support for PHYLIP

# Phylip flags
phylip_opts := 

# TODO: add support for IQTREE

# IQTree flags
iqtree_opts := --alisim output/iqtree/$(if $(PREFIX),$(PREFIX),result) -alninfo
iqtree_opts += -af $(OUTFMT) -T $(THREADS) -m $(MODEL) --seed $(SEED)

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
	@echo ""
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 8)"
	@echo ""
	@echo "Alignment settings"
	@echo "  MSA               path to alignment file for generating trees"
	@echo "  MODEL             specify evolutionary mode to use (default: GTRGAMMA)"
	@echo "  N                 number of starting trees (ML) or replicates (bootstrap)"
	@echo "  PREFIX            add a prefix to the output file"
	@echo "  RAPID_BOOTSTRAP   perform rapid bootstrapping [true|false](default: false)"
	@echo "  SEED              provide a seed for ML/bootstrapping (default: 12345)"
	@echo "  TOOL              specify tools for performing tree inference (default: raxml)"
	@echo ""
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-phylo)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo ""

# Display evolutionaryu models
models:
	@echo "Evolutionary models:"
	@echo "  - CAT"
	@echo "  - GAMMA"

# Perform phylogenetic inference using tool of preference
ML: $(MSA)
ifeq ($(TOOL),raxml)
	@mkdir -p output/raxml
	raxml-ng --msa $< $(raxml_ml_opts)
else ifeq ($(TOOL),iqtree)
	@mkdir -p output/iqtree
	iqtree -s $< $(iqtree_opts)
endif

# Peform bootstrapping
bootstrap: $(MSA)
ifeq ($(TOOL),raxml)
	@mkdir -p output/raxml
	raxml-ng --msa $< $(raxml_bootstrap_opts)
endif
