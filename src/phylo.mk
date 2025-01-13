# 
# Conduct phylogenetic tree inference from an alignment
#

# import config variables
include src/_config.mk

# import global variables
include src/_globals.mk

.PHONY: help params init clean

# Project root
ROOT_DIR = $(shell dirname $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST)))))

# Conda environment
ENV := bf-phylo

# Path to conda environment
ENV_DIR = $(shell $(ENV_MANAGER) info | grep "envs directories" | cut -d ":" -f 2 | xargs)/$(ENV)

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
N ?= 10

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
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo
	@echo "Sequence alignment"
	@echo "  MSA               path to alignment file for generating trees"
	@echo "  MODEL             specify evolutionary mode to use (default: GTRGAMMA)"
	@echo "  N                 number of starting trees (ML) or replicates (bootstrap)"
	@echo "  PREFIX            add a prefix to the output file"
	@echo "  RAPID_BOOTSTRAP   perform rapid bootstrapping [true|false](default: false)"
	@echo "  SEED              provide a seed for ML/bootstrapping (default: 12345)"
	@echo "  TOOL              specify tools for performing tree inference (default: raxml)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 8)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-phylo)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Display evolutionaryu models
models:
	@echo
	@echo "Nucleotide substitutions models:"
	@echo "  JC     - Jukes and Cantor (1969)"
	@echo "  K80    - Kimura (1980)"
	@echo "  K81    - Kimura [equal freq.](1981)"
	@echo "  K81uf  - Kimura [unequal freq.](1981)"
	@echo "  F81    - Felsenstein (1981)"
	@echo "  HKY    - Hasegawa et al. (1985)"
	@echo "  TN93   - Tamura and Nei [unequal freq.](1993)"
	@echo "  TN93ef - Tamura and Nei [equal freq.](1993)"
	@echo "  TVM    - Posada [unequal freq.](2003)"
	@echo "  TVMef  - Posada [equal freq.](2003)"
	@echo "  GTR    - Tavare (1986)"
	@echo
	@echo "Rate heterogeneity models:"
	@echo "  G      - discrete GAMMA with 4 categories"
	@echo "  GA     - discrete GAMMA with median categories"
	@echo "  Gn     - discrete GAMMA with n categories"
	@echo "  Rn     - free rate with n categories"
	@echo "  Rn     - free rate with n categories"
	@echo

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

clean:
	rm -rf output/raxml output/iqtree
