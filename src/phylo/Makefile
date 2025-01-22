# 
# Conduct phylogenetic tree inference from an alignment
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# Absoulate path to conda environment files
ENV_ROOT := ${BIOFLOW_PREFIX}/envs

# Conda environment for module
ENV := bf-phylo

# Execute within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/phylo.mk

# Display help message.
help:
	@echo
	@echo "phylo.mk: infer phylogenies and construct trees from sequence alignments"
	@echo
	@echo "Usage:"
	@echo "  bf-phylo <command> [options]"
	@echo
	@echo "Commands:"
	@echo "  bootstrap - perform a simple bootstrap analysis"
	@echo "  ML        - perform a maximum likelihood search"
	@echo "  models    - list available evolutionary models"
	@echo "  init      - download all dependencies"
	@echo "  params    - display complete list of parameters"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -f $(ENV_ROOT)/bf-phylo.yml

# Display available parameters
params:
	@echo
	@echo "Sequence alignment"
	@echo "  TOOL              specify tools for performing tree inference (default: raxml)"
	@echo "  MSA               path to alignment file for generating trees"
	@echo "  MODEL             specify evolutionary mode to use (default: GTR+G)"
	@echo "  N                 number of starting trees (ML) or replicates (bootstrap)"
	@echo "  PREFIX            add a prefix to the output file"
	@echo "  BOOTSTRAPS        number of bootstraps to perform (default: 1001)"
	@echo "  SEED              provide a seed for ML/bootstrapping (default: 12345)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 8)"
	@echo "  OUTDIR            path to output directory"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-phylo)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Display evolutionary models.
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

# RAXML options.
raxml_opts := MSA=$(MSA)
raxml_opts += $(if $(N),N=$(N))
raxml_opts += $(if $(MODEL),MODEL=$(MODEL))
raxml_opts += $(if $(SEED),SEED=$(SEED))
raxml_opts += $(if $(THREADS),THREADS=$(THREADS))

# IQTREE options.
iqtree_opts := MSA=$(MSA)
iqtree_opts += $(if $(SEED),SEED=$(SEED))
iqtree_opts += $(if $(MODEL),MODEL=$(MODEL))
iqtree_opts += $(if $(PREFIX),PREFIX=$(PREFIX))
iqtree_opts += $(if $(OUTDIR),OUTDIR=$(OUTDIR))

# Perform phylogenetic inference using tool of preference.
ML: $(MSA)
ifeq ($(TOOL),raxml)
	@make -f $(MK_ROOT)/run_raxml.mk $(raxml_ml_opts) METHOD=ml
else ifeq ($(TOOL),iqtree)
	@make -f $(MK_ROOT)/run_iqtree.mk $(iqtree_opts) METHOD=ml
else
	echo "Error: specified tool is not currently supported"
endif

# Peform bootstrapping.
bootstrap: $(MSA)
ifeq ($(TOOL),raxml)
	@make -f $(MK_ROOT)/run_raxml.mk $(raxml_ml_opts) METHOD=bootstrap
else ifeq ($(TOOL),iqtree)
	@make -f $(MK_ROOT)/run_iqtree.mk $(iqtree_opts) METHOD=bootstrap \
		$(if $(BOOTSTRAPS),BOOTSTRAPS=$(BOOTSTRAPS))
else
	echo "Error: specified tool is not currently supported"
endif

# Delete all generated output.
clean:
	rm -rf $(OUTDIR)

.PHONY: help params init clean
