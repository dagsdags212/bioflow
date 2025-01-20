# 
# Perform pairwise and multiple sequence alignment.
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# Absoulate path to conda environment files
ENV_ROOT := ${BIOFLOW_PREFIX}/envs

# Conda environment for module
ENV := bf-alignment

# Execute within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/alignment.mk

.PHONY: help params init clean

# Display help message
help:
	@echo
	@echo "alignment.mk: perform sequence alignment"
	@echo
	@echo "Usage:"
	@echo "  bfalign <command> [options]"
	@echo
	@echo "Commands:"
	@echo "  pairwise   - perform alignment on two given sequences"
	@echo "  msa        - perform iterative alignment of multiple sequences"
	@echo "  list       - display list of supported sequence aligners"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo
	@echo "Alignment settings"
	@echo "  FA                path to FASTA file to align"
	@echo "  DIR               path to directory containing FASTA files to align"
	@echo "  ALIGNER           tool of choice for conducting sequence alignment (default: mafft)"
	@echo "  OUTNAME           filename of output alignment (default: aln)"
	@echo "  GOP            	 gap opening penalty (default: 1.53)"
	@echo "  GEP               gap extension penalty (default: 0.0)"
	@echo "  ITER              number of iterations for iterative refinemane (default: 3)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 8)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-alignment)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Display supported aligners
list:
	@echo
	@echo "Supported aligners:"
	@echo "  - ClustalO"
	@echo "  - MAFFT"
	@echo "  - MUSCLE"
	@echo

# MAFFT options
mafft_opts := $(if $(THREADS), -t $(THREADS))
mafft_opts := $(if $(OP), -o $(OP))
mafft_opts := $(if $(EP), -e $(EP))

# MUSCLE options
muscle_opts := $(if $(PERMS), -p $(PERMS))
muscle_opts += $(if $(SEED), -s $(SEED))

# Clustal-Omega options
clustalo_opts := $(if $(PERMS), -p $(PERMS)) -v -v

# Perform pairwise alignment
pairwise:
ifndef S1
	@echo "Error: reference sequence (S1) not provided"
	exit 1
else ifndef S2
	@echo "Error: query sequence (S2) not provided"
	exit 1
endif
	@bio align $(S1) $(S2)

# Perform multiple sequence alignment
msa:
# Align with MAFFT
ifeq ($(ALIGNER),mafft)
	@mkdir -p output/mafft/
	@. $(MK_ROOT)/run_mafft.sh $(mafft_opts) $(if $(FA),$(FA),$(DIR))
# Align with MUSCLE
else ifeq ($(ALIGNER),muscle)
	@mkdir -p output/muscle/
	@. $(MK_ROOT)/run_muscle.sh $(muscle_opts) $(if $(FA),$(FA),$(DIR))
# Align with ClustalO
else ifeq ($(ALIGNER),clustalo)
	@mkdir -p output/clustalo
	@. $(MK_ROOT)/run_clustalo.sh $(clustalo_opts) $(if $(FA),$(FA),$(DIR))
else
	echo "Error: specified aligner not currently supported"
endif

# Remove all generated output
clean:
	rm -rf output/mafft/ output/muscle/ output/clustlo/
