#
# Perform phylogenetic tree inference with raxml-ng.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Absolute path of parent directory.
ROOT_PATH = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Micromamba environment.
ENV = bf-phylo

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 4

# Method for tree construction.
METHOD ?= ml

# Number of ML searches.
N ?= 10

# Set initial seed.
SEED ?=

# Nucleotide substitution model
MODEL ?= GTR+G

# Input alignment.
ALN ?= 

# Basename of alignment file.
BASENAME = $(basename ALN)

# Output directory.
OUTDIR ?= raxml

# Target output.
OUT ?= $(OUTDIR)/$(BASENAME).bestTree

# Raxml-ng options.
FLAGS := --threads $(THREADS) --model $(MODEL) 
FLAGS += $(if $(SEED),--seed $(SEED))
FLAGS += $(if $(SEED),--seed $(SEED))

ifeq ($(METHOD),ml)
	FLAGS += --tree pars{$(N)}
else ifeq ($(METHOD),bootstrap)
	FLAGS += --boostrap --bs-trees $(N)
else
	@echo "Error: invalid tree construction method (METHOD=$(METHOD))"
	@echo "Set value to METHOD=[ml|bootstrap]"
	@exit -1
endif

# Compose raxml-ng command.
RAXML_CMD = raxml-ng $(FLAGS) --msa $(ALN)


help::
	@echo "raxml-ng.mk: perform tree inference with raxml-ng"
	@echo ""
	@echo "Input:"
	@echo "  ALN=$(ALN)"
	@echo "  OUTDIR=$(OUTDIR)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform phylogenetic inference and construct a tree"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all raxml-ng output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(RAXML_CMD)"
	@echo ""

# Alignment file must exist.
$(ALN):
	@if [ ! -f $(ALN) ]: then
		echo "Error: alignment file does not exist (ALN=$(ALN))"
	fi
	@exit -1

# Generate phylogenetic tree.
$(OUT): $(ALN)
	# Create output directory.
	mkdir -p $(dir $@)

	# Peform tree inference
	$(ENV_RUN) $(RAXML_CMD)

	# Move output files to target directory.
	mv $(ALN).raxml.* $(dir $@)

# Invoked tree inference command.
run: $(OUT)
	ls -lh $(dir $^)

# Remove raxml-ng output.
run!::
	rm -rf $(dir $(OUT))

# Alias for run!
clean: run!

# Print command for installing dependencies.
install::
	@echo "micromamba install raxml-ng"

.PHONY: help clean install
