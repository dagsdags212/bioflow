#
# Perform phylogenetic tree inference with iqtree.
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

# Set initial seed.
SEED ?=

# Nucleotide substitution model
MODEL ?= GTR+G

# Input alignment.
ALN ?=

# Basename of alignment file.
BASENAME = $(basename $(ALN))

# Output directory.
OUTDIR ?= iqtree

# Target output.
OUT ?= $(OUTDIR)/$(BASENAME).iqtree

# IQTREE options.
FLAGS := -nt $(THREADS) -m $(MODEL) -pre $(BASENAME)

# Compose iqtree command.
IQTREE_CMD = iqtree $(FLAGS) -s $(ALN)


help::
	@echo "iqtree.mk: perform tree inference with iqtree"
	@echo ""
	@echo "Input:"
	@echo "  ALN=$(ALN)"
	@echo "  OUTDIR=$(OUTDIR)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform phylogenetic inference and construct trees"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all IQTREE output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(IQTREE_CMD)"
	@echo ""

# Alignment file must exist.
$(ALN):
	@if [ ! -f $(ALN) ]; then
		echo "Error: alignment file does not exist (ALN=$(ALN))"
	fi
	@exit -1

OUT_SFX = iqtree bionj ckp.gz log mldist treefile
OUT_FILES = $(addprefix $(BASENAME).,$(OUT_SFX))

# Generate phylogenetic tree.
$(OUT): $(ALN)
	# Create output directory.
	mkdir -p $(dir $@)

	# Peform tree inference
	$(ENV_RUN) $(IQTREE_CMD)

	# Move output to target directory.
	mv $(BASENAME).* $(OUTDIR)

# Invoked tree inference command.
run: $(OUT)
	ls -lh $(dir $^)

# Remove iqtree output.
run!:: $(OUT)
	rm -rf $(dir $^)

# Alias for run!
clean: run!

# Print command for installing dependencies.
install::
	@echo "micromamba install iqtree"

.PHONY: help clean install
