#
# Assemble reads into contigs using minia.
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
ENV = bf-assembly

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 4

# Minia only accepts a single FASTQ file.
R1 ?=

# Basename of FASTQ file.
BASENAME = $(basename $(R1))

# Output directory.
OUTDIR ?= minia

# Target output.
OUT = $(OUTDIR)/contigs.fa

# K-mer size(s).
K ?= 31

# Minia options.
FLAGS := -nb-cores $(THREADS) -kmer-size $(K)

# Compose minia command.
MINIA_CMD = minia $(FLAGS) -in $(R1) -out $(OUTDIR)/$(BASENAME) -out-dir $(OUTDIR)


help::
	@echo "minia.mk: assembled reads into contigs using minia"
	@echo ""
	@echo "Input:"
	@echo "  R1=$(R1)"
	@echo "  OUTDIR=$(OUTDIR)"
	@echo ""
	@echo "  THREADS=$(THREADS)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform sequence assembly"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all minia output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(MINIA_CMD)"
	@echo ""

# R1 must exist.
$(R1):
	@if [ ! -f $(R1) ]; then
		echo "Error: R1 not found (R1=$(R1))"
	fi
	@exit -1

$(OUT): $(R1)
	# Remove output directory if it already exists.
	if [ -d $(dir $@) ]; then
		rm -rf $(dir $@)
	fi

	# Perform sequence assembly with minia.
	$(ENV_RUN) $(MINIA_CMD)

	# 

# Invoke minia command.
run: $(OUT)
	ls -lh $(dir $^)
	
run!::
	rm -rf $(OUTDIR)

# Alternative target to run!
clean: run!

# Print command for installing dependencies.
install::
	@echo "micromamba install minia"

.PHONY: help run run! clean install
