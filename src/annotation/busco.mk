# 
# Annotate genomes with BUSCO.
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
ENV = bf-annotation

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 2

# Input sequence file or directory.
IN ?= genome.fa

# BUSCO analysis mode.
MODE ?= genome

# Basename for input file.
BASE = $(basename $(notdir $(IN)))

# Output directory for prokka output.
OUTDIR ?= busco

# Output annotation file.
OUT ?= $(OUTDIR)/$(BASE).gff

# E-value threshold.
CUTOFF ?= 0.000001

# Prokka flags.
FLAGS := -c $(THREADS) -e $(CUTOFF)

# Compose prokka command.
BUSCO_CMD = busco $(FLAGS) -i $(IN) -o $(OUTDIR)


help::
	@echo "busco.mk: annotate genomes with BUSCO"
	@echo ""
	@echo "Input:"
	@echo "  IN=$(IN)"
	@echo "  MODE=$(MODE)"
	@echo "  OUTDIR=$(OUTDIR)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform genome annotation of input file"
	@echo "  list       display available BUSCO databases"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all BUSCO output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(BUSCO_CMD)"
	@echo ""

# Input file must exist.
$(IN):
	@if [ ! -f $(IN) ]; then
		@echo "# Error: FASTA file not found (IN=$(IN))."
	fi

$(OUTDIR): $(IN)
	# Reinitialize output directory.
	@if [ -d $(dir $@) ]; then rm -rf $(dir $@); fi
	mkdir -p $(dir $@)

	# Perform gene prediction.
	$(ENV_RUN) $(BUSCO_CMD)

# Invoked gene prediction command.
run: $(OUTDIR)
	ls -lh $^

# Remove busco output and datasets.
run!:: $(OUTDIR)
	rm -rf $^
	rm -rf busco_downloads/

# List prokka databases.
list::
	@$(ENV_RUN) busco --list-datasets

# Alias for run!
clean: run!

# Display installation command.
install::
	@echo "micromamba install busco"

# Rules with no target files.
.PHONY: help run run! list clean install
