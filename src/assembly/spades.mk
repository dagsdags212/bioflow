#
# Assemble reads into contigs using megahit.
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

# Set memory allocation to 50% of machine's capacity.
MEM ?= 8

# An SRA accession.
ACC ?=

# First of read pair.
R1 ?=

# Second of read pair. If not specified, perform assembly in SE mode.
R2 ?=

# Output directory.
OUTDIR ?= spades

# Target output.
OUT = $(OUTDIR)/scaffolds.fasta

# Sequencing platform.
PLATFORM ?=

# K-mer size(s).
K ?= auto

# Megahit options.
FLAGS := -t $(THREADS) -m $(MEM) -k $(K) --gfa11


# Compose megahit command.
ifeq ($(R2),)
	SPADES_CMD = spades.py $(FLAGS) -s $(R1) -o $(OUTDIR)
else
	SPADES_CMD = spades.py $(FLAGS) -1 $(R1) -2 $(R2) -o $(OUTDIR)
endif


help::
	@echo "spades.mk: assembled reads into contigs using spades"
	@echo ""
	@echo "Input:"
	@echo "  R1=$(R1)"
	@echo "  R2=$(R2)"
	@echo "  OUTDIR=$(OUTDIR)"
	@echo ""
	@echo "  THREADS=$(THREADS)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform sequence assembly"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all spades output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(SPADES_CMD)"
	@echo ""

# R1 must exist.
$(R1):
	@if [ ! -f $(R1) ]; then
		echo "Error: R1 not found (R1=$(R1))"
	fi
	@exit -1

# If set, R2 must exist.
ifneq ((R2),)
$(R2):
	@if [ ! -f $(R2) ]; then
		echo "Error: R2 not found (R2=$(R2))"
	fi
	@exit -1
endif

$(OUT): $(R1) $(if $(R2),$(R2))
	# Remove output directory if it already exists.
	if [ -d $(dir $@) ]; then
		rm -rf $(dir $@)
	fi

	# Perform sequence assembly with megahit.
	$(ENV_RUN) $(SPADES_CMD)

# Invoke spades command.
run: $(OUT)
	ls -lh $(dir $^)
	
run!::
	rm -rf $(OUTDIR)

# Alternative target to run!
clean: run!

# Print command for installing dependencies.
install::
	@echo "micromamba install spades"

.PHONY: help run run! clean install
