# 
# Perform ab initio gene prediction with Prodigal.
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
THREADS ?= 4

# Supported formats of prodigal output.
FORMATS ?= gbk gff sco

# Format of prodigal output.
FORMAT ?= gbk

# FASTA/Genbank file input.
IN ?= 

# Basename for input file.
BASE = $(basename $(notdir $(IN)))

# Output directory for prodigal output.
DIR ?= prodigal

# Prodigal output file.
OUT ?= $(DIR)/$(BASE)_prodigal_out.$(FORMAT)

# File containing predicted genes from prodigal.
GENES = $(DIR)/$(BASE).genes.out

# File containing sequences of predicted genes.
SEQS = $(DIR)/$(BASE).genes.fa

# Prodigal flags.
FLAGS := -s $(GENES) -d $(SEQS)
FLAGS += $(if $(filter $(FORMATS),$(FORMAT)),-f $(FORMAT))

# Compose prodigal command.
PRODIGAL_CMD = prodigal $(FLAGS) -i $(IN) -o $(OUT)


help::
	@echo "prodigal.mk: ab initio gene prediction with Prodigal"
	@echo ""
	@echo "Input:"
	@echo "  IN=$(IN)"
	@echo "  FORMAT=$(FORMAT)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform ab initio gene prediction in input file"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all Prodigal output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(PRODIGAL_CMD)"
	@echo ""

# Input file must exist.
$(IN):
	@if [ ! -f $(IN) ]; then
		@echo "# Error: FASTA file not found (IN=$(IN))."
	fi

$(OUT): $(IN)
	# Create output directory.
	mkdir -p $(dir $@)

	# Perform gene prediction.
	$(ENV_RUN) $(PRODIGAL_CMD)

# Invoked gene prediction command.
run: $(OUT)
	ls -lh $<

# Remove prodigal output.
run!::
	rm -f $(OUT) $(GENES) $(SEQS)

# Alias for run!
clean: run!

# Display installation command.
install::
	@echo "micromamba install prodigal"

# Rules with no target files.
.PHONY: help run run! clean install
