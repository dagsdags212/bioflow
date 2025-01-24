#
# Align multiple sequences with clustal-omega.
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
ENV = bf-alignment

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 4

# List of supported output formats.
FORMATS = fasta clustal msf phylip selex stockholm vienna

# An optional output format for MSA.
FORMAT ?= fasta

# Clustalo options.
FLAGS := --threads $(THREADS) $(filter $(FORMATS),--outfmt=$(FORMAT))

# A FASTA file for alignment.
FA ?= 

# The output MSA file.
MSA ?= msa/$(basename $(FA)).clustalo.aln


help::
	@echo "#"
	@echo "# clustalo.mk: align sequences with clustalo"
	@echo "#"
	@echo "# FA=$(FA)"
	@echo "#"
	@echo "# MSA=$(MSA)"
	@echo "#"
	@echo "# make run|align|install|clean"
	@echo "#"

# Fasta file must exist.
$(FA):
	@if [ ! -f $(FA) ]; then
		@echo "# Error: FASTA file not found (FA=$(FA))."
	fi

# Compose clustalo command.
CLUSTALO_CMD := clustalo $(FLAGS) -i $(FA)

$(MSA): $(FA)
	# Create output directory for MSA.
	mkdir -p $(dir $@)

	# Generate MSA file.
	$(ENV_RUN) $(CLUSTALO_CMD) > $@

# Invoke MSA generation.
align: $(MSA)

# Alternative target for align.
run: align

# Remove MSA file.
run!::
	rm -f $(MSA)

# Alternative target for run!
clean: run!

# Display dependencies.
install::
	@echo "micromamba install clustalo"

.PHONY: help run run! clean install

