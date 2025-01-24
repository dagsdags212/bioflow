#
# Align multiple sequences with MAFFT.
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
FORMATS = fasta clustal

# An optional output format for MSA.
FORMAT ?= fasta

# Alignment method (global, affine, local).
METHOD ?= local

# Clustalo options.
FLAGS := --thread $(THREADS)
ifeq ($(FORMAT),clustal)
	FLAGS += --clustalout
endif

# Conditionally set alignment method.
ifeq ($(METHOD),local)
	FLAGS += --localpair
else ifeq ($(METHOD),affine)
	FLAGS += --genafpair
else ifeq ($(METHOD),global)
	FLAGS += --globalpair
else
	@echo "# Error: unsupported alignment method, falling back to default local."
	FLAGS += --localpair
endif

# A FASTA file for alignment.
FA ?= 

# The output MSA file.
MSA ?= msa/$(basename $(FA)).mafft.aln


help::
	@echo "#"
	@echo "# mafft.mk: align sequences with mafft"
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
MAFFT_CMD := mafft $(FLAGS) $(FA)

$(MSA): $(FA)
	# Create output directory for MSA.
	mkdir -p $(dir $@)

	# Generate MSA file.
	$(ENV_RUN) $(MAFFT_CMD) > $@

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
	@echo "micromamba install mafft"

.PHONY: help run run! clean install
