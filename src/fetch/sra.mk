#
# Download sequencing reads from the SRA.
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
ENV = bf-fetch

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Store the reads in this directory.
DIR ?= reads

# Sequencing run accession.
SRR ?= SRR1554325

# Number of spots to download.
N ?= ALL

# Name of the unpacked reads.
P1 ?= $(DIR)/$(SRR)_1.fastq
P2 ?= $(DIR)/$(SRR)_2.fastq

# Your custom name the the reads
R1 ?= $(P1)
R2 ?= $(P2)

help:
	@echo "#"
	@echo "# sra.mk: download FASTQ files from SRA"
	@echo "#"
	@echo "# SRR=$(SRR)"
	@echo "# N=$(N) (use N=ALL to download all reads)"
	@echo "#"
	@echo "# R1=$(R1)"
	@echo "# R2=$(R2)"
	@echo "#"
	@echo "# make run|test|aria|clean"
	@echo "#"

ifeq ($(N),ALL)
FLAGS ?= -F --split-files
else
FLAGS ?= -F --split-files -X $(N)
endif

$(R1):
	# Create output directory.
	mkdir -p $(DIR)

	# Download the reads
	$(ENV_RUN) fastq-dump $(FLAGS) -O $(DIR) $(SRR)

	# Rename the first in pair if final name is different.
	if [ "$(P1)" != "$(R1)" ]; then mv -f $(P1) $(R1); fi
	
	# Rename the second in pair if it exists and is different.
	if [ -f "$(P2)" ] && [ "$(P2)" != "$(R2)" ]; then mv -f $(P2) $(R2); fi


# List the data to check if it is paired.
run: $(R1)
	@if [ -f $(R2) ]; then
		@ls -lh $(R1) $(R2)
	else
		@ls -lh $(R1)
	fi

# Download using aria2c.
# This process may be more reliable than fastq-dump.
aria:
	# Extract the ftp links for gzipped fastq files and download.
	$(ENV_RUN) bio search $(SRR) | jq -r '.[].fastq_url[]' | \
		parallel -j 1 --lb make -f $(ROOT_PATH)/aria.mk URL={} DIR=$(DIR) run

# Remove downloaded reads.
clean:
	rm -f $(P1) $(P2) $(R1) $(R2)

# Alternative rule for clean.
run!: clean

# Run the test suite.
test: clean run

# Display dependencies.
install::
	@echo micromamba install sra-tools jq
	@echo pip install bio --upgrade

# Non-file targets.
.PHONY: help run run! test install
