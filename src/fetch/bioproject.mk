#
# Downloads NCBI run information based on a bioproject number
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Micromamba environment.
ENV = bf-fetch-bioproject

# Run command within environment.
ENV_RUN = micromamba run -n ${ENV}

# Project number
ID ?= PRJNA257197

# Output directory.
OUTDIR ?= metadata

# Project runinfo file.
METADATA ?= ${OUTDIR}/${ID}.csv

# List of accessions from the runinfo file.
ACC_LIST ?= ${OUTDIR}/${ID}_accessions.txt


# General usage information.
help::
	@echo ""
	@echo "bioproject.mk: downloads runinfo for an SRA bioproject"
	@echo ""
	@echo "Usage:"
	@echo "  bf-bioproject [options] <command> ID=<ID>"
	@echo ""
	@echo "Commands:"
	@echo "  run            download SRA bioproject metadata"
	@echo "  get            an alias for 'run'"
	@echo "  accessions     get a list of SRA accessions"
	@echo "  install        initialize conda environment"
	@echo "  clean          remove all output generated by this program"
	@echo ""
	@echo "Options:"
	@echo "  ID             a bioproject accession (required)"
	@echo "  OUTDIR         a directory path for storing output"
	@echo ""

# Print out example usage.
example:
	@echo ""
	@echo "# Retrieve the runinfo associated with a bioproject accession"
	@echo "make -f bioproject.mk ID=PRJNA257197 run"
	@echo ""
	@echo "# Get a list of read accessions from the same project"
	@echo "make -f bioproject.mk ID=PRJNA257197 accessions"
	@echo ""

# Project run information.
${METADATA}:
	# Create output directory.
	mkdir -p $(dir $@)

	# Retrieve run metadata.
	${ENV_RUN} bio search ${ID} --header --csv > $@

# Extract accessions from project runinfo.
${ACC_LIST}: ${METADATA}
	@cat $< | cut -d, -f1 | tail +2 > $@

# Target to download all the data.
run:: ${METADATA}
	@ls -lh ${METADATA}

# Remove bioproject
run!::
	rm -rf $(dir ${METADATA})

# For backward compatibility.
get: run
clean: run!

accessions: ${ACC_LIST}

# Run test suite.
test: clean run

# Installation instructions
install::
	micromamba create -n ${ENV} --yes
	${ENV_RUN} micromamba install pip --yes
	${ENV_RUN} pip install bio --upgrade

.PHONY: help example accessions run run! get clean test install
