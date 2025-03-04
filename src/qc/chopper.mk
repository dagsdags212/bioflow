#
# Perform read and quality filtering for long-reads using chopper.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Micromamba environment.
ENV = bf-qc-chopper

# Run command within environment.
ENV_RUN = micromamba run -n ${ENV}

# Read file to process.
FQ ?=

# Output directory to store filtered reads.
OUT ?= $(if ${FQ}, $(shell dirname ${FQ}))

# Filtered FASTQ file.
FILTERED = ${OUT}/filtered_$(notdir ${FQ})

# Minimum read length to keep.
MINLEN ?= 500

# Minimum read quality to keep.
MINQUAL ?= 10

# Number of worker threads.
THREADS ?= 4


# Display usage.
help::
	@echo ""
	@echo "chopper.mk: filter long-reads by length and quality"
	@echo ""
	@echo "Usage:"
	@echo "  bf-chopper FQ=<FQ> [options]"
	@echo ""
	@echo "Commands:"
	@echo "  run            filter long-reads by length and quality"
	@echo "  filter         an alias for 'run'"
	@echo "  install        initialize conda environment"
	@echo "  clean          remove all output generated by this program"
	@echo ""
	@echo "Options:"
	@echo "  FQ             a directory path containing the target reads"
	@echo "  MINLEN         minimum read length (default: 500)"
	@echo "  MINQUAL        minimum read quality (default: 10)"
	@echo "  OUT            a directory path for storing reports (default: .)"
	@echo "  THREADS        number of worker threads to use (default: 4)"
	@echo ""

# chopper options.
chopper_opts ?= -q ${MINQUAL} -l ${MINLEN} --threads ${THREADS}

# Input FASTQ file must exist.
${FQ}:
	echo "Error: FASTQ file not found (FQ=${FQ})"
	exit -1

# Run chopper.
${FILTERED}: ${FQ}
	# Create directory for filtered reads.
	mkdir -p $(dir $@)

	# Filter reads.
	${ENV_RUN} chopper ${chopper_opts} -i ${FQ} | gzip > $@.gz

# Invoke chopper.
run: ${FILTERED}
	@ls -lh ${FILTERED}

# An alias for 'run'.
filter: run

# Remove filtered reads.
run!:
	rm -rf ${FILTERED}

# An alias for 'run!'.
clean: run!

# Run test suite.
test::
	make -f $${BIOFLOW}/src/fetch/sra.mk SRR=SRR8848097 N=1000 MODE=SE run
	make -f $${BIOFLOW}/src/qc/chopper.mk FQ=reads/SRR8848097.fastq run! run

DEPS := chopper
# Command for installing dependencies.
install::
	micromamba create -n ${ENV}
	${ENV_RUN} micromamba install ${DEPS} --yes

# Rules without target files.
.PHONY: help run run! clean test install
