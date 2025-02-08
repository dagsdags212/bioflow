#
# Perform quality control and generate summary reports for FASTQ files.
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
ENV = bf-qc-fastqc

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Directory containing reads to include in QC.
DIR ?= reads

# Contaminant file.
CONTAMINANT_FILE ?=

# Output directory to store report.
OUT ?= fastqc

# Number of worker threads.
THREADS ?= 2

# Display usage.
help::
	@echo ""
	@echo "fastqc.mk: generate read quality reports with FASTQC"
	@echo ""
	@echo "Usage:"
	@echo "  make -f fastqc.mk [options] DIR=<dir>"
	@echo ""
	@echo "Commands:"
	@echo "  run                 generate FASTQC reports on target sequence(s)"
	@echo "  install             initialize conda environment"
	@echo "  clean               remove all output generated by this program"
	@echo ""
	@echo "Options:"
	@echo "  DIR                 a directory path containing the target reads"
	@echo "  OUT                 a directory path for storing reports"
	@echo "  CONTAMINANT_FILE    a file containing contaminant sequences"
	@echo "  THREADS             number of worker threads to use"
	@echo ""

# Fastqc options.
fastqc_opts ?= $(if ${CONTAMINANT_FILE}, -c ${CONTAMINANT_FILE}) \
							 -t ${THREADS} -o ${OUT}

# Directory storing reads must exist.
${DIR}:
	@if [ ! -d $@ ]; then \
		echo "Error: input directory containing reads was not found (DIR=$(DIR))" \
		exit -1; \
	fi

# Discover FASTQ files within given directory.

# Run fastqc.
${OUT}: ${DIR}
	# Create directory for reports.
	mkdir -p $@

	# Generate quality reports.
	${ENV_RUN} fastqc ${fastqc_opts} $(shell find ${DIR} -name "*.fastq" -o -name "*.fastq.gz") \
		|| echo "Error: no FASTQ files detected in ${DIR}"; exit -1

# Generate fastqc reports.
run: ${OUT}
	@ls -lh ${OUT}

# Remove generated reports.
run!:
	rm -rf ${OUT}

# Alternative rule for run!
clean: run!

# Run test suite.
test::
	make -f $(dir ${ROOT_PATH})fetch/sra.mk SRR=SRR1553425 N=1000 MODE=PE run
	make -f ${ROOT_PATH}/fastqc.mk DIR=${DIR} OUT=${OUT} run! run

DEPS := fastqc
# Command for installing dependencies.
install::
	micromamba create -n ${ENV}
	${ENV_RUN} micromamba install ${DEPS} --yes

# Rules without target files.
.PHONY: help run run! clean test install
