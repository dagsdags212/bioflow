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
ENV = bf-qc

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Directory containing reads to include in QC.
DIR ?= reads

# Discover FASTQ files within given directory.
READS = $(shell find $(DIR) -type f -name "*.fastq" -o -name "*.fastq.gz")

# Contaminant file.
CONTAMINANT_FILE ?=

# Output directory to store report.
OUT ?= fastqc

# Number of worker threads.
THREADS ?= 2

# Display usage.
help::
	@echo "#"
	@echo "# fastqc.mk: generate read quality reports with FASTQC"
	@echo "#"
	@echo "# READS=$(READS)"
	@echo "#"
	@echo "# OUT=$(OUT)"
	@echo "#"
	@echo "# make run|clean|install"

# Fastqc options.
fastqc_opts ?= $(if $(CONTAMINANT_FILE),-c $(CONTAMINANT_FILE)) \
							 -t $(THREADS) -o $(OUT)

# Compose fastqc command.
fastqc_cmd = fastqc $(fastqc_opts) $(READS)

# Directory storing reads must exist.
$(DIR):
	@echo "# Error: input directory containing reads was not found (DIR=$(DIR))"
	@exit -1

# Run fastqc.
$(OUT): $(READS) $(DIR)
	mkdir -p $(OUT)
	$(ENV_RUN) $(fastqc_cmd)

# Generate fastqc reports.
run: $(OUT)
	@ls -lh $(OUT)

# Remove generated reports.
run!:
	rm -rf $(OUT)

# Alternative rule for run!
clean: run!

# Run test suite.
test::
	make -f $(dir $(ROOT_PATH))fetch/sra.mk SRR=SRR1553425 N=1000
	make -f $(ROOT_PATH)/fastqc.mk run! run DIR=$(DIR) OUT=$(OUT)

# Command for installing dependencies.
install::
	@echo micromamba install fastqc

# Rules without target files.
.PHONY: help run run! clean install
