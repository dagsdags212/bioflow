#
# Aggregate FASTQC reports into a single interface.
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

# Directory containing FASTQC reports.
DIR ?= fastqc

# Output directory to store report.
OUT ?= multiqc

# Display usage.
help::
	@echo "#"
	@echo "# multiqc.mk: consolidate multiple FASTQC reports"
	@echo "#"
	@echo "# DIR=$(DIR)"
	@echo "#"
	@echo "# OUT=$(OUT)"
	@echo "#"
	@echo "# make run|clean|install"

# Fastqc options.
multiqc_opts ?= --clean-up --verbose

# Compose fastqc command.
multiqc_cmd = multiqc $(multiqc_opts) -o $(OUT) $(DIR)

# Directory storing reads must exist.
$(DIR):
	@echo "# Error: input directory containing FASTQC reports not found (DIR=$(DIR))"
	@exit -1

# Run fastqc.
$(OUT): $(DIR)
	mkdir -p $(OUT)
	$(ENV_RUN) $(multiqc_cmd)

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
	make -f $(dir $(ROOT_PATH))fetch/sra.mk SRR=SRR1553425 N=1000 DIR=reads
	make -f $(ROOT_PATH)/fastqc.mk run! run DIR=reads OUT=$(DIR)
	make -f $(ROOT_PATH)/multiqc.mk run! run DIR=$(DIR) OUT=$(OUT)

# Command for installing dependencies.
install::
	@echo micromamba install multiqc

# Rules without target files.
.PHONY: help run run! clean install
