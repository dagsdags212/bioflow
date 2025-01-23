#
# Generates SNP calls with freebayes.
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

# Number of worker threads.
THREADS ?= 2

# The alignment file.
BAM ?= bam/results.bam

# The reference genome.
REF ?= refs/genome.fa

# Name of the VCF with all variants.
VCF ?= vcf/$(notdir $(basename $(BAM))).freebayes.vcf.gz

# Additional flags passed to freebayes
freebayes_opts = --pvar 0.5


# The first target is always the help.
help::
	@echo "#"
	@echo "# freebayes.mk: calls variants using freebayes"
	@echo "#"
	@echo "# REF=$(REF)"
	@echo "# BAM=$(BAM)"
	@echo "# VCF=$(VCF)"
	@echo "#"
	@echo "# make run"
	@echo "#"

# Reference file must exist.
$(REF):
	@echo "# Error: Reference file not found (REF=$(REF))"
	@exit -1

# Read 1 must exist.
$(BAM):
	@echo "# Error: BAM alignment file not found (BAM=$(BAM))"
	@exit -1

# Call SNPs with freebayes.
$(VCF): $(BAM) $(REF)
	# Create output directory for VCF file.
	mkdir -p $(dir $@)

	# Call SNPs and normalize.
	freebayes $(freebayes_opts) -f $(REF) $(BAM) | bcftools norm -f $(REF) -d all -O z  > $(VCF)

# The VCF index file.
$(VCF).tbi: $(VCF)
	# Generate index for VCF file.
	bcftools index -t -f $<

# Invoke SNP calling.
run:: $(VCF).tbi
	@ls -lh $(VCF)

# Delete VCF file.
run!::
	rm -rf $(VCF)

# TODO: Test the entire pipeline.
# test:

# Print installation instructions.
install::
	@echo micromamba install bcftools freebayes

# Targets that are not files.
.PHONY: run run! install help test clean
