#
# Generates SNP calls with bcftools.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Absolute path of parent directory.
ROOT_PATH = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Micromamba environment.
ENV = bf-qc

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# A root to derive output default names from.
SRR = SRR1553425

# Number of worker threads.
THREADS ?= 2

# Genbank accession number.
ACC ?= AF086833

# The reference genome.
REF ?= refs/$(ACC).fa

# The alignment file.
BAM ?= bam/$(SRR).bam

# The variant file.
VCF ?= vcf/$(notdir $(basename $(BAM))).vcf.gz

# Additional bcf flags for pileup annotation.
PILE_FLAGS =  -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'

# Additional bcf flags for calling process.
CALL_FLAGS = --ploidy 2 --annotate 'FORMAT/GQ'


# The first target is always the help.
help::
	@echo "#"
	@echo "# bcftools.mk: calls variants using bcftools"
	@echo "#"
	@echo "# REF=$(REF)"
	@echo "# BAM=$(BAM)"
	@echo "# VCF=$(VCF)"
	@echo "#"
	@echo "# make run|install|clean"
	@echo "#"

$(VCF): $(BAM) $(REF)
	# Create directory for VCF output.
	mkdir -p $(dir $@)

	# Call SNPs, normalize, and sort.
	bcftools mpileup $(PILE_FLAGS) -O u -f $(REF) $(BAM) | \
		bcftools call $(CALL_FLAGS) -mv -O u | \
		bcftools norm -f $(REF) -d all -O u | \
		bcftools sort -O z > $@

# Generate index for VCF file.
$(VCF).tbi: $(VCF)
	bcftools index -t -f $<

# Generat summary statistics for VCF file.
$(VCF).stats: $(VCF)
	bcftools stats $< > $@

# Invoke VCF indexing.
run:: $(VCF).tbi $(VCF).stats
	@ls -lh $(VCF)

# Remove VCF file and its index.
run!::
	rm -rf $(VCF) $(VCF).tbi $(VCF).stats

# Alternative rule for run!
clean:: run!

# Print installation instructions.
install::
	@echo micromamba install bcftools

# TODO: Test the entire pipeline.
# test:

# Targets that are not files.
.PHONY: run run! install help test clean
