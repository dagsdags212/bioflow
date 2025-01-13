# 
# Detect variation from reads based on a reference genome
#

# import config variables
include src/_config.mk

# import global variables
include src/_globals.mk

.PHONY: help params init clean split

# Project root
ROOT_DIR = $(shell dirname $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST)))))

# Conda environment
ENV := bf-vc

# Path to conda environment
ENV_DIR = $(shell $(ENV_MANAGER) info | grep "envs directories" | cut -d ":" -f 2 | xargs)/$(ENV)

# Check if dependencies are installed
dependencies := freebayes samtools bamtools bedtools bcftools vcflib rtg-tools 

# Path to reference genome
REF ?=

# Path to binary alignment file 
BAM ?=

# First word in BAM filename
BAM_ROOT = $(notdir $(word 1,$(subst $(dot),$(space),$(BAM))))

# Minimum quality for filtering variants
MINQUAL ?= 20

# Flags for rtg vcffilter subcommand
vcffilter_opts := --no-gzip -q $(MINQUAL)

# Freebayes flags
PLOIDY ?= 2

freebayes_opts := -p $(PLOIDY)

# Display help message
help:
	@echo
	@echo "variant_calling.mk: call variants based on a reference"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/variant_calling.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  call       - perform variant calling of mapped reads based on a reference genome"
	@echo "  filter     - only retain variants with scores above threshold"
	@echo "  stats      - generate variant calling metrics"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo
	@echo "Variant calling settings"
	@echo "  REF               path to reference genome"
	@echo "  BAM               path to binary alignment file"
	@echo "  MINQUAL           set variant score threshold for filtering (default: 10)"
	@echo "  PLOIDY            ploidy number of reference (default: 1)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-mapping)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Create index for reference genome
$(REF).fai: $(REF)
	@echo "Indexing reference genome: $<"
	samtools faidx $(REF)

# Index BAM file
$(BAM).bai: $(BAM)
	@echo "Indexing BAM file: $<"
	bamtools index -in $(BAM)

# Call variants
call: $(REF).fai $(BAM).bai
	@echo "Calling variants using freebayes"
	@mkdir -p output/freebayes
	freebayes -f $(REF) $(freebayes_opts) $(BAM) > output/freebayes/$(BAM_ROOT).vcf
	@echo "DONE"
	@echo "Output: output/freebayes/$(BAM_ROOT).vcf"

# Produce summary metrics for VCF file and graph
stats: output/freebayes/$(BAM_ROOT).vcf
	@echo "Generating summary statistics for $<"
	bcftools stats -F $(REF) -s - $< > $<.stats
	@echo "DONE"
	@echo "Output: $<.stats"

filter: output/freebayes/$(BAM_ROOT).vcf
	@mkdir -p output/rtg/
	@echo "Filtering $<: only retaining variants with scores of at least $(MINQUAL)"
	rtg vcffilter -i $< -o output/rtg/$(BAM_ROOT).q$(MINQUAL).vcf $(vcffilter_opts)
	@echo "DONE"
	@echo "Output: output/rtg/$(BAM_ROOT).q$(MINQUAL).vcf"

clean:
	rm -rf output/freebayes/ output/rtg/
	rm -rf $(REF).fai
