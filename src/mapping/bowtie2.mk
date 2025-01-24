#
# Map reads against a reference using bowtie2.
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

# Number of worker threads.
THREADS ?= 4

# Bowtie2 options.
FLAGS := --sensitive-local

# First in read pair.
R1 ?=

# Second in read pair, if any.
R2 ?=

# Set MODE to pair-end (PE) if R2 is provided.
ifeq ($(R2),)
	MODE ?= SE
else
	MODE ?= PE
endif

# The reference genome.
REF ?= refs/genome.fa

# The alignment file.
BAM ?= bam/aln.bam

# A directory to stores file for the bowtie index.
IDX_DIR ?= $(dir ${REF})idx

# The bowtie2 index prefix.
IDX ?= $(IDX_DIR)/$(notdir $(REF))

# File in the index directory.
IDX_FILE = $(IDX).1.bt2

# Set attributes for the read group.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA
RG ?= "@RG\tID:$(ID)\tSM:$(SM)\tLB:$(LB)\tPL:$(PL)"

# Print the help message.
help::
	@echo "#"
	@echo "# bowtie2.mk: aligns read using bowtie2 "
	@echo "#"
	@echo "# REF=${REF}"
	@echo "# IDX=${IDX}"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "#"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# make index|run|test|clean"
	@echo "#"

# Read 1 must exist.
$(R1):
	@echo "# Error: Read 1 not found (R1=$(R1))"
	@exit -1

# If set, read 2 must also exist.
ifneq ($(R2),)
$(R2):
	@echo "# Error: Read 2 not found (R2=$(R2))"
endif

# Reference file must exist.
$(REF):
	@if [ ! -f $(REF) ]; then
		echo "# Error: Reference file not found (REF=$(REF))"
		exit -1
	fi

# Generate bowtie index for reference.
$(IDX_FILE):
	# Create output directory for bowtie index.
	mkdir -p $(dir $(IDX_FILE))

	# Generate bowtie2 index for the reference.
	bowtie2-build $(REF) $(IDX)

# Invoke reference indexing.
index: $(REF) $(IDX_FILE)
	@echo "# bowtie2 index: $(IDX)"

# Remove the index.
index!:
	rm -f $(IDX_FILE)

# Compose bowtie2 command.
ifeq ($(R2),)
	BOWTIE2_CMD ?= bowtie2 $(FLAGS) -p $(THREADS) -x $(IDX) -U $(R1)
else
	BOWTIE2_CMD ?= bowtie2 $(FLAGS) -p $(THREADS) -x $(IDX) -1 $(R1) -2 $(R2)
endif

# Generate a sorted alignment file.
$(BAM): $(R1) $(R2)
	# Output directory for BAM files.
	mkdir -p $(dir $@)

	# Perform read mapping.
	$(ENV_RUN) $(BOWTIE2_CMD) | samtools sort -@ $(THREADS) > $@

# Create the BAM index.
$(BAM).bai: $(BAM)
	samtools index $(BAM)

# Invoke the read mapping rule.
align: $(BAM).bai
	@ls -lh $(BAM)

# Alternative rule for align.
run: align

# Remove output BAM files.
run!:
	rm -rf $(BAM) $(BAM).bai

# Alternative rule for run!
clean: run!

# Show installation command.
install::
	@echo "micromamba install bowtie2 samtools"

# Targets that are not files.
.PHONY: run run! install help index test
