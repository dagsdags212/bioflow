#
# Map reads against a reference using bwa.
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
ENV = bf-mapping

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 4

# Bowtie2 options.
FLAGS := -t $(THREADS)

# First in read pair.
R1 ?=

# Second in read pair, if any.
R2 ?=

# The reference genome.
REF ?=

# The alignment file.
BAM ?=

# A directory to stores file for the bowtie index.
IDX_DIR ?= $(dir ${REF})idx

# The bowtie2 index prefix.
IDX ?= $(IDX_DIR)/$(notdir $(REF))

# File in the index directory.
IDX_FILE = $(IDX).acc

# Set attributes for the read group.
ID ?= run1
SM ?= sample1
LB ?= library1
PL ?= ILLUMINA
RG ?= "@RG\tID:$(ID)\tSM:$(SM)\tLB:$(LB)\tPL:$(PL)"

# The name of the stats file.
STATS = $(basename $(BAM)).stats

# Print the help message.
help::
	@echo "#"
	@echo "# bwa.mk: align reads using BWA"
	@echo "#"
	@echo "# R1=${R1}"
	@echo "# R2=${R2}"
	@echo "# REF=${REF}"
	@echo "# IDX=${IDX} (optional)"
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

# Generate bwa index for reference.
$(IDX_FILE): $(REF)
	# Create output directory for bwa index.
	mkdir -p $(dir $@)

	# Generate bwa index for the reference.
	bwa index -p $(IDX) $<

# Invoke reference indexing.
index: $(REF) $(IDX_FILE)
	@echo "# bwa index: $(IDX)"

# Remove the index.
index!:
	rm -f $(IDX_FILE)

# Compose bowtie2 command.
BWA_CMD = bwa mem $(FLAGS) -R $(RG) $(R1) $(R2)

# Generate a sorted alignment file.
$(BAM): $(R1) $(R2)
	# Output directory for BAM files.
	mkdir -p $(dir $@)

	# Perform read mapping.
	$(ENV_RUN) $(BWA_CMD) | samtools sort -@ $(THREADS) > $@

# Create the BAM index.
$(BAM).bai: $(BAM)
	samtools index $(BAM)

# Invoke the read mapping rule.
align: $(BAM).bai
	@ls -lh $(BAM)

# Alternative rule for align.
run: align

# Generate alignment statistics.
$(STATS): $(BAM).bai
	samtools flagstat $(BAM) > $(STATS)

# Trigger stats generation.
stats: $(STATS)
	@echo "# $(STATS)"
	@cat $(STATS)

# Remove output BAM files.
run!:
	rm -rf $(BAM) $(BAM).bai

# Alternative rule for run!
clean: run!

# Show installation command.
install::
	@echo "micromamba install bwa samtools"

# Targets that are not files.
.PHONY: run run! install help index test
