# 
# Perform reference-based mapping of sequencing reads.
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
.PHONY: help params init clean

# Formatting variables
dot := .
comma := ,
empty := 
space := $(empty) $(empty)

# Number of cores
THREADS ?= 8

# Environment manager
ENV_MANAGER := micromamba

# Conda environment
ENV := bwf-mapping

# Check if dependencies are installed
dependencies := bwa bowtie2 samtools qualimap

# Tool of choice for read mapping
MAPPER ?= bwa

# Directory containing reads
READ_DIR ?=

# Path of reference file
REF ?=

# Specify if reads are pair-end
PE ?= false

# Name and extension for output alignment file
OUTPUT ?=

# Alignment file format
OUTFMT ?= sam

# Name of alignment file
OUTNAME := output/bwa/$(if $(OUTPUT),$(OUTPUT)$(dot)sorted$(dot)$(OUTFMT),aln$(dot)sorted$(dot)$(OUTFMT))

# Separate reads if pair-end, otherwise store as single FASTQ file
ifeq ($(PE),true)
R1 = $(shell find $(READ_DIR) -type f -name "*_1.fastq*")
R2 = $(shell find $(READ_DIR) -type f -name "*_2.fastq*")
else
R = $(shell find $(READ_DIR) -type f -not -name "*_[12].fastq*")
endif

# MAPQ score threshold
MINQUAL ?= 20

# Display help message
help:
	@echo
	@echo "mapping.mk: align reads to a reference genome"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/mapping.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  map        - map reads to a reference using your mapper of choice"
	@echo "  stats      - generate mappings statistics for all alignment files"
	@echo "  evaluate   - compute MAPQ scores of mapped reads using qualimap"
	@echo "  visualize  - visualize SAM/BAM alignments in IGV"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Display available parameters
params:
	@echo
	@echo "Mapping settings"
	@echo "  READ_DIR          path of directory containing reads"
	@echo "  REF               path to reference file"
	@echo "  PE                if true, align pair-end reads (default: false)"
	@echo "  SORT         	 	 sort reads in alignment file"
	@echo "  MAPPER         	 program of choice to perform read mapping (default: bwa)"
	@echo "  MINQUAL         	 minimum MAPQ score for filtering (default: 20)"
	@echo "  OUTPUT         	 filename for output SAM file"
	@echo "  OUTFMT         	 file format for alignment file (default: sam)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-mapping)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

$(REF).bwt:
	# Index reference file
	bwa index $(REF)

map: $(REF).bwt
ifeq ($(MAPPER),bwa)
	@mkdir -p output/bwa
ifeq ($(PE),true)
	# Map, sort, and convert pair-end reads to reference
	bwa mem $(REF) $(R1) $(R2) | samtools sort | samtools view -h -O $(OUTFMT) > $(OUTNAME)
else
	# Map, sort, and convert single-end reads to reference
	bwa mem $(REF) $(R) | samtools sort | samtools view -h -O $(OUTFMT) > $(OUTNAME)
endif
endif

stats:
	# Print mapping statistics for all SAM/BAM files
	for aln in $(shell find output/ -maxdepth 2 -name "*.sam" -o -name "*.bam"); do \
		echo; \
		echo "=============== $$aln ==============="; \
		samtools flagstat $$aln; \
	done

evaluate:
	@mkdir -p output/qualimap
	for bam in $(shell find output/ -maxdepth 2 -name "*sorted.bam"); do \
		qualimap bamqc -outdir output/qualimap -outformat HTML -bam $$bam; \
	done

clean:
	rm -rf output/bwa output/qualimap
