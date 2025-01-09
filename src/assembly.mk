# 
# Perform de novo and reference-based sequence assembly and annotation.
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
.PHONY: help init assemble evaluate visualize

# Directory containing FASTQ reads
READ_DIR ?=

# Path to reference genomes
REF ?=

# Specify if reads are pair-end
PE ?= true

# Separate reads if pair-end, otherwise store as single FASTQ file
ifeq ($(PE),true)
R1 = $(shell find $(READ_DIR) -type f -name "*_1*fastq*")
R2 = $(shell find $(READ_DIR) -type f -name "*_2*fastq*")
else
R = $(shell find $(READ_DIR) -type f -name "*fastq*")
endif

# Number of cores
THREADS ?= 8

# Environment manager
ENV_MANAGER ?= micromamba

# Conda environment
ENV ?= bwf-assembly

# Check if dependencies are installed
dependencies := megahit spades quast

# Specify assembly program
ASSEMBLER ?=
ASSEMBLERS := megahit spades

# Minimum contig length to keep
MIN_CONTIG_LEN ?= 200

# Spades flags
spades_opts := --careful -t $(THREADS)

# Megahit flags
megahit_opts := -t $(THREADS) --min-contig-len $(MIN_CONTIG_LEN) --cleaning-rounds 5

# Quast flags
quast_opts := -t $(THREADS) --circos

# Display help message
help:
	@echo
	@echo "assembly.mk: assemble contigs from sequencing reads"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/qc.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  assemble   - generate sequence assembly using program of choice"
	@echo "  evaluate   - compute contig statistics and compare assemblies"
	@echo "  visualize  - visualize assembly with bandage"

# Display available parameters
params:
	@echo "Global settings"
	@echo "  READ_DIR          path to directory containing reads"
	@echo "  REF               path to reference genomes"
	@echo "  THREADS           number of cores (default: 4)"
	@echo "Assembly settings"
	@echo "  ASSEMBLER         program of choice to perform assembly (default: none)"
	@echo "  MIN_CONTIG_LEN    minimum contig length to retain (default: 200)"
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-qc)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Generate contigs from a set of reads
assemble:
# Assemble reads using megahit
ifeq ($(ASSEMBLER),megahit)
	mkdir -p assembly/megahit
	megahit -o assembly/megahit $(megahit_opts) -1 $(R1) -2 $(R2)
endif

# Assemble reads using spades
ifeq ($(ASSEMBLER),spades)
	mkdir -p assembly/spades
	spades.py -o assembly/spades $(spades_opts) -1 $(R1) -2 $(R2)
endif

# Assemble reads using multiple programs
ifndef ASSEMBLER
	# Run megahit
	mkdir -p assembly/megahit
	megahit -o assembly/megahit $(megahit_opts) -1 $(R1) -2 $(R2)
	# Run spades
	mkdir -p assembly/spades
	spades.py -o assembly/spades $(spades_opts) -1 $(R1) -2 $(R2)
endif

# Compute assembly statistics for evaluation
evaluate: $(foreach asm,$(ASSEMBLERS),assembly/$(asm)/*contigs*)
	# Discover contig files in assembly directory and run quast
ifeq ($(ASSEMBLER),megahit)
	quast $(quast_opts) $(shell find assembly/megahit/ -maxdepth 1 -type f -regex '.*contigs.fa')
endif

ifeq ($(ASSEMBLER),spades)
	quast $(quast_opts) $(shell find assembly/spades/ -maxdepth 1 -type f -regex '.*contigs.fasta')
endif

ifndef ASSEMBLER
	quast $(quast_opts) --labels "megahit, spades" \
		$(shell find assembly -maxdepth 2 -type f -regex '.*contigs.fa' -o -regex '.*contigs.fasta')
endif

# Generate assembly graphs using bandage
visualize:
	mkdir -p assembly/bandage
ifeq ($(ASSEMBLER),megahit)
	# Convert contig FASTA to fastg file
	megahit_toolkit contig2fastg 99 assembly/megahit/final.contigs.fa > assembly/megahit/assembly_graph.fastg
	Bandage image assembly/megahit/assembly_graph.fastg assembly/bandage/megahit_assembly.png
endif

ifeq ($(ASSEMBLER),spades)
	Bandage image assembly/spades/assembly_graph.fastg assembly/bandage/spades_assembly.png
endif

ifndef ASSEMBLER
	megahit_toolkit contig2fastg 99 assembly/megahit/final.contigs.fa > assembly/megahit/assembly_graph.fastg
	Bandage image assembly/spades/assembly_graph.fastg assembly/bandage/megahit_assembly.png
	Bandage image assembly/spades/assembly_graph.fastg assembly/bandage/spades_assembly.png
endif

