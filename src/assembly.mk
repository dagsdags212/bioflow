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
R1 = $(shell find $(READ_DIR) -type f -name "*_1.fastq*")
R2 = $(shell find $(READ_DIR) -type f -name "*_2.fastq*")
else
R = $(shell find $(READ_DIR) -type f -not -name "*_[12].fastq*")
endif

# Number of cores
THREADS ?= 8

# Environment manager
ENV_MANAGER ?= micromamba

# Conda environment
ENV ?= bwf-assembly

# Check if dependencies are installed
dependencies := megahit spades quast minia

# List of support assemblers (maintain lexicographic order)
ASSEMBLERS := minia megahit spades

# Specify assembly program
ASSEMBLER ?=

# Minimum contig length to keep
MIN_CONTIG_LEN ?= 200

# Spades flags
spades_opts := --careful -t $(THREADS)

# Megahit flags
megahit_opts := -t $(THREADS) --min-contig-len $(MIN_CONTIG_LEN) --cleaning-rounds 5

# Minia flags
minia_opts := -nb-cores $(THREADS) -traversal contig -kmer-size 31

# Quast flags
comma := ,
null :=
space := $(null) $(null)
quast_opts := -t $(THREADS) --circos --labels "$(subst $(space),$(comma),$(ASSEMBLERS))"

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
	@echo "  PE                download reads in pair-end mode (default: true)"
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
	if [ -d "output/megahit"]; then rm -rf output/megahit; fi
	megahit -o output/megahit $(megahit_opts) -1 $(R1) -2 $(R2)
endif

# Assemble reads using spades
ifeq ($(ASSEMBLER),spades)
	mkdir -p output/spades
	spades.py -o output/spades $(spades_opts) -1 $(R1) -2 $(R2)
endif

# Assemble reads using minia
ifeq ($(ASSEMBLER),minia)
ifeq ($(PE),false)
	mkdir -p output/minia
	minia -in $(R) -out output/minia/final $(minia_opts)
else
	@echo "minia only supports single-ended reads as input"
	@echo "download SE reads using the fetch.mk workflow"
endif
endif

# Assemble reads using multiple programs
ifndef ASSEMBLER
ifeq ($(PE),true)
	# Run megahit
	if [ -d "output/megahit"]; then rm -rf output/megahit; fi
	megahit -o output/megahit $(megahit_opts) -1 $(R1) -2 $(R2)
	# Run spades
	mkdir -p output/spades
	spades.py -o output/spades $(spades_opts) -1 $(R1) -2 $(R2)
endif
	# Run minia
ifeq ($(PE),false)
	mkdir -p output/minia
	minia -in $(R) -out output/minia/final $(minia_opts)
else
	@echo "minia only supports single-ended reads as input"
	@echo "download SE reads using the fetch.mk workflow"
endif
endif

# Compute assembly statistics for evaluation
evaluate: $(foreach asm,$(ASSEMBLERS),output/$(asm)/*contigs*)
	# Discover contig files in assembly directory and run quast
ifeq ($(ASSEMBLER),megahit)
	quast $(quast_opts) $(shell find output/megahit/ -maxdepth 1 -type f -regex '.*contigs.fa')
endif

ifeq ($(ASSEMBLER),spades)
	quast $(quast_opts) $(shell find output/spades/ -maxdepth 1 -type f -regex '.*contigs.fasta')
endif

ifndef ASSEMBLER
	quast -o output/quast $(quast_opts) $(shell find output -maxdepth 2 -type f -name '*contigs.fasta' -o -name '*contigs.fa')
endif

# Generate assembly graphs using bandage
visualize:
	mkdir -p output/bandage
ifeq ($(ASSEMBLER),megahit)
	# Convert contig FASTA to fastg file
	megahit_toolkit contig2fastg 99 output/megahit/final.contigs.fa > output/megahit/assembly_graph.fastg
	Bandage image output/megahit/assembly_graph.fastg output/bandage/megahit_assembly.png
endif

ifeq ($(ASSEMBLER),spades)
	Bandage image output/spades/assembly_graph.fastg output/bandage/spades_assembly.png
endif

ifndef ASSEMBLER
	megahit_toolkit contig2fastg 99 output/megahit/final.contigs.fa > output/megahit/assembly_graph.fastg
	Bandage image output/spades/assembly_graph.fastg output/bandage/megahit_assembly.png
	Bandage image output/spades/assembly_graph.fastg output/bandage/spades_assembly.png
endif

