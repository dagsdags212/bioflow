# 
# Perform de novo and reference-based sequence assembly and annotation.
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# Absoulate path to conda environment files
ENV_ROOT := ${BIOFLOW_PREFIX}/envs

# Conda environment for module
ENV := bf-assembly

# Execute within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

.PHONY: help params init assemble visualize

# List of support assemblers (maintain lexicographic order)
ASSEMBLERS := minia megahit spades

# Display help message
help:
	@echo
	@echo "assembly.mk: build contigs from sequencing reads"
	@echo ""
	@echo "Usage:"
	@echo "  bf-assemble <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  assemble   - generate sequence assembly using program of choice"
	@echo "  params     - list all available options"
	@echo "  init       - download conda environment"
	@echo "  clean      - delete all output produced by this module"
	@echo

# Display available parameters
params:
	@echo
	@echo "Assembly settings"
	@echo "  ASSEMBLER         program of choice to perform assembly (default: megahit)"
	@echo "  PE1               path to FASTQ file for pair-end #1"
	@echo "  PE2               path to FASTQ file for pair-end #2"
	@echo "  SE                path to single-end FASTQ file"
	@echo "  K                 k-mer length to use"
	@echo "  MIN_CONTIG_LEN    minimum contig length to retain (default: 200)"
	@echo "  PLATFORM          specify sequencing platform for spades"
	@echo "  PREFIX            prepend prefix to output files of minia (default: final)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo "  OUTDIR            path to directory for storing output files"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-assembly)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -f $(ENV_ROOT)/bf-assembly.yml

# Perform assembly with a specified tool.
assemble:
# Assemble reads using megahit
ifeq ($(ASSEMBLER),megahit)
ifdef SE
	make -f $(MK_ROOT)/run_megahit.mk SE=$(SE) \
		THREADS=$(THREADS) MIN_CONTIG_LEN=$(MIN_CONTIG_LEN) OUTDIR=$(OUTDIR)
else ifdef PE1
	make -f $(MK_ROOT)/run_megahit.mk PE1=$(PE1) PE2=$(PE2) \
		THREADS=$(THREADS) MIN_CONTIG_LEN=$(MIN_CONTIG_LEN) OUTDIR=$(OUTDIR)
endif
endif

# Assemble reads using spades
ifeq ($(ASSEMBLER),spades)
ifdef SE
	make -f $(MK_ROOT)/run_spades.mk SE=$(SE) \
		THREADS=$(THREADS) K=$(K) PLATFORM=$(PLATFORM) OUTDIR=$(OUTDIR)
else ifdef PE1
	make -f $(MK_ROOT)/run_spades.mk PE1=$(PE1) PE2=$(PE2) \
		THREADS=$(THREADS) K=$(K) PLATFORM=$(PLATFORM) OUTDIR=$(OUTDIR)
endif
endif

# Assemble reads using minia
ifeq ($(ASSEMBLER),minia)
ifdef SE
	make -f $(MK_ROOT)/run_minia.mk SE=$(SE) \
		THREADS=$(THREADS) K=$(K) PREFIX=$(PREFIX) OUTDIR=$(OUTDIR)
else ifdef PE1
	make -f $(MK_ROOT)/run_minia.mk PE1=$(PE1) PE2=$(PE2) \
		THREADS=$(THREADS) K=$(K) PREFIX=$(PREFIX) OUTDIR=$(OUTDIR)
endif
endif

clean:
	rm -rf output/megahit output/spades output/minia
