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
	@echo

# Display available parameters
params:
	@echo
	@echo "Assembly settings"
	@echo "  READ_DIR          path to directory containing reads"
	@echo "  REF               path to reference genomes"
	@echo "  PE                download reads in pair-end mode (default: true)"
	@echo "  ASSEMBLER         program of choice to perform assembly (default: none)"
	@echo "  MIN_CONTIG_LEN    minimum contig length to retain (default: 200)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-assembly)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Perform assembly with a specified tool.
assemble:
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

clean:
	rm -rf output/megahit output/spades output/minia
