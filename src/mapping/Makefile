# 
# Perform reference-based mapping of sequencing reads.
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# Absoulate path to conda environment files
ENV_ROOT := ${BIOFLOW_PREFIX}/envs

# Conda environment for module
ENV := bf-mapping

# Execute within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/mapping.mk

# Display help message
help:
	@echo
	@echo "mapping.mk: align reads to a reference genome"
	@echo ""
	@echo "Usage:"
	@echo "  bf-map <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  map        - map reads to a reference using your mapper of choice"
	@echo "  stats      - generate mappings statistics for all alignment files"
	@echo "  evaluate   - compute MAPQ scores of mapped reads using qualimap"
	@echo "  visualize  - visualize SAM/BAM alignments in IGV"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -f $(ENV_ROOT)/bf-mapping
	@# Extract bbtool scripts and add to env path
	bbmap_tar=$(ROOT_DIR)/tools/tar/BBMap_39.14.tar.gz
	tar -xzf $$bbmap_tar -C $(ROOT_DIR)/tools
	mv $(ROOT_DIR)/tools/bbmap/* $(ENV_DIR)/bin/
	rm -rf $(ROOT_DIR)/tools/bbmap/

# Display available parameters
params:
	@echo
	@echo "Mapping settings"
	@echo "  MAPPER         	 choose a mapping tool (default: bwa)"
	@echo "  REF               path to reference file"
	@echo "  PE1               path to FASTQ file for read-pair #1"
	@echo "  PE2               path to FASTQ file for read-pair #2"
	@echo "  SE                path to single-end FASTQ file"
	@echo "  MINQUAL         	 minimum MAPQ score for filtering (default: 20)"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo "  PREFIX            prepend to filename of resulting alignment file"
	@echo "  OUTDIR            path to directory for storing output"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-mapping)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Additional mapper options.
opts := $(if $(OUTDIR),OUTDIR=$(OUTDIR))
opts += $(if $(PREFIX),PREFIX=$(PREFIX))
opts += $(if $(THREADS),THREADS=$(THREADS))


# Perform read alignment with mapper of choice.
map:
# Map reads using bwa.
ifeq ($(MAPPER),bwa)
ifdef SE
	make -f $(MK_ROOT)/run_bwa.mk REF=$(REF) SE=$(SE) $(opts)
else ifdef PE1
	make -f $(MK_ROOT)/run_bwa.mk REF=$(REF) PE1=$(PE1) PE2=$(PE2) $(opts)
else
	@echo "Error: path to FASTQ reads not specified"
endif

# Map reads using bowtie.
else ifeq ($(MAPPER),bowtie)
ifdef SE
	make -f $(MK_ROOT)/run_bowtie.mk REF=$(REF) SE=$(SE) $(opts)
else ifdef PE1
	make -f $(MK_ROOT)/run_bowtie.mk REF=$(REF) PE1=$(PE1) PE2=$(PE2) $(opts)
else
	@echo "Error: path to FASTQ reads not specified"
endif
else
	@echo "Error: specified MAPPER is currently not supported"
endif

# Print mapping statistics for all SAM/BAM files
stats:
	@for aln in $(shell find $(OUTDIR) -type f -name "*.sam" -o -name "*.bam"); do \
		printf "\n=============== %s ===============\n\n" $$aln; \
		samtools flagstat $$aln; \
	done

evaluate:
	@mkdir -p output/qualimap
	for bam in $(shell find . -type f -name "*sorted.bam"); do \
		outdir=output/qualimap/$$(basename $${bam%%.sorted.bam}); \
		mkdir -p $${outdir}; \
		qualimap bamqc -outdir $${outdir} -outformat HTML -bam $$bam; \
	done

clean:
	rm -rf $(OUTDIR)

.PHONY: help params init clean
