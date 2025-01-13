# 
# Perform reference-based mapping of sequencing reads.
#

# import config variables
include src/_config.mk

# import global variables
include src/_globals.mk

.PHONY: help params init clean

# Project root
ROOT_DIR = $(shell dirname $(shell dirname $(realpath $(MAKEFILE_LIST))))

# Conda environment
ENV := bf-mapping

# Path to conda environment
ENV_DIR = $(shell $(ENV_MANAGER) info | grep "envs directories" | cut -d ":" -f 2 | xargs)/$(ENV)

# Check if dependencies are installed
dependencies := bwa bowtie2 samtools qualimapMSG

# Tool of choice for read mapping
MAPPER ?= bwa

# Directory containing reads
READ_DIR ?=

# Path of reference file
REF ?=

# Specify if reads are pair-end
PE ?= false

# Name and extension for output alignment file
OUTPUT ?= aln

# Alignment file format
OUTFMT ?= sam

# Name of alignment file
OUTNAME := $(if $(OUTPUT),$(OUTPUT).sorted.$(OUTFMT),aln.sorted.$(OUTFMT))

# Set required index files based on specified mapping tool
BWA_SFX = amb ann bwt pac sa
BOWTIE_SFX = 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2

ifeq ($(MAPPER),bwa)
IDX = $(addprefix $(REF).,$(BWA_SFX))
else ifeq ($(MAPPER),bowtie)
IDX = $(addprefix $(REF).,$(BOWTIE_SFX))
endif

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
	$(ENV_MANAGER) create -n $(ENV) $(dependencies) --yes
	@# Extract bbtool scripts and add to env path
	bbmap_tar=$(ROOT_DIR)/tools/tar/BBMap_39.14.tar.gz
	tar -xzf $$bbmap_tar -C $(ROOT_DIR)/tools
	mv $(ROOT_DIR)/tools/bbmap/* $(ENV_DIR)/bin/
	rm -rf $(ROOT_DIR)/tools/bbmap/

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

$(IDX): $(REF)
	# Index reference file
ifeq ($(MAPPER),bwa)
	bwa index $(REF)
else ifeq ($(MAPPER),bowtie)
	bowtie2-build $(REF) $(basename $(REF))
endif

map: $(IDX)
ifeq ($(MAPPER),bwa)
	@mkdir -p output/bwa
ifeq ($(PE),true)
	# Map, sort, and convert pair-end reads to reference
	bwa mem $(REF) $(R1) $(R2) | samtools sort - | samtools view -h -O $(OUTFMT) -o output/bwa/$(OUTNAME)
else
	# Map, sort, and convert single-end reads to reference
	bwa mem $(REF) $(R) | samtools sort - | samtools view -h -O $(OUTFMT) -o output/bwa/$(OUTNAME)

endif

else ifeq ($(MAPPER),bowtie)
	@mkdir -p output/bowtie
ifeq ($(PE),true)
	bowtie2 -x $(basename $(REF))	-1 $(R1) -2 $(R2) | samtools sort - | samtools view -h -O $(OUTFMT) -o output/bowtie/$(OUTNAME)
else
	bowtie2 -x $(basename $(REF)) -U $(R) | samtools sort - | samtools view -h -O $(OUTFMT) -o output/bowtie/$(OUTNAME)
endif

endif

# Print mapping statistics for all SAM/BAM files
stats:
	@for aln in $(shell find output/ -maxdepth 2 -name "*.sam" -o -name "*.bam"); do \
		echo; \
		echo "=============== $$aln ==============="; \
		samtools flagstat $$aln; \
	done

evaluate:
	@mkdir -p output/qualimap
	for bam in $(shell find output/ -maxdepth 2 -name "*sorted.bam"); do \
		outdir=output/qualimap/$$(basename $${bam%%.*}); \
		mkdir -p $${outdir}; \
		qualimap bamqc -outdir $${outdir} -outformat HTML -bam $$bam; \
	done

clean:
	rm -rf output/bwa/ output/bowtie/ output/qualimap/
