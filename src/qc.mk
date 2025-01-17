# 
# Perform quality checks and filtering on sequencing data.
#

# Absoluate path for the `bioflow/src` directory
SRC_ROOT := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

# Absolute path for `bioflow` directory
PROJECT_ROOT := $(shell dirname $(SRC_ROOT))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := $(PROJECT_ROOT)/config

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/_tools.mk

.PHONY: help params init clean

# Project root
ROOT_DIR = $(shell dirname $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST)))))

# Conda environment
ENV = bf-qc

# Path to conda environment
ENV_DIR = $(shell $(ENV_MANAGER) info | grep "envs directories" | cut -d ":" -f 2 | xargs)/$(ENV)

# Check if dependencies are installed
dependencies = fastqc multiqc trimmomatic fastp

# Run a command within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

# Target directory to run fastqc (default: READ_DIR)
TARGET ?= 

# fastp parameters
MINLEN ?= 30
MAXLEN ?= 150
MINQUAL ?= 20
PE ?= true

ifeq ($(PE),true)
R1 = $(shell find $(READ_DIR) -type f -name "*_1.fastq*")
R2 = $(shell find $(READ_DIR) -type f -name "*_2.fastq*")
else
R = $(shell find $(READ_DIR) -type f -not -name "*_[12].fastq*")
endif

# fastqc options
fastqc_opts := --svg --threads $(THREADS)

# fastp options
fastp_opts := --overrepresentation_analysis --correction --cut_right
fastp_opts += --length_required $(MINLEN) --length_limit $(MAXLEN) -q $(MINQUAL)
fastp_opts += --html fastp/reports/fastp_report_$(shell date +%y%m%d%H%M%S).html
fastp_opts += --json fastp/reports/fastp_report_$(shell date +%y%m%d%H%M%S).json

# Display help message
help:
	@echo
	@echo "qc.mk: perform quality control on sequencing data"
	@echo
	@echo "Usage:"
	@echo "  make -f src/qc.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  fastqc  - generate FASTQC report for a set of reads"
	@echo "  multiqc - consolidate FASTQC output into a single report"
	@echo "  fastp   - perform adapter trimming and quality filtering"
	@echo

# Display available parameters
params:
	@echo
	@echo "Trimming and filtering"
	@echo "  TARGET         directory path containing read data"
	@echo "  MINLEN         minimum read length (default: 30)"
	@echo "  MAXLEN         maximum read length (default: 150)"
	@echo "  MINQUAL        minimum acceptable quality score (default: 20)"
	@echo "  PE             treat data as pair-end reads (default: true)"
	@echo
	@echo "Global settings"
	@echo "  THREADS        number of cores (default: 4)"
	@echo
	@echo "Environment settings"
	@echo "  ENV            environment name (default: bwf-qc)"
	@echo "  ENV_MANAGER    environment manager (default: micromamba)"
	@echo

# Create new self-contained environment
init:
	@$(ENV_MANAGER) create -n $(ENV) $(dependencies) --yes
	# Extract bbtool scripts and add to env path
	bbmap_tar=$(ROOT_DIR)/tools/tar/BBMap_39.14.tar.gz
	tar -xzf $$bbmap_tar -C $(ROOT_DIR)/tools
	mv $(ROOT_DIR)/tools/bbmap/* $(ENV_DIR)/bin/
	rm -rf $(ROOT_DIR)/tools/bbmap/

output/fastqc output/multiqc output/fastp:
	@mkdir -p $@

# Discover FASTQ files in specified directory to be used by FASTQC
fastqc: READS = $(shell find $(if $(TARGET),$(TARGET),$(READ_DIR)) -type f -name "*.fastq" -o -name "*.fastq.gz")

# Generate fastqc report for a set of reads
fastqc: output/fastqc
	@outdir=output/fastqc/report$(shell ls output/fastqc | wc -l)
	@echo $${outdir} > /tmp/most_recent_fastqc_report.txt
	@mkdir -p $${outdir}
	@$(RUN) fastqc -o $${outdir} --threads $(THREADS) $(READS)

# Consolidate fastqc files into a single report
multiqc: /tmp/most_recent_fastqc_report.txt output/multiqc
ifeq ($(TARGET),all)
	@$(RUN) multiqc -o multiqc/all output/fastqc
else
	@$(RUN) multiqc -o $(subst fastqc,multiqc,$(shell cat $<)) $(shell cat $<)
endif

fastp: output/fastp $(wildcard $(READ_DIR)/*.fastq*)
	@mkdir -p output/fastp/reads
	@mkdir -p output/fastp/reports
ifeq ($(PE),true)
	$(RUN) fastp -i $(R1) -I $(R2)	$(fastp_opts) -o output/fastp/reads/trimmed_$(shell basename $(R1)) -O output/fastp/reads/trimmed_$(shell basename $(R2))
else
	$(RUN) fastp -i $(R)	$(fastp_opts) -o output/fastp/reads/trimmed_$(shell basename $(R))
endif

clean:
	rm -rf output/fastqc/ output/multiqc/ output/fastp/
