# 
# Perform quality checks and filtering on sequencing data.
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
.PHONY: help init fastqc fastp

READ_DIR ?=
THREADS ?= 4

# Discover FASTQ files in specified directory to be used by FASTQC
fastqc: READS = $(shell find $(READ_DIR) -type f -name "*fastq*")

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

fastp_opts = --overrepresentation_analysis --correction --cut_right
fastp_opts += --length_required $(MINLEN) --length_limit $(MAXLEN) -q $(MINQUAL)
fastp_opts += --html fastp/reports/fastp_report_$(shell date +%y%m%d%H%M%S).html
fastp_opts += --json fastp/reports/fastp_report_$(shell date +%y%m%d%H%M%S).json

# Environment manager
ENV_MANAGER ?= micromamba

# Conda environment
ENV ?= bwf-qc

# Check if dependencies are installed
dependencies := fastqc multiqc trimmomatic fastp

# Display help message
help:
	@echo ""
	@echo "qc.mk: perform quality control on sequencing data"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/qc.mk <command> [options]"
	@echo ""
	@echo "COMMANDS:"
	@echo "  fastqc  - generate FASTQC report for a set of reads"
	@echo "  multiqc - consolidate FASTQC output into a single report"
	@echo "  fastp   - perform adapter trimming and quality filtering"

# Display available parameters
params:
	@echo "Global settings"
	@echo "  READ_DIR       directory path containing read data"
	@echo "  PE       			specify reads are pair-end (default: true)"
	@echo "  THREADS        number of cores (default: 4)"
	@echo "Environment settings"
	@echo "  ENV            environment name (default: bwf-qc)"
	@echo "  ENV_MANAGER    environment manager (default: micromamba)"
	@echo "Trimming and filtering"
	@echo "  MINLEN         minimum read length (default: 30)"
	@echo "  MAXLEN         maximum read length (default: 150)"
	@echo "  MINQUAL        minimum acceptable quality score (default: 20)"
	@echo "  PE             treat data as pair-end reads (default: true)"


# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Generate fastqc report for a set of reads
fastqc:
	mkdir -p output/fastqc/
	fastqc -o output/fastqc/ --threads $(THREADS) $(READS)

# Consolidate fastqc files into a single report
multiqc: output/fastqc/
	mkdir -p output/multiqc/
	multiqc -o output/multiqc/ output/fastqc/

fastp: $(wildcard $(READ_DIR)/*.fastqc*)
	mkdir -p output/fastp
	mkdir -p output/fastp/reads
	mkdir -p output/fastp/reports
ifeq ($(PE),true)
	fastp -i $(R1) -I $(R2)	$(fastp_opts) -o output/fastp/reads/trimmed_$(shell basename $(R1)) -O output/fastp/reads/trimmed_$(shell basename $(R2))
else
	fastp -i $(R)	$(fastp_opts) -o output/fastp/reads/trimmed_$(shell basename $(R))
endif
