# 
# Perform quality checks and filtering on sequencing data.
#
# Parameters:
# 	Global settings
# 		READ_DIR - directory path containing read data
# 		THREADS  - number of cores
#
# 	Environment settings
# 		ENV         - environment name
# 		ENV_MANAGER - environment manager
#
# 	Trimming and filtering
# 		MINLEN  - minimum read length
# 		MAXLEN  - maximum read length
# 		MINQUAL - minimum acceptable quality score
# 		PE      - if true, treat data as pair-end reads
#

.DELETE_ON_ERROR:
.ONESHELL:
.PHONY: help init fastqc fastp

READ_DIR ?=
THREADS ?= 4

# fastqc parameters
READS = $(shell find $(READ_DIR) -type f -name "*fastq*")

# fastp parameters
MINLEN ?= 30
MAXLEN ?= 150
MINQUAL ?= 20
PE ?= true

ifeq ($(PE),true)
R1 = $(shell find $(READ_DIR) -type f -name "*_1*fastq*")
R2 = $(shell find $(READ_DIR) -type f -name "*_2*fastq*")
else
R = $(shell find $(READ_DIR) -type f -name "*fastq*")
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
	@echo "  THREADS        number of cores"
	@echo "Environment settings"
	@echo "  ENV            environment name"
	@echo "  ENV_MANAGER    environment manager"
	@echo "Trimming and filtering"
	@echo "  MINLEN         minimum read length"
	@echo "  MAXLEN         maximum read length"
	@echo "  MINQUAL        minimum acceptable quality score"
	@echo "  PE             if true, treat data as pair-end reads"


# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Generate fastqc report for a set of reads
fastqc:
	mkdir -p fastqc/
	fastqc -o fastqc/ --threads $(THREADS) $(READS)

# Consolidate fastqc files into a single report
multiqc: fastqc/
	mkdir -p multiqc/
	multiqc -o multiqc/ fastqc/

fastp: $(READ_DIR)
	mkdir -p fastp
	mkdir -p fastp/reads
	mkdir -p fastp/reports
ifeq ($(PE),true)
	fastp -i $(R1) -I $(R2)	$(fastp_opts) -o fastp/reads/trimmed_$(shell basename $(R1)) -O fastp/reads/trimmed_$(shell basename $(R2))
else
	fastp -i $(R)	$(fastp_opts) -o fastp/reads/trimmed_$(shell basename $(R))
endif
