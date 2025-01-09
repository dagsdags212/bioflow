# 
# Retrieve read data from the Sequence Read Archive (SRA).
#

SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
.PHONY: help init sra ref pdb

# Parameters for sequencing reads
PRJNA ?=
SRR ?=
X ?=
PE ?= false

# Parameters for reference genomes
ACC ?=
INCLUDE_GFF ?= false

# Parameters for PDB files
PDB ?= 

# Environment manager
ENV_MANAGER ?= micromamba

# Conda environment
ENV ?= fetch

# Check if dependencies are installed
dependencies := entrez-direct sra-tools ncbi-datasets-cli

# fastq-dump parameters
fastq_dump_opts := --origfmt
ifeq ($(PE),true)
fastq_dump_opts += --split-3
endif

# ncbi-datasets parameters
ifeq ($(INCLUDE_GFF),true)
	ncbi_datasets_opts := --include genome,gff3
endif

# Display help message
help:
	@echo ""
	@echo "fetch.mk: retrieve biological data of all forms"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/fetch.mk <command> [options]"
	@echo ""
	@echo "COMMANDS:"
	@echo "  sra - retrieve sequencing data from the SRA"
	@echo "  ref - retrieve genome and annotation from NCBI"
	@echo "  pdb - retrieve structure data from PDB"

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies)

# Retrive sequencing reads from the SRA
sra:
	mkdir -p reads
ifdef PRJNA
	# Fetch project metadata and extract a list of SRR accessions
	@echo "Fetching runinfo for $(PRJNA)"
	esearch -db sra -query $(PRJNA) | efetch -format runinfo \
		| cut -d, -f 1 | tail -n +2 > /tmp/fetch_acc.txt
	# Create directoy for each SRR accession
	cat /tmp/fetch_acc.txt | parallel -- mkdir -p reads/{}
ifdef X
	# Limit number of spots to X
	@echo "Downloading $(X) spots for each sample"
	cat /tmp/fetch_acc.txt | parallel -j 3 -- 'fastq-dump -O reads/{} -X $(X) $(fastq_dump_opts) {}'
else
	# Download all spots
	@echo "Downloading all spots for each sample"
	cat /tmp/fetch_acc.txt | parallel -j 3 -- 'fastq-dump -O reads/{} $(fastq_dump_opts) {}'
endif
endif

ifdef SRR
	mkdir -p reads/$(SRR)
ifdef X
	@echo "Downloading $(X) spots for $(SRR)"
	fastq-dump -O reads/$(SRR) -X $(X) $(fastq_dump_opts) $(SRR)
else
	@echo "Downloading all spots for $(SRR)"
	fastq-dump -O reads/$(SRR) $(fastq_dump_opts) $(SRR)
endif
endif

# Retrieve reference genomes from NCBI
ref:
	mkdir -p ref
	datasets download genome accession $(ACC) $(ncbi_datasets_opts) --filename "$(ACC).zip"
	# Decompress and extract relevant files, delete the rest
	unzip $(ACC).zip
	find ncbi_dataset/ -name "*.fna" -exec mv -t ref/ {} +
	find ncbi_dataset/ -name "*.gff" -exec mv -t ref/ {} +
	rm -rf ncbi_dataset/ $(ACC).zip *.txt *.md

# Retrieve structures files from PDB
pdb:
	mkdir -p pdb
	@echo "Fetching $(PDB) from PDB"
	pdb_fetch $(PDB) > pdb/$(PDB).pdb

# TODO: retrieve annotation data from genbank
