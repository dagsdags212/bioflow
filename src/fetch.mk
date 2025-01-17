# 
# Retrieve read data from the Sequence Read Archive (SRA).
#

# Absoluate path for the `bioflow/src` directory
SRC_ROOT := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

# Absolute path for `bioflow` directory
PROJECT_ROOT := $(shell dirname $(SRC_ROOT))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := $(PROJECT_ROOT)/config

# import config variables
include $(CONFIG_ROOT)/_config.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

.PHONY: help params init clean

# Conda environment
ENV := bf-fetch

# Path to conda environment
ENV_DIR = $(shell $(ENV_MANAGER) info | grep "envs directories" | cut -d ":" -f 2 | xargs)/$(ENV)

# Check if dependencies are installed
dependencies := entrez-direct sra-tools ncbi-datasets-cli

# Parameters for sequencing reads
PRJNA ?=
SRR ?=
X ?=
PE ?= false

# Convert MGNify ID to SRA accession
MGYS = $(if $(findstring MGYS,$(PRJNA)),$(shell . $(PROJECT_ROOT)/scripts/mgnify2sra.sh $(PRJNA)))

# Parameters for reference genomes
ACC ?=
INCLUDE_GFF ?= false

# Parameters for PDB files
PDB ?= 

# Query string of pubmed search
QUERY ?=

# fastq-dump parameters
fastq_dump_opts := --origfmt
ifeq ($(PE),true)
fastq_dump_opts += --split-3
endif

# ncbi-datasets parameters
ncbi_datasets_opts :=
ifeq ($(INCLUDE_GFF),true)
ncbi_datasets_opts += --include genome,gff3
endif

# Display help message
help:
	@echo
	@echo "fetch.mk: retrieve a variety of biological file formats"
	@echo
	@echo "Usage:"
	@echo "  bf-fetch <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  genbank - retrive a GenBank record from NCBI"
	@echo "  pdb     - retrieve structure data from PDB"
	@echo "  pubmed  - retrieve a list of Pubmed journals"
	@echo "  sra     - retrieve sequencing data from the SRA"
	@echo "  ref     - retrieve genome and annotation from NCBI"
	@echo

# Display available parameters
params:
	@echo
	@echo "For retrieving SEQUENCE RECORDS"
	@echo "  ACC               sequence accession identifier"
	@echo "  INCLUDE_GFF       download annotation file for fetched reference genome"
	@echo
	@echo "For retrieving READS"
	@echo "  PRJNA             sequencing PROJECT identifier"
	@echo "  SRR               sequencing RUN identifier"
	@echo "  PE                if true, download reads in pair-end mode (default: false)"
	@echo "  X                 number of spots to download"
	@echo
	@echo "For retrieving PROTEIN STRUCTURES"
	@echo "  PDB               4-character protein identifier"
	@echo
	@echo "For retrieving JOURNAL METADATA"
	@echo "  QUERY             search string used to query Pubmed"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bf-fetch)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -n $(ENV) $(dependencies) --yes

# Retrive sequencing reads from the SRA
sra:
	@mkdir -p reads/
ifdef PRJNA
	# Fetch project metadata and extract a list of SRR accessions
	@echo "Fetching runinfo for $(PRJNA)"
ifdef MGYS
	@echo "Metagenomic samples detected"
	@echo "Converting MGNify accesion to SRA accession"
endif
	@esearch -db sra -query $(if $(MGYS),$(MGYS),$(PRJNA)) | efetch -format runinfo \
		| cut -d, -f 1 | tail -n +2 > /tmp/fetch_acc.txt
	# Create directoy for each SRR accession
	@cat /tmp/fetch_acc.txt | parallel -- mkdir -p reads/{}
ifdef X
	# Limit number of spots to X
	@echo "Downloading $(X) spots for each sample"
	@cat /tmp/fetch_acc.txt | parallel -j 3 -- 'fastq-dump -O reads/{} -X $(X) $(fastq_dump_opts) {}'
else
	# Download all spots
	@echo "Downloading all spots for each sample"
	@cat /tmp/fetch_acc.txt | parallel -j 3 -- 'fastq-dump -O reads/{} $(fastq_dump_opts) {}'
endif
endif

ifdef SRR
	@mkdir -p reads/$(SRR)
ifdef X
	@echo "Downloading $(X) spots for $(SRR)"
	@fastq-dump -O reads/$(SRR) -X $(X) $(fastq_dump_opts) $(SRR)
else
	@echo "Downloading all spots for $(SRR)"
	@fastq-dump -O reads/$(SRR) $(fastq_dump_opts) $(SRR)
endif
endif

# Retrieve reference genomes from NCBI
ref:
	@mkdir -p ref/
	datasets download genome accession $(ACC) $(ncbi_datasets_opts) --filename "$(ACC).zip"
	# Decompress and extract relevant files, delete the rest
	unzip $(ACC).zip
	find ncbi_dataset/ -name "*.fna" -exec mv -t ref/ {} +
	find ncbi_dataset/ -name "*.gff" -exec mv -t ref/ {} +
	rm -rf ncbi_dataset/ $(ACC).zip *.txt *.md

# Retrieve a genbank file from NCBI
genbank:
ifdef ACC
	@printf "Fetching GenBank record for $(ACC)\n" 1>&2
	@efetch -db nuccore -id ${ACC} -format gb
else
	@printf "Error: accession not provided\n" 1>&2
endif

# Retrieve structures files from PDB
pdb:
ifdef PDB
	@printf "Fetching structure file for $(PDB)\n" 1>&2
	pdb_fetch $(PDB)
else
	@printf "Error: PDB ID not provided" 1>&2
endif

# Query Pubmed for a list of jounrnals base on a search string
pubmed:
ifdef QUERY
	@. $(PROJECT_ROOT)/scripts/query_pubmed.sh '$(QUERY)'
else
	@printf "Error: query string not provided" 1?	
endif

clean:
	rm -rf reads/ ref/
	rm /tmp/fetch_acc.txt
