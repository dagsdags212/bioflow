# 
# Retrieve read data from the Sequence Read Archive (SRA).
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# Absoulate path to conda environment files
ENV_ROOT := ${BIOFLOW_PREFIX}/envs

# Conda environment for module
ENV := bf-fetch

# Execute within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/fetch.mk

.PHONY: help params init clean

# Display help message
help:
	@echo
	@echo "bf-fetch: a module for retrieving biological data"
	@echo
	@echo "Usage:"
	@echo "  bf-fetch <command> [options]"
	@echo
	@echo "Commands:"
	@echo "  reads   - retrieve sequencing data from the SRA"
	@echo "  seq     - retrieve genomic sequence data from NCBI"
	@echo "  pdb     - retrieve protein structures from PDB"
	@echo "  pubmed  - query PubMed for a list of related journals"
	@echo


# Display available parameters
params:
	@echo
	@echo "For retrieving SEQUENCE RECORDS"
	@echo "  ACC               sequence accession identifier"
	@echo "  INCLUDE_GFF       download annotation file for the accession (if available)"
	@echo "  INCLUDE_GB        download GenBank record for the accession"
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
	$(ENV_MANAGER) create -f $(ENV_ROOT)/$(ENV).yml

# Retrive sequencing reads from the SRA
fetch_reads_opts := $(if $(filter true, $(PE)), -p)
fetch_reads_opts += $(if $(X), -x $(X))

reads:
# Fetch project metadata and extract a list of SRR accessions
ifdef PRJNA
	@. $(MK_ROOT)/fetch_runinfo.sh $(PRJNA) | cut -d, -f1 | tail -n +2 \
		| parallel ". $(MK_ROOT)/fetch_reads.sh {} $(fetch_reads_opts)"
else ifdef SRR
	@. $(MK_ROOT)/fetch_reads.sh $(SRR) $(fetch_reads_opts)
else
	@echo "Error: please provide either a project (PRJNA) or read (SRR) accession"
endif

fetch_sequence_opts := $(if $(filter true,$(WITH_GFF)),-g)
fetch_sequence_opts += $(if $(filter true,$(WITH_GB)),-b)

# Retrieve sequence files and records from NCBI
seq:
ifdef ACC
	@. $(MK_ROOT)/fetch_sequence.sh $(ACC) $(fetch_sequence_opts)
else
	@echo "Error: sequence accession not provided"
endif

# Retrieve structures files from PDB
pdb:
ifdef PDB
	@. $(MK_ROOT)/fetch_structure.sh $(PDB)
else
	@echo "Error: protein accession not provided"
endif

# Query Pubmed for a list of jounrnals base on a search string
pubmed:
ifdef QUERY
	@. $(MK_ROOT)/query_pubmed.sh '$(QUERY)'
else
	@echo "Error: query string not provided"
endif

clean:
	rm -rf reads/ data/
