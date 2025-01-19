# 
# Retrieve read data from the Sequence Read Archive (SRA).
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/fetch.mk

.PHONY: help params init clean

# ncbi-datasets parameters
ncbi_datasets_opts :=
ifeq ($(INCLUDE_GFF),true)
ncbi_datasets_opts += --include genome,gff3
endif

%.runinfo.csv:
	@. $(MK_ROOT)/fetch_runinfo.sh $(PRJNA)


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
fetch_reads_opts := $(if $(filter true, $(PE)), -p)
fetch_reads_opts += $(if $(X), -x $(X))


reads:
# Fetch project metadata and extract a list of SRR accessions
ifdef PRJNA
	@. $(MK_ROOT)/fetch_runinfo.sh $(PRJNA) | cut -d, -f1 | tail -n +2 \
		| parallel ". $(MK_ROOT)/fetch_reads.sh $(fetch_reads_opts) -a {}"
else ifdef SRR
	@. $(MK_ROOT)/fetch_reads.sh $(fetch_reads_opts) -a $(SRR)
else
	@echo "Error: please provide either a project (PRJNA) or read (SRR) accession"
endif


# Retrieve sequence files from NCBI
seq:
	@. $(MK_ROOT)/fetch_sequence.sh -a $(ACC) $(if $(filter true, $(WITH_GFF)), -g)

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
