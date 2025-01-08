# 
# Retrieve read data from the Sequence Read Archive (SRA).
#
# Parameters:
# 	ENV_MANAGER - environment manager installed in system
# 	ENV - environment name
# 	PRJNA - Sequencing project ID
# 	SRR - Sequencing run accession
# 	X - number of spots
#
# Usage:
#   # print usage message
#   make -f src/fetch.mk help
#
# 	# create new environment named 'fetch' and install dependencies using conda
#   make -f src/fetch.mk init ENV_MANAGER=conda ENV=fetch
#
# 	# download reads from all samples, specifying number of spots
# 	make -f src/fetch.mk sra PRJNA=PRJNA1066786 X=100000
#
# 	# download reads (all spots) from a single sequencing run
# 	make -f src/fetch.mk sra SRR=SRR27644850

.DELETE_ON_ERROR:
.ONESHELL:
.PHONY: help init sra ref pdb

# Parameters for sequencing reads
PRJNA ?=
SRR ?=
X ?=

# Parameters for reference genomes
ACC ?=

# Parameters for PDB files
PDB ?= 

# Environment manager
ENV_MANAGER ?= micromamba

# Conda environment
ENV ?= fetch

# Check if dependencies are installed
dependencies := entrez-direct sra-tools ncbi-datasets-cli

# fastq-dump parameters
fastq_dump_opts = --split-3 --origfmt

# ncbi-datasets parameters
ncbi_datasets_opts = --include genome,gff3

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
	datasets download genome accession $(ACC) $(ncbi_datasets_opts) --filename $(ACC)
	# Decompress and extract relevant files, delete the rest
	unzip $(ACC).zip
	mv $(ACC)/data/**/*.{fna,faa,gff3} ref/
	rm -rf $(ACC) $(ACC).zip

# Retrieve structures files from PDB
pdb:
	mkdir -p pdb
	@echo "Fetching $(PDB) from PDB"
	pdb_fetch $(PDB) > pdb/$(PDB).pdb

# TODO: retrieve annotation data from genbank
