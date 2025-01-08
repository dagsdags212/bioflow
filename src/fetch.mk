# 
# Retrieve read data from the Sequence Read Archive (SRA).
#
# Parameters:
# 	PRJNA - Sequencing project ID
# 	SRR - Sequencing run accession
# 	X - number of spots
#
# Usage:
# 	# download reads from all samples, specifying number of spots
# 	make -f src/fetch.mk PRJNA=PRJNA1066786 X=100000
#
# 	# download reads (all spots) from a single sequencing run
# 	make -f src/fetch.mk SRR=SRR27644850

.DELETE_ON_ERROR:
.ONESHELL:
.PHONY: sra

# Conda environment
ENV = fetch

# Check if dependencies are installed
dependencies := entrez-direct sra-tools

PRJNA ?=
SRR ?=
X ?=

# fastq-dump parameters
fastq_dump_opts = --split-3 --origfmt

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
