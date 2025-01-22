.RECIPEPREFIX = >

# Project accession.
PRJN ?=

# Sequencing run accession.
SRX ?=

# Number of spots to download.
X ?=

# Enable pair-end mode.
PE ?= false

# Save reads in this directory.
OUTDIR ?= data/reads

# ===== fastq-dump options =====
fastq_dump_opts := --origfmt -v
fastq_dump_opts += $(if $(X),-X $(X))
ifeq ($(PE),true)
fastq_dump_opts += --split-3
endif
fastq_dump_opts += -O $(OUTDIR)
# ===== fastq-dump options =====

# Target for storing runinfo.
RUNINFO = $(PRJN)_runinfo.csv

# Targt for saving accession list.
ACC_LIST = $(PRJN)_accessions.txt

reads:
> @mkdir -p $(OUTDIR)
ifdef PRJN
# Fetch metadata associated with accession.
> @echo "=== Project accession detected: $(PRJN) ==="
> @echo "=== Fetching runinfo for $(PRJN) ==="
> esearch -db sra -query $(PRJN) | efetch -format runinfo > $(RUNINFO)
> @echo "=== Extracting accessions for individual runs ==="
> cat $(RUNINFO) | cut -d, -f1 | tail -n +2 > $(ACC_LIST)
> @echo "=== Downloading reads ==="
> cat $(ACC_LIST) | parallel "fastq-dump $(fastq_dump_opts) {}"

# Downloading single sequencing run.
else ifdef SRX
> @echo "=== Downloading reads for $(SRX) ==="
> fastq-dump $(fastq_dump_opts) $(SRX)

# No accession provided, do nothing.
else
> @echo "Error: accession not provided"
endif
