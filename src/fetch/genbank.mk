#
# Download sequence data from NCBI.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Micromamba environment.
ENV = bf-fetch-genbank

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# NCBI accession number.
ACC ?= AF086833

# The accession as a fasta file.
REF ?= refs/$(ACC).fa

# The accession as a genbank file.
GBK ?= refs/$(ACC).gb

# The accession as a annotation file.
GFF ?= refs/$(ACC).gff


# Display usage.
help:
	@echo ""
	@echo "genbank.mk: download sequence data from GenBank"
	@echo ""
	@echo "Usage:"
	@echo "  bf-genbank [options] ACC=<ACC>"
	@echo ""
	@echo "Commands:"
	@echo "  fasta          retrieve a sequence file in FASTA format"
	@echo "  gff            retrieve an annotation file in GFF3 format"
	@echo "  genbank        retrieve a genbank record"
	@echo "  all            retrieve all supported records associatiated"
	@echo "                 with a given accession"
	@echo ""
	@echo "Options:"
	@echo "  ACC            an accession identifier (required)"
	@echo "  REF            a path for storing the FASTA file (optional)"
	@echo "  GBK            a path for storing the GenBank record (optional)"
	@echo "  GFF            a path for storing the annotation file (optional)"
	@echo ""

example:
	@echo ""
	@echo "# Download the sequence and annotation files of Ebola"
	@echo "make -f genbank.mk ACC=AF086833 fasta gff"
	@echo ""
	@echo "# Specify the output path"
	@echo "make -f genbank.mk ACC=AF086833 GBK=records/ebola.gb genbank"
	@echo ""
	@echo "# Download all records associated with an accession"
	@echo "make -f genbank.mk ACC=AF086833 all"
	@echo ""

# Obtain a fasta file from NCBI.
$(REF):
	mkdir -p $(dir $@)
	$(ENV_RUN) efetch -db nuccore -id $(ACC) -format fasta > $@

# Obtain a genbank records from NCBI.
$(GBK):
	mkdir -p $(dir $@)
	$(ENV_RUN) efetch -db nuccore -id $(ACC) -format gb > $@

# Obtain an annotation file from NCBI.
$(GFF):
	mkdir -p $(dir $@)
	$(ENV_RUN) bio fetch $(ACC) --format gff > $@

# Download fasta file.
fasta:: $(REF)
	@ls -lh $(REF)

# Download genbank file.
genbank:: $(GBK)
	@ls -lh $(GBK)

# Download annotation file.
gff:: $(GFF)
	@ls -lh $(GFF)

# Remove fasta file.
fasta!::
	rm -rf $(REF)

# Remove genbank file.
genbank!::
	rm -rf $(GBK)

# Remove annotation file.
gff!::
	rm -rf $(GFF)

# Set run target to reference file.
run:: fasta

# Remove the run target.
run!:: fasta!

# Generate all three outputs.
all:: fasta genbank gff

# Remove downloaded files.
clean:
	rm -f $(REF) $(GBK) $(GFF)

# Run the test suite.
test:
	make -f $(ROOT_PATH)/genbank.mk fasta! fasta ACC=$(ACC) REF=$(REF)
	make -f $(ROOT_PATH)/genbank.mk gff! gff ACC=$(ACC) GFF=$(GFF)
	make -f $(ROOT_PATH)/genbank.mk genbank! genbank ACC=$(ACC) GBK=$(GBK)

DEPS := entrez-direct

# Display dependencies.
install::
	micromamba create -n ${ENV}
	${ENV_RUN} micromamba install ${DEPS}
	${ENV_RUN} pip install bio --upgrade

# Non-file targets.
.PHONY: help example run run! test install
