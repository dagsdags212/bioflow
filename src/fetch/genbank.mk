#
# Download sequence data from NCBI.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Absolute path of parent directory.
ROOT_PATH = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Micromamba environment.
ENV = bf-fetch

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

help:
	@echo "#"
	@echo "# genbank.mk: download sequence data from GenBank"
	@echo "#"
	@echo "# ACC=$(ACC)"
	@echo "# REF=$(REF)"
	@echo "# GBK=$(GBK)"
	@echo "# GFF=$(GFF)"
	@echo "#"
	@echo "# make fasta|gff|genbank|all"
	@echo "#"

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

# Display dependencies.
install::
	@echo micromamba install entrez-direct
	@echo pip install bio --upgrade

# Non-file targets.
.PHONY: help run run! test install
