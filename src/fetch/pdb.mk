#
# Retrieve protein structures from the Protein Data Bank.
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

# PDB identifier.
ID ?= 7FCD

# Output diretory.
DIR ?= data

# Target path to PDB file.
PDB = $(DIR)/$(ID).pdb

help:
	@echo "#"
	@echo "# pdb.mk: download protein structures from the PDB"
	@echo "#"
	@echo "# ID=$(ID)"
	@echo "#"
	@echo "# make run|test|clean|install"
	@echo "#"

$(PDB):
	# Create output directory.
	mkdir -p $(DIR)

	# Dowload protein structure.
	$(ENV_RUN) pdb_fetch $(ID) > $(PDB)

run: $(PDB)
	@ls -lh $(PDB)

# Run the test suite.
test:
	make -f $(ROOT_PATH)/pdb.mk ID=$(ID) run

# Delete PDB file.
run!: $(PDB)
	rm -f $(PDB)

# Alternative rule for run!
clean: run!

# Display dependencies.
install::
	@echo micromamba install pdb-tools

# Non-file targets.
.PHONY: help run test run! clean install
