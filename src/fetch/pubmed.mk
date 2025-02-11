#
# Perform queries against the PubMed literature database.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Micromamba environment.
ENV = bf-fetch-pubmed

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Pubmed identifier.
PMID ?= 22388286 19451168


help:
	@echo ""
	@echo "pubmed.mk: query the PubMed database"
	@echo ""
	@echo "Usage:"
	@echo "  bf-pubmed [options] PMID=<PMID>"
	@echo ""
	@echo "Commands:"
	@echo "  abstract         display journal abstract"
	@echo "  apa              display APA citation"
	@echo "  medline          retrieve medline metadata"
	@echo "  citations        retrieve a list of cited articles"
	@echo "  install          initialize conda environment"
	@echo ""
	@echo "Options:"
	@echo "  PMID             a PubMed accession"
	@echo ""

# Retrieve abstract.
abstract:: $(PMID)
	@$(ENV_RUN) efetch -db pubmed -id $(PMID) -format abstract

# Retrieve APA citation.
apa:: $(PMID)
	@$(ENV_RUN) efetch -db pubmed -id $(PMID) -format apa | cut -f2

# Retrieve medline record.
medline:: $(PMID)
	@$(ENV_RUN) efetch -db pubmed -id $(PMID) -format medline

# Retrieve a list of citations.
citations:: $(PMID)
	@$(ENV_RUN) efetch -db pubmed -id $(PMID) -format xml | \
		$(ENV_RUN) xtract -pattern PubmedArticle -element PubmedData/ReferenceList \
		-group Reference -deq "\n- " -element Citation

# Run the test suite.
test:
	make -f $(ROOT_PATH)/pubmed.mk PMID=$(PMID) abstract
	make -f $(ROOT_PATH)/pubmed.mk PMID=$(PMID) apa
	make -f $(ROOT_PATH)/pubmed.mk PMID=$(PMID) medline
	make -f $(ROOT_PATH)/pubmed.mk PMID=$(PMID) citations

DEPS := entrez-direct
# Initialize micromamba environment.
install::
	micromamba create -n ${ENV}
	${ENV_RUN} micromamba install ${DEPS}

# Non-file targets.
.PHONY: help $(PMID) abstract apa medline citations
