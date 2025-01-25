# 
# Annotate prokaryotic genes with Prokka.
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
ENV = bf-annotation

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 2

# Supported formats of prodigal output.
FORMATS ?= gbk gff sco

# Format of prodigal output.
FORMAT ?= gbk

# FASTA file input.
IN ?= 

# Basename for input file.
BASE = $(basename $(notdir $(IN)))

# Output directory for prokka output.
OUTDIR ?= prokka

# Output annotation file.
OUT ?= $(OUTDIR)/$(BASE).gff

# E-value threshold.
CUTOFF ?= 0.000001

# Annotation mode.
KINGDOM ?= bacteria

# Specify taxonomic details of source organism.
GENUS ?=
SPECIES ?=
STRAIN ?=
PLASMID ?=

# Prokka flags.
FLAGS := --cpus $(THREADS) --prefix $(BASE) --gffver 3 --evalue $(CUTOFF) --force
FLAGS += $(if $(KINGDOM),--kingdom $(KINGDOM))
FLAGS += $(if $(GENUS),--genus $(GENUS))
FLAGS += $(if $(SPECIES),--species $(SPECIES))
FLAGS += $(if $(STRAIN),--strain $(STRAIN))
FLAGS += $(if $(PLASMID),--plasmid $(PLASMID))

# Compose prokka command.
PROKKA_CMD = prokka $(FLAGS) --outdir $(OUTDIR) $(IN) 


help::
	@echo "prokka.mk: annotate prokaryotic genomes with Prokka"
	@echo ""
	@echo "Input:"
	@echo "  IN=$(IN)"
	@echo "  OUTDIR=$(OUTDIR)"
	@echo "  FORMAT=$(FORMAT)"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform genome annotation of input file"
	@echo "  list       display available Prokka databases"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all Prokka output"
	@echo "  cite       print tool citation"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  $(PROKKA_CMD)"
	@echo ""

# Input file must exist.
$(IN):
	@if [ ! -f $(IN) ]; then
		@echo "# Error: FASTA file not found (IN=$(IN))."
	fi

$(OUT): $(IN)
	# Create output directory.
	mkdir -p $(dir $@)

	# Perform gene prediction.
	$(ENV_RUN) $(PROKKA_CMD)

OUT_EXT := gff gbk fna faa ffn sqn fsa tbl err log txt tsv
OUT_FILES = $(addprefix $(OUTDIR)/$(BASE).,$(OUT_EXT))

# Invoked gene prediction command.
run: $(OUT_FILES)
	ls -lh $^

# Remove prodigal output.
run!:: $(OUTFILES)
	rm -f $^

# List prokka databases.
list::
	@$(ENV_RUN) prokka --listdb

cite::
	@$(ENV_RUN) prokka --citation

# Alias for run!
clean: run!

# Display installation command.
install::
	@echo "micromamba install prokka"

# Rules with no target files.
.PHONY: help run run! list clean cite install
