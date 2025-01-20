# 
# Perform gene prediction to annotate genomes.
#

MK_ROOT = $(shell dirname $(abspath $(firstword $(MAKEFILE_LIST))))

# Absoluate path for the `bioflow/config` directory
CONFIG_ROOT := ${BIOFLOW_PREFIX}/config

# Absoulate path to conda environment files
ENV_ROOT := ${BIOFLOW_PREFIX}/envs

# Conda environment for module
ENV := bf-annotation

# Execute within the module environment
RUN = $(ENV_MANAGER) run -n $(ENV)

# import Make-specific configuration
include $(CONFIG_ROOT)/_preamble.mk

# import global variables
include $(CONFIG_ROOT)/_globals.mk

# import module configuration
include $(CONFIG_ROOT)/annotation.mk

.PHONY: help params init clean

# Display help message
help:
	@echo
	@echo "annotation.mk: predict genes and provide annotation for genomes"
	@echo ""
	@echo "Usage:"
	@echo "  make -f src/annotation.mk <command> [options]"
	@echo
	@echo "COMMANDS:"
	@echo "  annotate   - run annotation pipeline on sequence file"
	@echo "  datasets   - list available BUSCO datasets"
	@echo "  predict    - run de novo gene prediction using Prodigal"
	@echo

# Create new self-contained environment
init:
	$(ENV_MANAGER) create -f $(ENV_ROOT)/annotation.yml

# Display available parameters
params:
	@echo
	@echo "Annotation settings"
	@echo "  FA                path to target FASTA file"
	@echo "  MODE              BUSCO analysis mode to run [genome|protein|transcriptome] (default: genome)"
	@echo "  DOMAIN            specify the domain where the input sequence belongs to [prokaryote|eukaryote]"
	@echo "  PRODIGAL_OUTNAME  filename for Prodigal output"
	@echo "  PRODIGAL_OUTFMT   specify output format for Prodigal [gff|gbk|sco] (default: gff)"
	@echo "  PROKKA_OUTDIR     directory for saving Prokka output"
	@echo
	@echo "Global settings"
	@echo "  THREADS           number of cores (default: 4)"
	@echo
	@echo "Environment settings"
	@echo "  ENV               environment name (default: bwf-annotation)"
	@echo "  ENV_MANAGER       environment manager (default: micromamba)"
	@echo

# List available BUSCO datasets
datasets:
	@busco --list-datasets

run_prodigal_opts := $(if $(filter gbk gff sco,$(PRODIGAL_OUTFMT)),-f $(PRODIGAL_OUTFMT))
run_prodigal_opts += $(if $(PRODIGAL_OUTNAME),-o $(PRODIGAL_OUTFMT))

# Run ab initio gene prediction using prodigal
predict:
ifdef FA
	@. $(MK_ROOT)/run_prodigal.sh $(run_prodigal_opts) $(FA)
	@prodigal -i $(FA) $(prodigal_opts)
else
	@echo "Error: input FASTA file not provided"
endif

# BUSCO flags
run_busco_opts := $(if $(filter genome transcriptome proteins,$(MODE)),-m $(MODE))
run_busco_opts += $(if $(filter prokaryote eukaryote,$(DOMAIN)),-d $(DOMAIN))
run_busco_opts += $(if $(PRODIGAL_OUTDIR),-o $(PRODIGAL_OUTDIR))

# Prokka flags
run_prokka_opts := $(if $(PROKKA_OUTDIR),-o $(PROKKA_OUTDIR))

# Run BUSCO annotation pipeline
annotate:
ifdef FA
	@# annotate using BUSCO
	@. $(MK_ROOT)/run_busco.sh $(FA) $(run_busco_opts)
	@# annotate using Prokka
	@. $(MK_ROOT)/run_prokka.sh $(run_prokka_opts) $(FA)
else
	@echo "Error: input FASTA file not provided"
endif

clean:
	rm -rf output/busco/ output/prodigal/ output/prokka/
	rm -rf busco_downloads/
