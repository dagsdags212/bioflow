#
# Assemble long reads into contigs using flye.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Absolute path of parent directory.
ROOT_PATH = ${shell dirname ${abspath ${firstword ${MAKEFILE_LIST}}}}

# Micromamba environment.
ENV = bf-assembly

# Run command within environment.
ENV_RUN = micromamba run -n ${ENV}

# Number of worker threads.
THREADS ?= 4

# Long read file.
R1 ?=

# Output directory.
OUTDIR ?= flye_out

# Target output.
CONTIGS = ${OUTDIR}/assembly.fasta

# Estimated genome size.
GENOME_SIZE ?= 0.2m

# Supported long-read sequencing platforms.
SUPPORTED_PLATFORMS := pacbio ont

# Selected platform.
PLATFORM ?= ont

# Number of polishing iterations.
ITER ?= 5

# Megahit options.
FLAGS := -t ${THREADS} -i ${ITER} -g ${GENOME_SIZE}

# Conditional flags.
ifeq (${PLATFORM},pacbio)
	FLAGS += --pacbio-raw ${R1}
else ifeq (${PLATFORM},ont)
	FLAGS += --nano-raw ${R1}
else
	@echo "Error: platform not specified or invalid, choose [pacbio|ont]"
	@exit -1
endif

# Compose flye command.
FLYE_CMD = flye ${FLAGS} -o ${OUTDIR}

help::
	@echo "flye.mk: assemble long reads into contigs using flye"
	@echo ""
	@echo "Input:"
	@echo "  R1=${R1}"
	@echo "  OUTDIR=${OUTDIR}"
	@echo ""
	@echo "  PLATFORM=${PLATFORM}"
	@echo "  ITERATIONS=${ITER}"
	@echo "  THREADS=${THREADS}"
	@echo ""
	@echo "Commands:"
	@echo "  run        perform sequence assembly"
	@echo "  install    print command for installing dependencies"
	@echo "  clean      remove all flye output"
	@echo ""
	@echo "Invoking 'run' will execute the following command:"
	@echo ""
	@echo "  ${FLYE_CMD}"
	@echo ""

# R1 must exist.
${R1}:
	@if [ ! -f ${R1} ]; then
		echo "Error: R1 not found {R1=${R1}}"
	fi
	@exit -1

${CONTIGS}: ${R1}
	# Remove output directory if it already exists.
	[ -z ${dir $@} ] && rm -rf ${dir $@}

	# Create output directory.
	mkdir -p $(dir $@)

	# Perform sequence assembly with flye.
	${ENV_RUN} ${FLYE_CMD}

# Invoked megahit command.
run: ${CONTIGS}
	ls -lh ${dir $^}

stats:
	seqkit stats ${CONTIGS}
	
run!::
	rm -rf ${OUTDIR}

# Alternative target to run!
clean: run!

# Print command for installing dependencies.
install::
	@echo "micromamba install flye seqkit"

.PHONY: help run stats run! clean install
