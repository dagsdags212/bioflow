#
# Map reads against a reference using minimap2.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Micromamba environment.
ENV = bf-mapping-minimap2

# Run command within environment.
ENV_RUN = micromamba run -n ${ENV}

# Number of worker threads.
THREADS ?= 4

# FASTQ file to align.
R1 ?=

# Root name of FASTA file.
BASENAME = $(notdir $(basename ${R1}))

# The reference genome.
REF ?= refs/${ACC}.fa

# Output directory for BAM files.
DIR ?= bam

# The alignment file.
BAM ?= ${DIR}/${BASENAME}.bam


# Print usage.
help::
	@echo ""
	@echo "minimap2.mk: map reads against a reference using minimap2"
	@echo ""
	@echo "Usage:"
	@echo "  make -f minimap2.mk [options] <command>"
	@echo ""
	@echo "Commands:"
	@echo "  align          map reads against a reference using minimap2"
	@echo "  run            alias for 'align'"
	@echo "  install        initialize conda environment"
	@echo "  clean          remove all output generated by this program"
	@echo ""
	@echo "Options:"
	@echo "  REF            a reference file"
	@echo "  R1             read file in FASTQ format"
	@echo "  DIR            a directory path for storing BAM files"
	@echo "  THREADS        number of worker threads (default: 4)"
	@echo ""

# Read 1 must exist.
${R1}:
	@echo "# Error: FASTQ file not found {R1=${R1}}"
	@exit -1

# Reference file must exist.
${REF}:
	@if [ ! -f ${REF} ]; then
		echo "# Error: Reference file not found {REF=${REF}}"
		exit -1
	fi

# Minimap2 options.
FLAGS := -a -t ${THREADS}

# Generate a sorted alignment file.
${BAM}: ${REF} ${R1}
	# Output directory for BAM files.
	mkdir -p $(dir $@)

	# Perform read mapping.
	${ENV_RUN} minimap2 ${FLAGS} ${REF} ${R1} | samtools view -b | samtools sort -@ ${THREADS} > $@

# Create the BAM index.
${BAM}.bai: ${BAM}
	samtools index ${BAM}

# Invoke the read mapping rule.
align: ${BAM} ${BAM}.bai
	@ls -lh $^

# Alternative rule for align.
run: align

# Remove output BAM files.
run!: ${BAM} ${BAM}.bai
	rm -rf $^

# Alternative rule for run!
clean: run!

# Generate alignment statistics.
stats:
	@echo "==================== MAPPING STATISTICS ===================="
	@samtools flagstat ${BAM}
	@echo "============================================================"

# Alias for 'stats'.
stat: stats

# Run test suite.
test:
	# Retrieve FASTQ file.
	make -f ${BIOFLOW}/src/fetch/sra.mk SRR=SRR31340505 MODE=SE run
	# Retrieve reference genome.
	make -f ${BIOFLOW}/src/fetch/genbank.mk ACC=ON963982 fasta
	# Map reads to genome with minimap2.
	make -f ${BIOFLOW}/src/minimap2.mk ref=refs/on963982.fa r1=reads/srr31340505.fastq align stats

	@echo "All tests ran successfully!"

	# Clean up files.
	make -f ${BIOFLOW}/src/fetch/sra.mk SRR=SRR31340505 MODE=SE clean
	make -f ${BIOFLOW}/src/fetch/genbank.mk ACC=ON963982 clean
	make -f ${BIOFLOW}/src/minimap2.mk ref=refs/on963982.fa r1=reads/srr31340505.fastq clean

DEPS := minimap2 samtools
# Show installation command.
install::
	micromamba create -n ${ENV}
	${ENV_RUN} micromamba install ${DEPS} --yes

# Targets that are not files.
.PHONY: help run run! align align! stats stat install index test
