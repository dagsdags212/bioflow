#
# Map reads against a reference using minimap2.
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
ENV = bf-mapping

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# Number of worker threads.
THREADS ?= 4

# Minimap2 options.
FLAGS := -a -t $(THREADS)

# FASTQ file containing reads for alignment.
FQ ?=

# Genbank accession of the reference file.
ACC ?= ON963982

# The reference genome.
REF ?= refs/$(ACC).fa

# The alignment file.
BAM ?= bam/aln.bam

# The name of the stats file.
STATS = $(basename $(BAM)).stats

# Print the help message.
help::
	@echo "#"
	@echo "# minimap2.mk: align long reads using minimap2 "
	@echo "#"
	@echo "# REF=${REF}"
	@echo "#"
	@echo "# FQ=${FQ}"
	@echo "#"
	@echo "# BAM=${BAM}"
	@echo "#"
	@echo "# make index|run|test|clean"
	@echo "#"

# Read 1 must exist.
$(FQ):
	@echo "# Error: FASTQ file not found (FQ=$(FQ))"
	@exit -1

# Reference file must exist.
$(REF):
	@if [ ! -f $(REF) ]; then
		echo "# Error: Reference file not found (REF=$(REF))"
		exit -1
	fi

# Compose minimap2 command.
MINIMAP2_CMD = minimap2 $(FLAGS) $(REF) $(FQ)

# Generate a sorted alignment file.
$(BAM): $(FQ)
	# Output directory for BAM files.
	mkdir -p $(dir $@)

	# Perform read mapping.
	$(ENV_RUN) $(MINIMAP2_CMD) | samtools sort -@ $(THREADS) --output-fmt=BAM > $@

# Create the BAM index.
$(BAM).bai: $(BAM)
	samtools index $<

# Invoke the read mapping rule.
align: $(BAM) $(BAM).bai
	@ls -lh $^

# Alternative rule for align.
run: align

# Remove output BAM files.
run!: $(BAM) $(BAM).bai
	rm -rf $^

# Alternative rule for run!
clean: run!

# Generate alignment statistics.
$(STATS): $(BAM).bai
	samtools flagstat $(BAM) > $(STATS)

# Trigger stats generation.
stats: $(STATS)
	@echo "# $(STATS)"
	@cat $(STATS)

# Run test suite.
test:
	# Retrieve FASTQ file.
	make -f $(shell dirname $(ROOT_PATH))/fetch/sra.mk SRR=SRR31340505 run
	# Retrieve reference genome.
	make -f $(shell dirname $(ROOT_PATH))/fetch/genbank.mk ACC=ON963982 fasta
	# Map reads to genome with minimap2.
	make -f $(root_path)/minimap2.mk ref=refs/on963982.fa fq=reads/srr31340505_1.fastq align
	# Generate mapping statistics.
	make -f $(ROOT_PATH)/minimap2.mk REF=refs/ON963982.fa FQ=reads/SRR31340505_1.fastq stats

	@echo "All tests ran successfully!"
	# Clean up files.
	make -f $(shell dirname $(ROOT_PATH))/fetch/sra.mk SRR=SRR31340505 clean
	make -f $(shell dirname $(ROOT_PATH))/fetch/genbank.mk ACC=ON963982 clean
	make -f $(ROOT_PATH)/minimap2.mk REF=refs/ON963982.fa FQ=reads/SRR31340505_1.fastq clean

# Show installation command.
install::
	@echo "micromamba install minimap2 samtools"

# Targets that are not files.
.PHONY: run run! install help index
