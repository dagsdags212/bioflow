# Number of CPU cores.
THREADS ?= 4

# Path(s) to FASTQ single-end files.
SE ?=

# Path(s) to FASTQ pair-end #1 files.
PE1 ?=

# Path(s) to FASTQ pair-end #2 files.
PE2 ?=

# List of k-mer sizes (must be odd and less than 128).
K ?= auto

# Sequencing platform [sanger|pacbio|nanopore].
PLATFORM ?=

# Path to output directory.
OUTDIR ?= output/spades

# Target FASTA file as output.
CONTIGS = $(OUTDIR)/scaffolds.fasta

# Run assembly pipeline with all its dependencies.
all: $(CONTIGS)

# spades input.
spades_input := $(if $(SE),-s $(SE))
spades_input += $(if $(PE1),-1 $(PE1) -2 $(PE2))

# spades options.
spades_opts := -t $(THREADS)
spades_opts += $(if $(filter sanger,$(PLATFORM)),--sanger)
spades_opts += $(if $(filter pacbio,$(PLATFORM)),--pacbio)
spades_opts += $(if $(filter nanopore,$(PLATFORM)),--nanopore)
spades_opts += -o $(OUTDIR)

# Run assembly with spades.
$(CONTIGS): $(if $(SE),$(SE),$(PE1) $(PE2))
	@echo "===================================="
	@echo "    Running assembly with SPADES    "
	@echo "===================================="
	@sleep 2
	@spades.py $(spades_input) $(spades_opts)

# Print help message.
help:
	@echo "Perform sequence assembly with spades" 
	@echo
	@echo "Usage:"
	@echo "  make -f run_spades.mk [help] [options] [PE1=<pe1> PE2=<pe2> | SE=<se>]"
	@echo "  <pe1>  a command-separated list of FASTQ pair-end #1 files"
	@echo "  <pe2>  a command-separated list of FASTQ pair-end #2 files"
	@echo "  <se>   a command-separated list of FASTQ single-end files"
	@echo
	@echo "Options:"
	@echo "  K=<int>          list of k-mer sizes (must be odd and less than 128)"
	@echo "  PLATFORM=<str>   specify source sequencing platform [sanger|pacbio|nanopore]"
	@echo "  THREADS=<int>    number of threads (default: 4)"
	@echo "  OUTDIR=<path>    directory to store output (default: output/spades)"

# Remove files associated with the spades run.
clean:
	rm -rf $(OUTDIR)
