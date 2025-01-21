# Number of CPU cores.
THREADS ?= 4

# Path(s) to FASTQ single-end files.
SE ?=

# Path(s) to FASTQ pair-end #1 files.
PE1 ?=

# Path(s) to FASTQ pair-end #2 files.
PE2 ?=

# K-mer size.
KMER_SIZE ?= 31

# Path to output directory.
OUTDIR ?= output/minia

# Add prefix to minia output.
PREFIX ?= final

# Target FASTA file as output.
CONTIGS = $(OUTDIR)/$(PREFIX).contigs.fa

# Run assembly pipeline with all its dependencies.
all: $(CONTIGS)

# minia input.
minia_input := $(if $(SE),-in $(SE))
minia_input := $(if $(PE1),-in $(PE1) -in $(PE2))

# minia options.
minia_opts += -nb-cores $(THREADS) -verbose 2
minia_opts += -kmer-size $(KMER_SIZE)
minia_opts += -out "$(PREFIX)"
minia_opts += -out-dir $(OUTDIR)

# Run assembly with minia.
$(CONTIGS): $(if $(SE),$(SE),$(PE1) $(PE2))
	@echo "====================================="
	@echo "     Running assembly with MINIA     "
	@echo "====================================="
	@sleep 2
	@minia $(minia_input) $(minia_opts)
	@mv $(PREFIX)* $(OUTDIR)

# Print help message.
help:
	@echo "Perform sequence assembly with minia" 
	@echo
	@echo "Usage:"
	@echo "  make -f run_minia.mk [help] [options] [PE1=<pe1> PE2=<pe2> | SE=<se>]"
	@echo "  <pe1>  a command-separated list of FASTQ pair-end #1 files"
	@echo "  <pe2>  a command-separated list of FASTQ pair-end #2 files"
	@echo "  <se>   a command-separated list of FASTQ single-end files"
	@echo
	@echo "Options:"
	@echo "  K=<int>          k-mer size (default: 31)"
	@echo "  THREADS=<int>    number of threads (default: 4)"
	@echo "  OUTDIR=<path>    directory to store output (default: output/minia)"

# Remove files associated with the minia run.
clean:
	rm -rf $(OUTDIR)
