# Number of CPU cores.
THREADS ?= 4

# Path(s) to FASTQ single-end files.
SE ?=

# Path(s) to FASTQ pair-end #1 files.
PE1 ?=

# Path(s) to FASTQ pair-end #2 files.
PE2 ?=

# Minimun contig length to filter.
MIN_CONTIG_LEN ?= 200

# Path to output directory.
OUTDIR ?= output/megahit

# Target FASTA file as output.
CONTIGS = $(OUTDIR)/final.contigs.fa

# Run assembly pipeline with all its dependencies.
all: $(CONTIGS)

# Megahit input.
megahit_input := $(if $(SE),-r $(SE))
megahit_input += $(if $(PE1),-1 $(PE1) -2 $(PE2))

# Megahit options.
megahit_opts := -t $(THREADS)
megahit_opts += --min-contig-len $(MIN_CONTIG_LEN)
megahit_opts += -o $(OUTDIR)

# Run assembly with megahit.
$(CONTIGS): $(if $(SE),$(SE),$(PE1) $(PE2))
	@rm -rf $(dir $(OUTDIR))
	@mkdir -p $(dir $(OUTDIR))
	@echo "====================================="
	@echo "    Running assembly with MEGAHIT    "
	@echo "====================================="
	@sleep 2
	@megahit $(megahit_input) $(megahit_opts)

# Print help message.
help:
	@echo "Perform sequence assembly with Megahit" 
	@echo
	@echo "Usage:"
	@echo "  make -f run_megahit.mk [help] [options] [PE1=<pe1> PE2=<pe2> | SE=<se>]"
	@echo "  <pe1>  a command-separated list of FASTQ pair-end #1 files"
	@echo "  <pe2>  a command-separated list of FASTQ pair-end #2 files"
	@echo "  <se>   a command-separated list of FASTQ single-end files"
	@echo
	@echo "Options:"
	@echo "  MIN_CONTIG_LEN=<int>    minimum contig length (default: 200)"
	@echo "  THREADS=<int>           number of threads (default: 4)"
	@echo "  OUTDIR=<path>           directory to store output (default: output/megahit)"

# Remove files associated with the Megahit run.
clean:
	rm -rf $(OUTDIR)
