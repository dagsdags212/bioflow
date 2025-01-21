# Path to reference file (FASTA).
REF ?= 
# Path to single-end FASTQ file.
SE ?=

# Path to pair-end #1 FASTQ file.
PE1 ?=

# Path to pair-end #2 FASTQ file.
PE2 ?=

# Path to directory for storing output.
OUTDIR ?= output/bowtie

# Filename for output alignment.
PREFIX ?= bowtie_aln

# Number of CPU cores.
THREADS ?= 2

# Path to final alignment file.
ALN = $(OUTDIR)/$(PREFIX).sorted.bam

# Expected suffixes for bwa-indexed reference.
BOWTIE_SFX = 1.bt2 2.bt2 3.bt2 4.bt2 rev.1.bt2 rev.2.bt2
IDX = $(addprefix $(REF).,$(BOWTIE_SFX))

all: $(ALN)

# Index reference file.
$(IDX): $(REF)
	bowtie2-build $(REF) $(basename $(REF))

# Map reads to reference using bwa.
$(ALN): $(IDX)
	@echo "=============================="
	@echo "  Aligning reads with bowtie  "
	@echo "=============================="
	@sleep 2
	@mkdir -p $(OUTDIR)
ifdef SE
	@bowtie2 -x $(basename $(REF)) -U $(SE) -p $(THREADS) | samtools sort - | samtools view -Sbh > $(OUTDIR)/$(PREFIX).sorted.bam
else ifdef PE1
	@bowtie2 -x $(basename $(REF)) -1 $(PE1) -2 $(PE2) -p $(THREADS) | samtools sort - | samtools view -Sbh > $(OUTDIR)/$(PREFIX).sorted.bam
else
	@echo "Error: path of reads not specified"
endif
