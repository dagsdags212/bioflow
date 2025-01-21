# Path to reference file (FASTA).
REF ?= 

# Path to single-end FASTQ file.
SE ?=

# Path to pair-end #1 FASTQ file.
PE1 ?=

# Path to pair-end #2 FASTQ file.
PE2 ?=

# Path to directory for storing output.
OUTDIR ?= output/bwa

# Filename for output alignment.
PREFIX ?= bwa_aln

# Number of CPU cores.
THREADS ?= 2

# Path to final alignment file.
ALN = $(OUTDIR)/$(PREFIX).sorted.bam

# Expected suffixes for bwa-indexed reference.
BWA_SFX = amb ann bwt pac sa
IDX = $(addprefix $(REF).,$(BWA_SFX))

all: $(ALN)

# Index reference file.
$(IDX): $(REF)
	bwa index $(REF)

# Map reads to reference using bwa.
$(ALN): $(IDX)
	@echo "==========================="
	@echo "  Aligning reads with bwa  "
	@echo "==========================="
	@sleep 2
	@mkdir -p $(OUTDIR)
ifdef SE
	@bwa mem $(REF) $(SE) -t $(THREADS) | samtools sort - | samtools view -Sbh > $(OUTDIR)/$(PREFIX).sorted.bam
else ifdef PE1
	@bwa mem $(REF) $(PE1) $(PE2) -t $(THREADS) | samtools sort - | samtools view -Sbh > $(OUTDIR)/$(PREFIX).sorted.bam
endif
