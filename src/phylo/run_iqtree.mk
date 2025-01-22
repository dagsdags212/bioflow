# Sequence alignment file.
MSA ?=

# Output directory.
OUTDIR ?= output/iqtree

# Name of output file.
PREFIX ?= result

# Infer output target from OUTDIR and PREFIX.
TREE = $(OUTDIR)/$(PREFIX)

# Tree construction method.
METHOD ?= ml

# Starting seed.
SEED ?= 12345

# Nucleotide substitution model.
MODEL ?= GTR+G

# Number of bootstraps.
BOOTSTRAPS ?= 1001

# Number of CPU cores.
THREADS ?= 4

iqtree_sfx := iqtree fa alninfo bionj log mldist treefile ckp.gz
iqtree_output := $(addprefix $(PREFIX).,$(iqtree_sfx))

all: $(firstword $(iqtree_output))

$(iqtree_output): $(MSA)
	@mkdir -p $(OUTDIR)
	
# Run maximum-likelihood algorithm.
ifeq ($(METHOD),ml)
	@echo "==================================="
	@echo "  Generating ML tree using IQTREE  "
	@echo "==================================="
	@sleep 1
	@iqtree -s $(MSA) --alisim $(OUTDIR)/$(PREFIX) --prefix $(PREFIX) \
		-T $(THREADS) -m $(MODEL) --seed $(SEED) -af fasta -alninfo
	@mv $(PREFIX)* $(OUTDIR)

# Run bootstrap algorithm.
else ifeq ($(METHOD),bootstrap)
	@echo "==============================================="
	@echo "  Generating tree from bootstrap using IQTREE  "
	@echo "==============================================="
	@sleep 1
	@iqtree -s $(MSA) --alisim $(OUTDIR)/$(PREFIX) -B $(BOOTSTRAPS) --prefix $(PREFIX) \
		-T $(THREADS) -m $(MODEL) --seed $(SEED) -af fasta -alninfo
	@mv $(PREFIX)* $(OUTDIR)
else
	@echo "Error: invalid tree construction method"
endif

version:
	@iqtree -v | head -n 3

clean:
	rm -rf $(OUTDIR)

.PHONY: clean version
