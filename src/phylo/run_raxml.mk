# Sequence alignment file.
MSA ?=

# Output directory.
OUTDIR ?= output/raxml

# Name of output file.
PREFIX ?= result

# Infer output target from OUTDIR and PREFIX.
TREE = $(OUTDIR)/$(PREFIX)

# Tree construction method.
METHOD ?= ml

# Number of ML searches.
N ?= 10

# Starting seed.
SEED ?= 12345

# Nucleotide substitution model.
MODEL ?= GTR+G

# Number of CPU cores.
THREADS ?= 4

raxml_sfx := bestModel bestTree bestTreeCollapsed log mlTrees rba reduced.phy startTree
raxml_output := $(addprefix $(MSA).,$(raxml_sfx))

all: $(firstword $(raxml_output))

$(raxml_output): $(MSA)
	@mkdir -p $(OUTDIR)
ifeq ($(METHOD),ml)
	@echo "====================================="
	@echo "  Generating ML tree using raxml-ng  "
	@echo "====================================="
	@sleep 1
	@raxml-ng --msa $(MSA) --threads $(THREADS) --model $(MODEL) \
		--tree pars{$(N)} --seed $(SEED) $(TREE)
	@mv $(MSA).raxml.* $(OUTDIR)
else ifeq ($(METHOD),bootstrap)
	@echo "================================================="
	@echo "  Generating tree from bootstrap using raxml-ng  "
	@echo "================================================="
	@sleep 1
	@raxml-ng --msa $(MSA) --threads $(THREADS) --model $(MODEL) \
		--bootstrap --bs-trees $(N) --seed $(SEED) $(TREE)
	@mv $(MSA).raxml.* $(OUTDIR)
else
	@echo "Error: invalid tree construction method"
endif

version:
	@raxml-ng -v

clean:
	rm -rf $(OUTDIR)

.PHONY: clean version
