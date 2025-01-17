# bioflow

Build your pipeline with a set of modular components for running common tasks in bioinformatics.

**bioflow** is a collection of [Makefiles](https://www.gnu.org/software/make/manual/make.html) for orchestrating workflows in bioinformatics. Each Makefile reprents an module that packages together software tools for conducting a specific task in genomics.

## Modules

| Makefile | Functionality |
| ---- | ------------ |
| `alignment.mk` | pairwise alignment, multiple sequence alignment |
| `annotation.mk` | gene prediction, assembly annotation |
| `assembly.mk` | sequence assembly, assembly evaluation |
| `fetch.mk` | download biological data |
| `mapping.mk` | read mapping, mapping assessment |
| `phylo.mk` | maximum likelihood inference, bootstrapping phylogenies |
| `qc.mk` | produce FASTQ summary reports, adapter trimming, quality filtering |
| `variant_calling.mk` | detecting and filtering variants |
