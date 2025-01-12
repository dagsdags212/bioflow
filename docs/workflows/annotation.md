---
downloads:
  - file: ../../src/annotation.mk
    title: Makefile
  - file: ../../envs/annotation.yml
    title: env.yml
---

(bf-annotation)=
# annotation.mk

## Overview

The `annotation.mk` workflow contains rules for conducting gene prediction and sequence annotation.

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
make -f src/annotation.mk init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-annotation
```
:::

## Rules

### annotate

Run a `busco` annotation pipeline on a target sequence file.

All output is stored in the directory at `output/busco`.

**{sc}`Parameters`**

- FA: path to target sequence file
- MODE: specify BUSCO analysis mode (default: genome)
- THREADS: number of cores (default: 8)

**{sc}`Example Usage`**

Run annotation pipeline on a draft assembly.
```bash
make -f src/annotation.mk annotate \
    FA=output/megahit/final.contigs.fa MODE=genome
```

Specify domain where the sequence file belongs to.
```bash
make -f src/annotation.mk annotate \
    FA=output/megahit/final.contigs.fa MODE=genome DOMAIN=prokaryote
```

### datasets

Display all available BUSCO datasets.

By default, it also downloads lineage data within the `busco_downloads` directory.

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**

List all datasets.
```bash
make -f src/annotation.mk datasets
```

Delete lineage data.
```bash
make -f src/annotation.mk clean
```

### predict

Conduct HMM-based gene prediction using `prodigal` and store output in GFF format by default.

**{sc}`Parameters`**

- FA: path to target sequence file
- PRODIGAL_OUTFMT: format of gene prediction output

**{sc}`Example Usage`**

Run gene prediction on a draft assembly.
```bash
make -f src/annotation.mk predict FA=output/megahit/final.contigs.fa
```

Save output in GenBank Flat File Format (gbk).
```bash
make -f src/annotation.mk predict FA=output/megahit/final.contigs.fa PRODIGAL_OUTFMT=gbk
```
