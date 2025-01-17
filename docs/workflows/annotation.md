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

The `annotation.mk` workflow contains rules for conducting gene prediction and sequence annotation. Its entry point is the `bf-annotate` command.

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bf-annotate init ENV_MANAGER=micromamba
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
bf-annotate annotate \
    FA=output/megahit/final.contigs.fa MODE=genome
```

Specify domain where the sequence file belongs to.
```bash
bf-annotate annotate \
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
bf-annotate datasets
```

Delete lineage data.
```bash
bf-annotate clean
```

### predict

Conduct HMM-based gene prediction using `prodigal` and store output in GFF format by default.

**{sc}`Parameters`**

- FA: path to target sequence file
- PRODIGAL_OUTFMT: format of gene prediction output

**{sc}`Example Usage`**

Run gene prediction on a draft assembly.
```bash
bf-annotate predict FA=output/megahit/final.contigs.fa
```

Save output in GenBank Flat File Format (gbk).
```bash
bf-annotate predict FA=output/megahit/final.contigs.fa PRODIGAL_OUTFMT=gbk
```

### view

Visualize annotated genome using Artemis.

Run `annotate` prior to invoking this command to generate the `prokka` annotation.

**{sc}`Parameters`**

- FA: path to target sequence file

**{sc}`Example Usage`**

Visualize annotation by `prokka` alongside the reference genome.
```bash
bf-annotate view FA=output/megahit/final.contigs.fa
```

