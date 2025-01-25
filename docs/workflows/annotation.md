---
downloads:
  - file: ../../src/annotation/busco.mk
    title: busco.mk
  - file: ../../src/annotation/prodigal.mk
    title: prodigal.mk
  - file: ../../src/annotation/prokka.mk
    title: prokka.mk
  - file: ../../envs/bf-annotation.yml
    title: env.yml
---

(bf-annotation)=
# annotation

## Overview

The `annotation` module contains rules for performing _de novo_ gene prediction and sequence annotation.

:::{hint} Environment Setup
:class: dropdown

Display the command to install all dependencies:
```bash
make -f <path/to/makefile> install
```

To download all dependencies, run the following command:
```bash
eval $(make -f <path/to/makefile> install)
```

Make sure to replace `<path/to/makefile>` with a path pointing to your recipe.
:::

## Recipes

### busco.mk

Run the annotation pipeline provided by `busco`.

**{sc}`Parameters`**

- IN: an input sequence file or directory.
- MODE: specify BUSCO analysis mode (default: genome)
- CUTOFF: an E-value threshold for predicted genes (default: 1e-06)
- OUTDIR: an output directory for storing BUSCO output (default: busco)
- THREADS: number of worker threads (default: 2)

**{sc}`Example Usage`**

Run annotation pipeline on a draft genome assembly:
```bash
make -f busco.mk IN=final.contigs.fa MODE=genome run
```

Specfy an output directory:
```bash
make -f busco.mk IN=final.contigs.fa MODE=genome OUTDIR=annot run
```

List available BUSCO databases:
```bash
make -f busco.mk list
```

### prodigal.mk

Perform _ab initio_ gene prediction with `prodigal`.

**{sc}`Parameters`**

- IN: an input sequence file or directory.
- CUTOFF: an E-value threshold for predicted genes (default: 1e-06).
- FORMAT: specify an output format (default: gbk).
- OUTDIR: an output directory for storing BUSCO output (default: prodigal).
- THREADS: number of worker threads (default: 2).

**{sc}`Example Usage`**

Predict genes from an assembly and store in _gff_ format:
```bash
make -f prodigal.mk IN=final.contigs.fa FORMAT=gff run
```

Specify an output directory:
```bash
make -f prodigal.mk IN=final.contigs.fa FORMAT=gff OUTDIR=genes run
```

### prokka.mk

Annotate a prokaryotic genome with `prokka`.

**{sc}`Parameters`**

- IN: an input sequence file or directory.
- KINGDOM: Prokka annotation mode (default: bacteria).
- GENUS: genus of source organism (optional).
- SPECIES: species of source organism (optional).
- STRAIN: strain of source organism (optional).
- PLASMID: plasmid of source organism (optional).
- CUTOFF: an E-value threshold for predicted genes (default: 1e-06).
- OUTDIR: an output directory for storing BUSCO output (default: busco).
- THREADS: number of worker threads (default: 2).

**{sc}`Example Usage`**

Run annotation pipeline on a draft genome assembly:
```bash
make -f prokka.mk IN=final.contigs.fa MODE=genome run
```

List available Prokka databases:
```bash
make -f prokka.mk list
```
