---
downloads:
  - file: ../../src/assembly/megahit.mk
    title: megahit.mk
  - file: ../../src/assembly/minia.mk
    title: minia.mk
  - file: ../../src/assembly/spades.mk
    title: spades.mk
  - file: ../../envs/bf-assembly.yml
    title: env.yml
---

(bf-assembly)=
# assembly

## Overview

The `assembly` module contains recipes for performing _de novo_ and reference-based sequence assembly.

:::{note} TODO
Rules for reference-based assembly is in progress.
:::

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

### megahit.mk

Run a _de novo_ assembly pipeline using `megahit`.

**{sc}`Parameters`**

- R1: first set of pair-end reads.
- R2: second set of pair-end reads, if any.
- KMIN: minimum size of k-mer list (default: 21).
- KMAX: maximum size of k-mer list (default: 21).
- KSTEP: step size between k-mers in k-mer list (default: 21).
- OUTDIR: path to directory for storing output (default: megahit).
- MEM: maximum memory to allocate (default: 0.5)
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Assemble pair-end reads into contigs using `megahit`:
```bash
make -f megahit.mk R1=reads_1.fq R2=reads_2.fq run
```

Allocate a fixed amount of memory and CPU cores for the program:
```bash
make -f megahit.mk R1=reads_1.fq R2=reads_2.fq THREADS=8 MEM=0.8 run
```

### minia.mk

Run a _de novo_ assembly pipeline using `minia`. 

By default, the `minia` assembler only accepts a single FASTQ file as input.

**{sc}`Parameters`**

- R1: FASTA file containing reads.
- K: specify k-mer size (default: 31).
- OUTDIR: path to directory for storing output (default: megahit).
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Assemble inter-leaved reads into contigs using `minia`:
```bash
make -f minia.mk R1=reads.fq run
```

Specify k-mer size:
```bash
make -f minia.mk R1=reads.fq K=31 run
```

### spades.mk

Run a _de novo_ assembly pipeline using `spades`.

**{sc}`Parameters`**

- R1: first set of pair-end reads.
- R2: second set of pair-end reads, if any.
- PLATFORM: sequencing platform used in reads (optional).
- K: specify k-mer size (default: auto).
- OUTDIR: path to directory for storing output (default: spades).
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Assemble pair-end reads into contigs using `spades`:
```bash
make -f spades.mk R1=reads_1.fq R2=reads_2.fq run
```

Specify k-mer size:
```bash
make -f spades.mk R1=reads_1.fq R2=reads_2.fq K=31 run
```
