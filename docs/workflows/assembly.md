---
downloads:
  - file: ../../src/assembly/Makefile
    title: Makefile
  - file: ../../envs/bf-assembly.yml
    title: env.yml
---

(bf-assembly)=
# assembly.mk

## Overview

The `assembly` workflow contains rules for performing _de novo_ and reference-based sequence assembly. Its entry point is the `bf-assemble` command.

:::{note} TODO
Rules for reference-based assembly is in progress.
:::

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bf-assemble init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-assembly
```
:::

## Rules

### assemble

Run a _de novo_ assembly pipeline using `spades` and `megahit`. Each tool is designated a directory under `assembly` to store all output.

**{sc}`Parameters`**

- ASSEMBLER: specified tool for running the assembly (default: megahit).
- PE1: path to reads for first pair-end file.
- PE2: path to reads for pair-end file.
- SE: path to single-end reads; PE1 and PE2 are ignored if SE is specified.
- MIN_CONTIG_LEN: minimum contig length for filtering contigs (default: 200).
- KMER_SIZE: k-mer length to use (default: 31).
- PLATFORM: specify the sequencing platform to be used by spades.
- PREFIX: prepend a prefix to the output files of minia.
- OUTDIR: path to directory for storing output.
- THREADS: number of cores (default: 8)

**{sc}`Example Usage`**

Assemble pair-end reads into contigs using `megahit`:
```bash
bf-assemble assemble PE1=read1.fq PE2=read2.fq ASSEMBLER=megahit
```

Utilize 8 threads to assemble single-end reads into contigs using `spades`. Save the output to the `output/spades` directory:
```bash
bf-assemble assemble SE=reads.fq ASSEMBLER=spades THREADS=8 OUTDIR=output/spades
```

Assemble single-end reads into contigs using `minia` and specify a prefix:
```bash
bf-assemble assemble SE=reads.fq PREFIX=final ASSEMBLER=minia
```
