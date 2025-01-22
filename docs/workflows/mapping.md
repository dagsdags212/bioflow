---
downloads:
  - file: ../../src/mapping/Makefile
    title: Makefile
  - file: ../../envs/mapping.yml
    title: env.yml
---

(bf-mapping)=
# mapping.mk

## Overview

The `mapping` module contains rules for performing read alignment against a reference genome. Its entry point is the `bf-map` command.

Supported aligners include:

- [bwa](https://github.com/lh3/bwa)
- [bowtie2](https://github.com/BenLangmead/bowtie2)

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bf-map init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-mapping
```
:::

## Rules

### map

Map reads against a reference using `bwa`. 

The reference file is first indexed prior to mapping. By default, the output alignment is sorted and in SAM format.

**{sc}`Parameters`**

- MAPPER: specified tool for performing read mapping (default: bwa).
- REF: path to reference file (FASTA)
- PE1: path to reads for first pair-end file.
- PE2: path to reads for pair-end file.
- SE: path to reads for single-end file.
- OUTDIR: path to directory for storing output.
- PREFIX: filename of resulting alignment file.
- THREADS: number of cores (default: 8)

**{sc}`Example Usage`**
Map pair-end reads against a reference using `bwa`:
```bash
bf-map map PE1=reads1.fq PE2=reads2.fq REF=ref.fa MAPPER=bwa
```

Map pair-end reads against a reference using `bowtie`. Specify the number of worker threads:
```bash
bf-map map PE1=reads1.fq PE2=reads2.fq REF=ref.fa MAPPER=bowtie \
    THREADS=8 PREFIX=bt OUTDIR=output/bowtie
```

Map single-end reads against a reference using `bowtie`:
```bash
bf-map map SE=reads.fq REF=ref.fa MAPPER=bowtie
```

Specify output name and directory for the alignment file:
```bash
bf-map map SE=reads.fq REF=ref.fq MAPPER=bowtie \
    PREFIX=bt OUTDIR=output/bam THREADS=8
```

### stats

Get an overview of the mapped reads metrics.

All sorted BAM files within the current working directory are automatically discovered. Mapping metrics for each file are printer to `stdout`.

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**
Print mapping metrics for all generated alignment files.
```bash
bf-map stats
```

### evaluate

Assess the quality of the alignment file using `qualimap`.

Similar to `stats`, this command checks the current working directory for alignment files. For each BAM file, a summary report is generated in HTML format.

```{note}
BAM files must be sorted prior to evaluation.
```

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**
Run `qualimap` on all sorted BAM files.
```bash
bf-map evaluate
```
