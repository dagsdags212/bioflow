---
downloads:
  - file: ../../src/mapping.mk
    title: Makefile
  - file: ../../envs/mapping.yml
    title: env.yml
---

(bf-mapping)=
# mapping.mk

## Overview

The `mapping.mk` workflow contains rules for performing read alignment against a reference genome. Its entry point is the `bf-map` command.

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

- READ_DIR: path to directory storing reads
- REF: path to reference file (FASTA)
- OUTFMT: alignment format to output (default: sam)
- MAPPER: specified tool for read mapping (default: bwa)
- THREADS: number of cores (default: 8)
- PE: treat reads as pair-end (default: true)

**{sc}`Example Usage`**
Map pair-end reads against a reference using `bwa`.
```bash
bf-map map \
    READ_DIR=reads/ REF=ref.fa PE=true \
    MAPPER=bwa THREADS=8
```

Map pair-end reads against a reference using `bowtie`.
```bash
bf-map map \
    READ_DIR=reads/ REF=ref.fa PE=true \
    MAPPER=bowtie THREADS=8
```

Map single-end reads against a reference and save as a BAM file.
```bash
bf-map map \
    READ_DIR=reads/ REF=ref.fa PE=false \
    OUTFMT=bam MAPPER=bwa THREADS=8
```

Specify output name for alignment file.
```bash
bf-map map \
    READ_DIR=reads/ REF=ref.fa PE=false \
    OUTFMT=bam OUTPUT=myaln \
    MAPPER=bwa THREADS=8
```

### stats

Get an overview of the mapped reads metrics.

All alignment files stored in the `output` directory are discovered and their mapping metrics are printed to stdout.

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**
Print mapping metrics for all generated alignment files.
```bash
bf-map stats
```

### evaluate

Assess the quality of the alignment file using `qualimap`.

Similar to `stats`, this command checks the `output` directory for alignment files. For each BAM file, a summary report is generated in HTML format.

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
