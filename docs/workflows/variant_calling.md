---
downloads:
  - file: ../../src/variant_calling.mk
    title: Makefile
  - file: ../../envs/vc.yml
    title: env.yml
---

(bf-vc)=
# variant_calling.mk

## Overview

The `variant_calling.mk` workflow contains rules for identifying variants from mapped reads based on a reference genome. Its entry point is the `bf-vc` command.

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bf-vc init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-vc
```
:::

## Rules

### call

Perform variant calling with `freebayes`.

This commands expects a reference file in FASTA format and an alignment file in BAM format. Both reference and alignment files are indexed prior to calling variants.

**{sc}`Parameters`**

- REF: path to reference genome (FASTA)
- BAM: path to alignment file (BAM)
- PLOIDY: ploidy number of reference genome (default: 1)
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Call variants from a BAM file.
```bash
bf-vc call \
    REF=ref/ref.fna BAM=output/bwa/aln.bam
```

Indicate ploidy number of source organism.
```bash
bf-vc call \
    REF=ref/ref.fna BAM=output/bwa/aln.bam PLOIDY=1
```

### stats

Generate metric file for called variants.

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**

Produce a summary file containing metrics on called variants.
```bash
bf-vc stats \
    REF=ref/ref.fna BAM=output/bwa/aln.bam PLOIDY=1
```

### filter

Filter called variants based on a score threshold.

**{sc}`Parameters`**

- MINQUAL: variant score threshold for filtering (default: 10)

**{sc}`Example Usage`**

Filter variants with a score less than 15.
```bash
bf-vc filter \
    REF=ref/ref.fna BAM=output/bwa/aln.bam \
    MINQUAL=15
```
