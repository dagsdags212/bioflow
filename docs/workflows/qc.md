---
downloads:
  - file: ../../src/qc/fastp.mk
    title: fastp.mk
  - file: ../../src/qc/fastqc.mk
    title: fastqc.mk
  - file: ../../src/qc/multiqc.mk
    title: multiqc.mk
  - file: ../../envs/bf-qc.yml
    title: env.yml
---

(bf-qc)=
# qc

## Overview

The `qc` module contains recipes for performing quality control on sequencing reads.

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

### fastp

Perform adapter trimming and quality score filtering with `fastp`.

By default, the processed set of reads are saved in the same directory as the raw reads. As of now, there is no option to specify the output directory.

**{sc}`Parameters`**

- `P1`: path to the first FASTQ file in pair-end reads.
- `P2`: path to the FASTQ file in pair-end reads; if not specified, `fastp` is run in single-end mode.

**{sc}`Example Usage`**

Trim reads in pair-end mode.
```bash
SRR=SRR1554324
make -f fastp.mk P1=${SRR}_1.fastq P2=${SRR}_2.fastq run
```

Trim reads in single-end mode.
```bash
SRR=SRR1554324
make -f fastp.mk P1=${SRR}.fastq run
```

### fastqc

Generate a summary report on read quality.

FASTQ files are automatically discovered from the path specified by the `TARGET` parameter. If not given, the target defaults to the `READ_DIR` value. Reports are stored in the `output/fastqc/report#` directory where `#` is an auto-incrementing integer.

**{sc}`Parameters`**

- `DIR`: path to directory containing the reads (default: reads).
- `OUT`: output directory for reports (default: fastqc).

**{sc}`Example Usage`**

Generate a FASTQC report on all reads stored in the `raw` directory:
```bash
make -f fastqc.mk DIR=raw run
```

Save reports in the `reports` directory:
```bash
make -f fastqc.mk DIR=raw OUT=reports run
```

### multiqc

Aggregate FASTQC reports into a single interactive interface. 

**{sc}`Parameters`**

- `DIR`: path to directory containing the FASTQC reports (default: fastqc).
- `OUT`: output directory for consolidate report (default: multiqc).

**{sc}`Example Usage`**

Consolidate all FASTQC reports within the `fastqc` directory:
```bash
make -f multiqc DIR=fastqc run
```
