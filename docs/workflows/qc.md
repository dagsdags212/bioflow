---
downloads:
  - file: ../../src/qc/fastp.mk
    title: fastp.mk
  - file: ../../src/qc/fastqc.mk
    title: fastqc.mk
  - file: ../../src/qc/multiqc.mk
    title: multiqc.mk
---

(bf-qc)=
# qc

## Overview

The `qc` module contains recipes for performing quality control on sequencing reads.

:::{hint} Environment Setup
:class: dropdown

Initialize the micromamba environment and install the dependencies:
```{code-cell} bash
bf-<recipe> install
```

where `<recipe>` is the name of the specific recipe.

:::

## Recipes

### bf-fastp

Perform adapter trimming and quality score filtering with `fastp`.

By default, the processed set of reads are saved in the same directory as the raw reads. As of now, there is no option to specify the output directory.

**{sc}`Parameters`**

- `P1`: path to the first FASTQ file in pair-end reads.
- `P2`: path to the FASTQ file in pair-end reads; if not specified, `fastp` is run in single-end mode.
- `ADAPTER`: an adapter file or sequence for trimming
- `MINLEN`: minimum read length to keep (default: 50)
- `MINQUAL`: minimum read quality to keep (default: 30)
- `OUT`: path to directory for storing fastp output
- `THREADS`: number of worker threads to use (default: 4).

**{sc}`Example Usage`**

Trim reads in pair-end mode.
```{code-cell} bash
SRR=SRR1554324
bf-fastp P1=${SRR}_1.fastq P2=${SRR}_2.fastq run
```

Trim reads in single-end mode.
```{code-cell} bash
SRR=SRR1554324
bf-fastp P1=${SRR}.fastq run
```

### bf-fastqc

Generate a summary report on read quality.

FASTQ files are automatically discovered from the path specified by the `DIR` parameter. If not given, FASTQ files are searched in the _reads_ directory. Reports are stored in the `fastqc` directory by default.

**{sc}`Parameters`**

- `DIR`: path to directory containing the reads (default: reads).
- `CONTAMINANT_FILE`: a file containing contaminant sequences for filtering
- `OUT`: output directory for reports (default: fastqc).
- `THREADS`: number of worker threads to use (default: 4).

**{sc}`Example Usage`**

Generate a FASTQC report on all reads stored in the `raw` directory:
```{code-cell} bash
bf-fastqc DIR=raw run
```

Save reports in the `reports` directory:
```{code-cell} bash
bf-fastqc DIR=raw OUT=reports run
```

### bf-multiqc

Aggregate FASTQC output into a single interactive report.

**{sc}`Parameters`**

- `DIR`: path to directory containing the FASTQC reports (default: fastqc).
- `OUT`: output directory for consolidate report (default: multiqc).

**{sc}`Example Usage`**

Consolidate all FASTQC reports within the `fastqc` directory:
```{code-cell} bash
bf-multiqc DIR=fastqc run
```
