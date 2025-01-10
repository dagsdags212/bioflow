---
downloads:
  - file: ../../src/qc.mk
    title: Makefile
  - file: ../../envs/qc.yml
    title: env.yml
---

(bwf-qc)=
# qc.mk

The `qc.mk` workflow contains rules for performing quality control on sequencing reads.

## fastp

Perform automatic adapter trimming and quality score filtering.

Generates a set of trimmed/filtered reads and summary reports as output. These are stored in the `fastp/reads` and `fastp/reports` directories, respectively.

**{sc}`Parameters`**

- READ_DIR: path to directory storing reads
- MINLEN: minimum read length
- MAXLEN: maximum read length
- MINQUAL: minimum quality score threshold
- PE: treat reads as pair-end (default: true)

**{sc}`Example Usage`**
```bash
make -f src/qc.mk fastp READ_DIR=reads/
```

## fastqc

Generate a FASTQC report on selected reads.

FASTQ files are automatically discovered from the given directory and all detected files are ran through the program. Reports are stored in a `fastqc` directory.

**{sc}`Parameters`**

- READ_DIR: path to directory storing reads

**{sc}`Example Usage`**
Generate reports from trimmed reads outputted by `fastp`.
```bash
make -f src/qc.mk fastqc READ_DIR=fastp/reads
```

## multiqc

Aggregate FASTQC output into a single interactive report. 

Expects a `fastqc` directory containing the summary report from FASTQC. Output is stored in a `multiqc` directory.

**{sc}`Example Usage`**
```bash
make -f src/qc.mk multiqc
```
