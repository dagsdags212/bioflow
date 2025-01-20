---
downloads:
  - file: ../../src/alignment/Makefile
    title: Makefile
  - file: ../../envs/bf-alignment.yml
    title: env.yml
---

(bf-alignment)=
# alignment.mk

## Overview

The `alignment` module contains rules for performing pairwise and multiple sequence alignment. Its entry point is the `bf-align` command.

Supported alignment tools:

- ClustalO
- MAFFT
- MUSCLE

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bf-align init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-alignment
```
:::

## Commands

### pairwise

Perform alignment on two input sequences.

**{sc}`Parameters`**

- S1: a reference sequence
- S2: a query sequence

**{sc}`Example Usage`**

Align two sequences given by `S1` and `S2`:
```bash
bf-align pairwise S1=ACGTAG S2=TAGTAC
```

### msa

Align sequences using one of the supported aligners.

A multi-sequence FASTA file can serve as input thru the `FA` argument. Alternatively, specify a directory path containing multiple FASTA files for alignment using the `DIR` argument.

The `align` command generates three output files containing the same alignment data but in different formats (fasta, clustal, phylip).

**{sc}`Parameters`**

- FA: path to a multi-sequence FASTA file
- DIR: path to a directory containing multiple FASTA files
- ALIGNER: specify tool for performing alignment (default: mafft)
- GOP: gap opening penalty (default: 1.53)
- GEP: gap extension penalty (default: 0.0)
- PERM: maximum number of iteration for refinement (default: 1)
- OUTNAME: filename of output alignment file (default: aln)
- THREADS: number of cores (default: 8)

**{sc}`Example Usage`**

Align a multi-sequence FASTA file using `mafft`.
```bash
bf-align msa FA=seqs.fa ALIGNER=mafft
```

Customize scoring function for gap openings and extensions.
```bash
bf-align msa FA=seqs.fa GOP=2 GEP=0.5 ALIGNER=MAFFT
```

Specify output filename to generate _aln.fasta_, _aln.phylip_, and _aln.clustal_.
```bash
bf-align msa FA=seqs.fa ALIGNER=mafft OUTNAME=aln
```

Align all FASTA files located within the `seqs/` directory using MUSCLE and set number of threads to 8:
```bash
bf-align msa DIR=seqs/ ALIGNER=muscle THREADS=8
```

### list

Display all currently supported alignment tools.

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**

List available aligners.
```bash
bf-align list
```
