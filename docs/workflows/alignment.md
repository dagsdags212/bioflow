---
downloads:
  - file: ../../src/alignment.mk
    title: Makefile
  - file: ../../envs/alignment.yml
    title: env.yml
---

(bwf-alignment)=
# alignment.mk

## Overview

The `alignment.mk` workflow contains rules for performing pairwise and multiple sequence alignment.

Supported alignment tools:

- ClustalO
- MAFFT
- MUSCLE

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
make -f src/alignment.mk init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bwf-alignment
```
:::

## Rules

### align

Align sequences using one of the supported aligners.

A multi-sequence FASTA file can serve as input thru the `FA` argument. Alternatively, specify a directory path containing multiple FASTA files for alignment using the `DIR` argument.

**{sc}`Parameters`**

- FA: path to a multi-sequence FASTA file
- DIR: path to a directory containing multiple FASTA files
- ALIGNER: specify tool for performing alignment (default: mafft)
- GOP: gap opening penalty (default: 1.53)
- GEP: gap extension penalty (default: 0.0)
- ITER: maximum number of iteration for refinement (default: 1)
- OUTFMT: format of output alignment file (default: fasta)
- OUTNAME: filename of output alignment file (default: aln)
- THREADS: number of cores (default: 8)

**{sc}`Example Usage`**

Align a multi-sequence FASTA file using `mafft`.
```bash
make -f src/alignment.mk align \
    FA=seqs.fa ALIGNER=mafft
```

Customize scoring function for gap openings and extensions.
```bash
make -f src/alignment.mk align \
    FA=seqs.fa ALIGNER=mafft \
    GOP=2 GEP=0.5
```

Specify output filename and format.
```bash
make -f src/alignment.mk align \
    FA=seqs.fa ALIGNER=mafft \
    OUTFMT=clustal OUTNAME=my_alignment
```

### list

Display all supported alignment tools.

**{sc}`Parameters`**

This command does not accept parameters.

**{sc}`Example Usage`**

List available aligners.
```bash
make -f src/alignment.mk list
```

### view

Visualize the resulting alignment using `jalview`.

The output is stored in the same directory of the alignment file in HTML format.

**{sc}`Parameters`**

- OUTFMT: format of output alignment file (default: fasta)
- OUTNAME: filename of output alignment file (default: aln)

**{sc}`Example Usage`**

Specify `OUTNAME` and `OUTFMT` to generate visualization for aligner output.
```bash
make -f src/alignment.mk view \
    OUTFMT=clustal OUTNAME=my_alignment
```
