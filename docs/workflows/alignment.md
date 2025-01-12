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

The `align` command generates three output files containing the same alignment data but in different formats (fasta, clustal, phylip).

**{sc}`Parameters`**

- FA: path to a multi-sequence FASTA file
- DIR: path to a directory containing multiple FASTA files
- ALIGNER: specify tool for performing alignment (default: mafft)
- GOP: gap opening penalty (default: 1.53)
- GEP: gap extension penalty (default: 0.0)
- ITER: maximum number of iteration for refinement (default: 1)
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

Specify output filename to generate _aln.fasta_, _aln.phylip_, and _aln.clustal_.
```bash
make -f src/alignment.mk align \
    FA=seqs.fa ALIGNER=mafft \
    OUTNAME=aln
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

Specify the `OUTNAME` parameter to only render the alignment files with the given prefix. The output is stored in the same directory in HTML format.

**{sc}`Parameters`**

- OUTNAME: filename of output alignment file (default: aln)

**{sc}`Example Usage`**

Specify `OUTNAME` and `OUTFMT` to generate visualization for aligner output.
```bash
make -f src/alignment.mk view OUTNAME=aln
```
