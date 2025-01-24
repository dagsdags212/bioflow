---
downloads:
  - file: ../../src/alignment/clustalo.mk
    title: clustalo.mk
  - file: ../../src/alignment/mafft.mk
    title: mafft.mk
  - file: ../../src/alignment/muscle.mk
    title: muscle.mk
  - file: ../../envs/bf-alignment.yml
    title: env.yml
---

(bf-alignment)=
# alignment

## Overview

The `alignment` module contains rules for performing pairwise and multiple sequence alignment.

Currently supported alignment tools include:

- [Clustal-Omega](doi:10.1038/msb.2011.75)
- [MAFFT](doi:10.1093/nar/gkf436)
- [MUSCLE](doi:10.1186/1471-2105-5-113)

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

### clustalo.mk

Perform multiple sequence alignment with `clustal-omega`.

**{sc}`Parameters`**

- `FA`: a FASTA file containing sequences to align.
- `FORMAT`: output MSA format (default: fasta).
- `MSA`: specify path of output MSA file (optional).
- `THREADS`: number of worker threads.

**{sc}`Example Usage`**

Align a FASTA file with `clustlo`:
```bash
make -f clustalo.mk FA=seqs.fa THREADS=8 align
```

Output the alignment file in clustal format:
```bash
make -f clustalo.mk FA=seqs.fa THREADS=8 FORMAT=clustal align
```

### mafft.mk

Perform multiple sequence alignment with `mafft`.

**{sc}`Parameters`**

- `FA`: a FASTA file containing sequences to align.
- `FORMAT`: output MSA format (default: fasta).
- `METHOD`: alignment method to conduct (default: local).
- `MSA`: specify path of output MSA file (optional).
- `THREADS`: number of worker threads.

**{sc}`Example Usage`**

Perform global alignment with `mafft`:
```bash
make -f mafft.mk FA=seqs.fa METHOD=global align
```

Perform local alignment with `mafft`:
```bash
make -f mafft.mk FA=seqs.fa METHOD=local align
```

Output the alignment file in clustal format:
```bash
make -f mafft.mk FA=seqs.fa THREADS=8 FORMAT=clustal align
```

### muscle.mk

Perform multiple sequence alignment with `muscle`.

**{sc}`Parameters`**

- `FA`: a FASTA file containing sequences to align.
- `MSA`: specify path of output MSA file (optional).
- `THREADS`: number of worker threads.

**{sc}`Example Usage`**

Align a FASTA file with `Muscle`:
```bash
make -f clustalo.mk FA=seqs.fa THREADS=8 align
```

:::{tip} TODO
Structural alignment with muscle.
:::
