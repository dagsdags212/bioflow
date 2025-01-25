---
downloads:
  - file: ../../src/phylo/iqtree.mk
    title: iqtree.mk
  - file: ../../src/phylo/raxml-ng.mk
    title: raxml-ng.mk
  - file: ../../envs/bf-phylo.yml
    title: env.yml
---

(bf-phylo)=
# phylo

## Overview

The `phylo` module contains rules for performing phylogenetic tree inference using `raxml-ng` and `iqtree`.

Inference is conducted using either [maximum likelihood](wiki:Maximum_likelihood_estimation) or [bootstrapping](wiki:Bootstrapping_(statistics)) methods. To generate the alignment, refer to the [alignment workflow](./alignment.md).

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

### iqtree.mk

Generate trees with `iqtree`.

**{sc}`Parameters`**

- ALN: a multiple alignment file.
- MODEL: substitution matrix or evolutionary model to use (default: GTR+G).
- SEED: set initial seed.
- OUTDIR: path to a directory for storing output (default: iqtree).
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Perform basic inference with `iqtree`.
```bash
make -f iqtree.mk ALN=seq.aln run
```

### raxml-ng.mk

Generate trees with `raxml-ng`.

**{sc}`Parameters`**

- ALN: a multiple alignment file.
- MODEL: substitution matrix or evolutionary model to use (default: GTR+G).
- METHOD: tree construction method (default: ml).
- N: number of ML searches to perform (default: 10)
- SEED: set initial seed.
- OUTDIR: path to a directory for storing output (default: raxml).
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Perform basic inference with `iqtree`.
```bash
make -f raxml-ng.mk ALN=seq.aln run
```
