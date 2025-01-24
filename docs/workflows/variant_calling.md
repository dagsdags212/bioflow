---
downloads:
  - file: ../../src/variant_calling/bcftools.mk
    title: bcftools.mk
  - file: ../../src/variant_calling/freebayes.mk
    title: freebayes.mk
  - file: ../../envs/bf-vc.yml
    title: env.yml
---

(bf-vc)=
# variant_calling

## Overview

The **variant_calling** module contains recipes for calling SNPs using `bcftools` and `freebayes`. Both recipes expect a BAM alignment file and a reference genome in FASTA format.

Invoking the `run` command generates three outputs: (1) a VCF file, (2) an index to the VCF file and (3) a flatfile containing summary statistics for the called variants. The last file is generated from running `bcftools stats` on the output VCF file.

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

### freebayes.mk

Generate SNP calls with `freebayes`.

This commands expects a reference file in FASTA format and an alignment file in BAM format. Both reference and alignment files are indexed prior to calling variants.

**{sc}`Parameters`**

- ACC: Genbank accession for the reference file.
- BAM: an alignment file in BAM format.
- SRR: an SRA accession to derive the BAM filename (optional).
- THREADS: number of cores (default: 2)

**{sc}`Example Usage`**

Call variants from a BAM file:
```bash
make -f freebayes.mk ACC=AF086833 BAM=aln.bam run
```

Specify the path to a reference file instead of providing an accession:
```bash
make -f freebayes.mk ACC=AF086833 REF=refs/genome.fa run
```

### bcftools.mk

Generate SNP calls with `bcftools`.

**{sc}`Parameters`**

- ACC: Genbank accession for the reference file.
- BAM: an alignment file in BAM format.
- SRR: an SRA accession to derive the BAM filename (optional).
- THREADS: number of cores (default: 2)

**{sc}`Example Usage`**

Call variants from a BAM file:
```bash
make -f bcftools.mk ACC=AF086833 BAM=aln.bam run
```

Specify the path to a reference file instead of providing an accession:
```bash
make -f bcftools.mk ACC=AF086833 REF=refs/genome.fa run
```

### snpeff.mk

Annotate called variants with `snpeff`.

:::{tip} In Progress
...
:::
