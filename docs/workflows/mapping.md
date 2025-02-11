---
downloads:
  - file: ../../src/mapping/bowtie2.mk
    title: bowtie2.mk
  - file: ../../src/mapping/bwa.mk
    title: bwa.mk
  - file: ../../src/mapping/minimap2.mk
    title: minimap2.mk
  - file: ../../envs/bf-mapping.yml
    title: env.yml
---

(bf-mapping)=
# mapping

## Overview

The `mapping` module contains rules for performing read alignment against a reference genome.

Currently supported aligners include:

- [bwa](doi:10.1093/bioinformatics/btp698)
- [bowtie2](doi:10.1038/nmeth.1923)
- [minimap2](doi:10.1093/bioinformatics/bty191)

All alignment files are sorted and converted into BAM format by default.

:::{hint} Environment Setup
:class: dropdown

Initialize the micromamba environment and install the dependencies:
```{code-cell} bash
bf-<recipe> install
```

where `<recipe>` is the name of the specific recipe.

:::

## Recipes

### bf-bowtie2

Map reads against a reference using `bowtie2`. 

The reference file is first indexed prior to mapping. By default, the output alignment is sorted and converted in BAM format.

**{sc}`Parameters`**

- `REF`: a reference file in FASTA format.
- `R1`: first of read-pair in FASTQ format.
- `R2`: second of read-pair in FASTQ format.
- `DIR`: a directory path for storing BAM files (default: bam).
- `COVPLOT`: a filepath for saving the coverage plot.
- `THREADS`: number of cores (default: 4)

**{sc}`Example Usage`**

Map pair-end reads against a reference using `bowtie2`:
```{code-cell} bash
bf-bowtie2 R1=SRR818231_1.fastq R2=SRR818231_2.fastq REF=refs/genome.fa align
```

Map single-end reads against a reference:
```{code-cell} bash
bf-bowtie2 R1=SRR818231.fastq REF=refs/genome.fa align
```

Generate alignment statistics for BAM file:
```{code-cell} bash
bf-bowtie2 R1=SRR818231.fastq REF=refs/genome.fa stats
```

Plot the read coverage:
```{code-cell} bash
bf-bowtie2 R1=SRR818231.fastq REF=refs/genome.fa COVPLOT=plots/cov.png coverage
```


### bf-bwa

Map reads against a reference using `bwa`. 

**{sc}`Parameters`**

- `REF`: a reference file in FASTA format.
- `R1`: first of read-pair in FASTQ format.
- `R2`: second of read-pair in FASTQ format.
- `DIR`: a directory path for storing BAM files (default: bam).
- `COVPLOT`: a filepath for saving the coverage plot.
- `THREADS`: number of cores (default: 4)

**{sc}`Example Usage`**

Map pair-end reads against a reference using `bwa`:
```{code-cell} bash
bf-bwa R1=SRR818231_1.fastq R2=SRR818231_2.fastq REF=refs/genome.fa align
```

Map single-end reads against a reference:
```{code-cell} bash
bf-bwa R1=SRR818231.fastq REF=refs/genome.fa align
```

Generate alignment statistics for BAM file:
```{code-cell} bash
bf-bwa R1=SRR818231.fastq REF=refs/genome.fa stats
```

Index the output BAM file:
```{code-cell} bash
bf-bwa R1=SRR818231.fastq REF=refs/genome.fa index
```

Plot the read coverage:
```{code-cell} bash
bf-bwa R1=SRR818231.fastq REF=refs/genome.fa COVPLOT=plots/cov.png coverage
```

### bf-minimap2

Align long reads against a reference using `minimap2`.

Note that `minimap2` does not allow pair-end reads as input. A single FASTQ file must be used for alignment.

**{sc}`Parameters`**

- `R1`: a read file in FASTQ format.
- `REF`: a reference file in FASTA format.
- `DIR`: a directory path for storing BAM files (default: bam).
- `COVPLOT`: a filepath for saving the coverage plot.
- `THREADS`: number of cores (default: 4)

**{sc}`Example Usage`**

Map long reads against a reference using `minimap2`:
```{code-cell} bash
bf-minimap2 R1=SRR31340505.fastq REF=ON963982.fa align
```

Generate alignment statistics:
```{code-cell} bash
bf-minimap2 R1=SRR31340505.fastq REF=ON963982.fa stats
```

Plot the read coverage:
```{code-cell} bash
bf-minimap2 R1=SRR31340505.fastq REF=ON963982.fa COVPLOT=plots/cov.png coverage
```
