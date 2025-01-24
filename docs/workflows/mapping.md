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
- [minimap2](doi:10.1093/bioinformatics/bty191) (for long reads)

All alignment files are sorted and converted into BAM format by default.

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

### bowtie2.mk

Map reads against a reference using `bowtie2`. 

The reference file is first indexed prior to mapping. By default, the output alignment is sorted and converted in BAM format.

**{sc}`Parameters`**

- R1: first set of reads for pair-end FASTQ file.
- R2: second set of reads for pair-end FASTQ file, if any.
- REF: path to FASTA reference file.
- BAM: path to output BAM file.
- STATS: path to alignment statistics file (optional). 
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Map pair-end reads against a reference using `bowtie2`:
```bash
make -f bowtie2.mk R1=SRR818231_1.fastq R2=SRR818231_2.fastq REF=refs/genome.fa align
```

Map single-end reads against a reference:
```bash
make -f bowtie2.mk R1=SRR818231.fastq REF=refs/genome.fa align
```

Generate alignment statistics for BAM file:
```bash
make -f bowtie2.mk R1=SRR818231.fastq REF=refs/genome.fa stats
```


### bwa.mk

Map reads against a reference using `bwa`. 

**{sc}`Parameters`**

- R1: first set of reads for pair-end FASTQ file.
- R2: second set of reads for pair-end FASTQ file, if any.
- REF: path to FASTA reference file.
- BAM: path to output BAM file.
- STATS: path to alignment statistics file (optional). 
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Map pair-end reads against a reference using `bwa`:
```bash
make -f bwa.mk R1=SRR818231_1.fastq R2=SRR818231_2.fastq REF=refs/genome.fa align
```

Map single-end reads against a reference:
```bash
make -f bwa.mk R1=SRR818231.fastq REF=refs/genome.fa align
```

Generate alignment statistics for BAM file:
```bash
make -f bwa.mk R1=SRR818231.fastq REF=refs/genome.fa stats
```

### minimap2.mk

Align long reads against a reference using `minimap2`.

Note that `minimap2` does not allow pair-end reads as input. A single FASTQ file must be used for alignment.

**{sc}`Parameters`**

- FQ: path to a FASTQ file.
- REF: path to FASTA reference file.
- BAM: path to output BAM file (optional). 
- ACC: an accession used to construct default filenames (optional).
- STATS: path to alignment statistics file (optional). 
- THREADS: number of cores (default: 4)

**{sc}`Example Usage`**

Map long reads against a reference using `minimap2`:
```bash
make -f minimap2.mk FQ=SRR31340505.fastq REF=ON963982.fa align
```

Generate alignment statistics:
```bash
make -f minimap2.mk FQ=SRR31340505.fastq REF=ON963982.fa stats
```
