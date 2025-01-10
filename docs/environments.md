---
downloads:
  - file: ../envs/assembly.yml
    title: assembly.yml
  - file: ../envs/fetch.yml
    title: fetch.yml
  - file: ../envs/qc.yml
    title: qc.yml
---
# Environments

::::{grid} 1 1 2 3

:::{card}
:header: assembly.mk
```yaml
name: bwf-assembly
channels:
- bioconda
dependencies:
- bandage
- megahit
- quast
- spades
```
:::

:::{card}
:header: fetch.mk
```yaml
name: bwf-fetch
channels:
- bioconda
- conda-forge
dependencies:
- entrez-direct
- ncbi-datasets-cli
- sra-tools
```
:::

:::{card}
:header: mapping.mk
```yaml
name: bwf-mapping
channels:
- bioconda
dependencies:
- bowtie2
- bwa
- qualimap
- samtools
```
:::

:::{card}
:header: qc.mk
```yaml
name: bwf-qc
channels:
- bioconda
dependencies:
- cutadapt
- fastp
- fastqc
- multiqc
- trimmomatic
```
:::

::::

