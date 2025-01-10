---
downloads:
  - file: ../../src/fetch.mk
    title: Makefile
  - file: ../../envs/fetch.yml
    title: env.yml
---

(bwf-fetch)=
# fetch.mk

The `fetch.mk` workflow can be used to download different types of biological data from online databases. Currently support data formats are listed in @supported-formats:

:::{table} Support data formats by `fetch.mk`
:label: supported-formats
:align: center

| Data Type | Format | Command |
| ------ | ------ | ------ |
| Sequencing reads | FASTQ | `sra` | 
| Reference genomes | FASTA | `ref` |
| Protein structures | PDB | `pdb` |

:::

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
make -f src/fetch.mk init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bwf-fetch
```
:::

## sra

Retrieve a set of sequencing reads from a project ID (PRJNA) or a single sequencing run (SRR). 

All reads are stored in the `reads` directory. When downloading multiple sets of reads, a subdirectory for each set labeled after its SRR accession is created under `reads`.

**{sc}`Parameters`**

- PRJNA: project identifier
- SRR: sequencing run identifier
- X: number of spots to download
- PE: download pair-end reads (default: `true`)

**{sc}`Example Usage`**

Download a set of complete sequencing reads.
```bash
make -f src/fetch.mk PRJNA=PRJNA1066786
```

Download 100000 reads from a single run.
```bash
make -f src/fetch.mk SRR=SRR27644850 X=100000
```

## ref

**{sc}`Parameters`**

- ACC: accession identifer of a nucleotide sequence
- INCLUDE_GFF: if `true`, downloads the annotation file in GFF3 format

**{sc}`Example Usage`**

Download the canonical reference genome of the African Swine Fever virus.
```bash
make -f src/fetch.mk ACC=GCF_003047755.2
```

Include the annotation file.
```bash
make -f src/fetch.mk ACC=GCF_003047755.2 INCLUDE_GFF=true
```

## pdb

Retrieve a structure file from the Protein Data Bank using its PDB ID.

```{note}
PDB identifiers are four-character alphanumerics such as _2hbs_.
```

**{sc}`Parameters`**

- PDB: a PDB identifier

**{sc}`Example Usage`**

Download the SARS-CoV-2 spike glycoprotein with PDB ID `7FCD`.
```bash
make -f src/fetch.mk PDB=7FCD
```
