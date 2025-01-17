---
downloads:
  - file: ../../src/fetch.mk
    title: Makefile
  - file: ../../envs/fetch.yml
    title: env.yml
---

(bf-fetch)=
# fetch.mk

## Overview

The `fetch.mk` workflow can be used to download different types of biological data from online databases. Its entry point is the `bf-fetch` command. Currently supported data formats are listed in @supported-formats:

:::{table} Supported data formats by `fetch.mk`
:label: supported-formats
:align: center

| Data Type | Source | Format | Command |
| ------ | ------ | ------ | ------ |
| Sample reads | SRA | FASTQ | `sra` | 
| Metagenomic reads | MGnify/SRA | FASTQ | `sra` | 
| Reference genomes | NCBI | FASTA | `ref` |
| Protein structures | RCSB | PDB | `pdb` |
| Journal metadata | PubMed | text | `pubmed` |
| Sequence records | NCBI |  genbank/gb | `genbank` |

:::

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bf-fetch init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-fetch
```
:::

## Rules

### sra

Retrieve a set of sequencing reads from a project ID (PRJNA) or a single sequencing run (SRR).

Downloading metagenomic reads using a MGnify identifier is supported. Behind the scenes, a script converts the MGnify-formated accession into an SRA/ENA accession which is then use to fetch the reads from NCBI.

All reads are stored in the `reads` directory. When downloading multiple sets of reads, a subdirectory for each set labeled after its SRR accession is created under `reads`.

**{sc}`Parameters`**

- PRJNA: project identifier
- SRR: sequencing run identifier
- X: number of spots to download
- PE: download pair-end reads (default: `true`)

**{sc}`Example Usage`**

Download a set of complete sequencing reads.
```bash
bf-fetch sra PRJNA=PRJNA1066786
```

Download 100000 reads from a single run.
```bash
bf-fetch sra SRR=SRR27644850 X=100000
```

Download reads derived from a metagenomic sample stored in MGnify.
```bash
bf-fetch sra PRJNA=MGYS00000259 X=10000
```

### ref

**{sc}`Parameters`**

- ACC: accession identifer of a nucleotide sequence
- INCLUDE_GFF: if `true`, downloads the annotation file in GFF3 format

**{sc}`Example Usage`**

Download the canonical reference genome of the African Swine Fever virus.
```bash
bf-fetch ref ACC=GCF_003047755.2
```

Include the annotation file.
```bash
bf-fetch ref ACC=GCF_003047755.2 INCLUDE_GFF=true
```

### pdb

Retrieve a structure file from the Protein Data Bank using its PDB ID.

```{note}
PDB identifiers are four-character alphanumerics such as _2hbs_.
```

**{sc}`Parameters`**

- PDB: a PDB identifier

**{sc}`Example Usage`**

Download the SARS-CoV-2 spike glycoprotein with PDB ID `7FCD`.
```bash
bf-fetch pdb PDB=7FCD
```

### pubmed

Retrieve a list of PubMed articles from a query string. 

By default, the query results are printed to the standard output. Use the redirect operator (`>`) to save the results to a text file.

**{sc}`Parameters`**

- QUERY: a search term used to query PubMed

**{sc}`Example Usage`**

Search for a list of articles on ASFV assemblies:
```bash
bf-fetch pubmed QUERY="African swine fever virus assemblies"
```

Save the results to a text file:
```bash
bf-fetch pubmed QUERY="African swine fever virus assemblies" > asfv_assemblies.journals.txt
```

### genbank

Retrieve a Genbank record from an accession ID.

Similar to the output of `pubmed`, query results are printed to stdout. Use the redirect operator (`>`) to save the results to a text file.

**{sc}`Parameters`**

- ACC: accession identifer of a nucleotide sequence

**{sc}`Example Usage`**

Fetch the Genbank record of the Wuhan isolate of SARS-CoV-2.
```bash
bf-fetch genbank ACC=NC_045512
```
