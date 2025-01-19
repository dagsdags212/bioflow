---
downloads:
  - file: ../../src/fetch/Makefile
    title: Makefile
  - file: ../../envs/bf-fetch.yml
    title: env.yml
---

(bf-fetch)=
# fetch.mk

## Overview

The `fetch` module can be used to download different types of biological data from online databases. Its entry point is the `bffetch` command. Currently supported data formats are listed in @supported-formats:

:::{table} Supported data formats by `fetch.mk`
:label: supported-formats
:align: center

| Data Type | Source | Format | Command |
| ------ | ------ | ------ | ------ |
| Sample reads | SRA | FASTQ | `reads` | 
| Metagenomic reads | MGnify/SRA | FASTQ | `reads` | 
| Reference genomes | NCBI | FASTA | `seq` |
| Sequence records | NCBI |  genbank/gb | `seq` |
| Protein structures | RCSB | PDB | `pdb` |
| Journal metadata | PubMed | text | `pubmed` |

:::

:::{hint} Environment Setup
:class: dropdown

Prior to using the workflow, download the dependencies within a virtual environment using your manager of choice:

```bash
bffetch init ENV_MANAGER=micromamba
```

Activate environment to expose dependencies:
```bash
micromamba activate bf-fetch
```
:::

## Commands

### reads

Retrieve a set of sequencing reads from a project ID (PRJNA) or a single sequencing run (SRR).

Downloading metagenomic reads using a MGnify identifier is supported. Behind the scenes, the MGnify-formated accession if converted into an SRA/ENA accession which is then used to download the reads from NCBI.

All reads are stored in the `data` directory. When downloading multiple sets of reads, a subdirectory for each set labeled after its SRR accession is created under `data`.

**{sc}`Parameters`**

- PRJNA: project identifier
- SRR: sequencing run identifier
- X: number of spots to download
- PE: download pair-end reads (default: `true`)

**{sc}`Example Usage`**

Download a set of complete sequencing reads:
```bash
bffetch reads PRJNA=PRJNA1066786
```

Download 100000 reads from a single run:
```bash
bffetch reads SRR=SRR27644850 X=100000
```

Download reads derived from a metagenomic sample stored in MGnify:
```bash
bffetch reads PRJNA=MGYS00000259 X=10000
```

### seq

Retrieve the sequence file of a gene or assembly. For genomes, associated annotation files in GFF3 format can be included in the downloadeby setting `INCLUDE_GFF` to `true`. Likewise, the GenBank record for genes can be retrieved by setting `INCLUDE_GB` to `true`.

**{sc}`Parameters`**

- ACC: accession identifer of a nucleotide sequence
- INCLUDE_GFF: if `true`, downloads the associated annotation file in GFF3 format (default: false)
- INCLUDE_GB: if `true`, downloads the associated GenBank record (default: false)

**{sc}`Example Usage`**

Download an assembly (GCF/GCA) and include its annotation file (if it exists):
```bash
bffetch seq ACC=GCF_003047755.2 INCLUDE_GFF=true
```

Download a single gene and include its GenBank record:
```bash
bffetch seq ACC=887105 INCLUDE_GFF=true
```

### pdb

Retrieve a structure file from the Protein Data Bank.

```{note}
PDB identifiers are four-character alphanumerics such as _2HBS_. By conventional, alphabetic characters are in uppercase.
```

**{sc}`Parameters`**

- PDB: a PDB identifier

**{sc}`Example Usage`**

Download the SARS-CoV-2 spike glycoprotein with PDB ID `7FCD`.
```bash
bffetch pdb PDB=7FCD
```

### pubmed

Retrieve a list of PubMed articles from a query string. 

By default, the query results are printed to the standard output. Use the redirect operator (`>`) to save the output to a text file.

**{sc}`Parameters`**

- QUERY: a search term used to query PubMed

**{sc}`Example Usage`**

Search for a list of articles on ASFV assemblies:
```bash
bffetch pubmed QUERY="African swine fever virus assemblies"
```

Save the results to a text file:
```bash
bffetch pubmed QUERY="African swine fever virus assemblies" > asfv_assemblies.journals.txt
```
