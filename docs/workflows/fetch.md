---
downloads:
  - file: ../../src/fetch/aria.mk
    title: aria.mk
  - file: ../../src/fetch/sra.mk
    title: sra.mk
  - file: ../../src/fetch/genbank.mk
    title: genbank.mk
  - file: ../../src/fetch/pubmed.mk
    title: pubmed.mk
  - file: ../../src/fetch/pdb.mk
    title: pdb.mk
  - file: ../../envs/bf-fetch.yml
    title: bf-fetch.yml
---

(bf-fetch)=
# fetch

## Overview

The `fetch` module contains recipes for downloading biological data of various formats. Currently supported data formats are listed in @supported-formats.

:::{table} Supported data formats by the `fetch` module.
:label: supported-formats
:align: center

| Data Type | Source | Format | Makefile |
| ------ | ------ | ------ | ------ |
| Sample reads | SRA | FASTQ | `sra.mk` | 
| Metagenomic reads | MGnify/SRA | FASTQ | `sra.mk` | 
| Sequence files | NCBI | FASTA/GenBank | `genbank.mk` |
| Annotation files | NCBI | GFF | `genbank.mk` |
| Protein structures | RCSB | PDB | `pdb.mk` |
| Journal metadata | PubMed | text | `pubmed.mk` |

:::

:::{hint} Environment Setup
:class: dropdown

Initialize the micromamba environment and install the dependencies:
```{code-cell} bash
make -f $BIOFLOW/src/fetch/<recipe>.mk install
```

where `<recipe>` is the makefile for the specific recipe.

:::

## Recipes

### sra.mk

Retrieve FASTQ reads from the SRA using an accession identifier (SRR).

**{sc}`Parameters`**

- `SRR`: an SRA accession
- `N`: number of reads to download; set to `ALL` to download all reads (default: ALL)
- `DIR`: output directory for storing reads (default: reads)
- `MODE`: download in pair-end (PE) or single-end (SE) mode (default: SE)

**{sc}`Example Usage`**

Retrieve the complete set of reads for **SRR1554325**:
```{code-cell} bash
make -f sra.mk SRR=SRR1554325 run
```

Download a subset of 10,000 reads for **SRR1554325**, specifying the output directory:
```{code-cell} bash
make -f sra.mk SRR=SRR1554325 N=10000 DIR=reads run
```

### aria.mk

Download data from a URL using the [aria2c](https://aria2.github.io/) command line utility.

**{sc}`Parameters`**

- `URL`: an HTTP/FTP link to download.
- `DIR`: output directory for the data (default: data). 

**{sc}`Example Usage`**

Download the bowtie index for the human reference genome (GRCh38):
```{code-cell} bash
URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz
make -f aria.mk URL=$URL run
```

### genbank.mk

Retrieve sequence files from GenBank. This recipe is used to fetch reference sequences, GenBank records, and annotation files.

**{sc}`Parameters`**

- `ACC`: an NCBI accession.
- `REF`: output path for sequence file (default: refs/<ACC>.fa).
- `GFF`: output path for annotation file (default: refs/<ACC>.gff).
- `GBK`: output path for annotation file (default: refs/<ACC>.gb)

**{sc}`Example Usage`**

Retrieve the sequence file for **AF086833**.
```{code-cell} bash
make -f genbank.mk ACC=AF086833 fasta
```

Retrieve the GenBank record for **AF086833**.
```{code-cell} bash
make -f genbank.mk ACC=AF086833 genbank
```

Retrieve the annotation file for **AF086833**.
```{code-cell} bash
make -f genbank.mk ACC=AF086833 gff
```

Retrieve the sequence, annotation, and GenBank files for **AF086833** in one command.
```{code-cell} bash
make -f genbank.mk ACC=AF086833 all
```

### pubmed.mk

Retrieve journal metadata by querying the PubMed literature database.

Given a [PMID or PMCID](https://nexus.od.nih.gov/all/2015/08/31/pmid-vs-pmcid-whats-the-difference/), a user can retrieve journal-related information such as its abstract, medline record, APA citation, and a list of cited articles. Multiple PMIDs can also be passed in a single invocation, just make to delimit each identifier with a space and enclose them within parentheses. See below for some examples. 

By default, query results are printed to the standard output of your terminal. Use the redirect operator (`>`) to save the results to a designated file.

**{sc}`Parameters`**

- `PMID`: a Pubmed accession.

**{sc}`Example Usage`**

Retrieve the abstract for the **Bowtie2** journal (22388286):
```{code-cell} bash
make -f pubmed.mk PMID=22388286 abstract
```

:::{seealso} See Output
:class: dropdown
```
1. Nat Methods. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923.

Fast gapped-read alignment with Bowtie 2.

Langmead B(1), Salzberg SL.

Author information:
(1)Center for Bioinformatics and Computational Biology, Institute for Advanced 
Computer Studies, University of Maryland, College Park, Maryland, USA. 
blangmea@jhsph.edu

As the rate of sequencing increases, greater throughput is demanded from read 
aligners. The full-text minute index is often used to make alignment very fast 
and memory-efficient, but the approach is ill-suited to finding longer, gapped 
alignments. Bowtie 2 combines the strengths of the full-text minute index with 
the flexibility and speed of hardware-accelerated dynamic programming algorithms 
to achieve a combination of high speed, sensitivity and accuracy.

DOI: 10.1038/nmeth.1923
PMCID: PMC3322381
PMID: 22388286 [Indexed for MEDLINE]
```
:::

Retrieve the APA citation for a couple of journals:
```{code-cell} bash
make -f pubmed.mk PMID="22388286 19451168" apa
```

:::{seealso} See Output
:class: dropdown
```
Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics Oxford England, 25(14), 1754-1760. doi:10.1093/bioinformatics/btp324
Langmead, B., & Salzberg, S. L. (2012). Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4), 357-359. doi:10.1038/nmeth.1923
```
:::

Get a complete list of cited articles for the **BWA** journal (19451168) and save to a file:
```{code-cell} bash
make -f pubmed.mk PMID=19451168 citations > references.txt
```

:::{seealso} See Output
:class: dropdown
```
- Burrows M, Wheeler DJ. Technical report 124. Palo Alto, CA: Digital Equipment Corporation; 1994. A block-sorting lossless data compression algorithm.
- Campagna D, et al. PASS: a program to align short sequences. Bioinformatics. 2009;25:967–968.
- Eaves HL, Gao Y. MOM: maximum oligonucleotide mapping. Bioinformatics. 2009;25:969–970.
- Ferragina P, Manzini G. Proceedings of the 41st Symposium on Foundations of Computer Science (FOCS 2000) IEEE Computer Society; 2000. Opportunistic data structures with applications; pp. 390–398.
- Grossi R, Vitter JS. Proceedings on 32nd Annual ACM Symposium on Theory of Computing (STOC 2000) ACM; 2000. Compressed suffix arrays and suffix trees with applications to text indexing and string matching; pp. 397–406.
- Hon W.-K, et al. A space and time efficient algorithm for constructing compressed suffix arrays. Algorithmica. 2007;48:23–36.
- Jiang H, Wong WH. SeqMap: mapping massive amount of oligonucleotides to the genome. Bioinformatics. 2008;24:2395–2396.
- Jung Kim Y, et al. ProbeMatch: a tool for aligning oligonucleotide sequences. Bioinformatics. 2009;25:1424–1425.
- Lam TW, et al. Compressed indexing and local alignment of DNA. Bioinformatics. 2008;24:791–797.
- Langmead B, et al. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol. 2009;10:R25.
- Lin H, et al. ZOOM! Zillions of oligos mapped. Bioinformatics. 2008;24:2431–2437.
- Lippert RA. Space-efficient whole genome comparisons with Burrows-Wheeler transforms. J. Comput. Biol. 2005;12:407–415.
- Li H, et al. Mapping short DNA sequencing reads and calling variants using mapping quality scores. Genome Res. 2008a;18:1851–1858.
- Li R, et al. SOAP: short oligonucleotide alignment program. Bioinformatics. 2008b;24:713–714.
- Malhis N, et al. Slider–maximum use of probability information for alignment of short sequence reads and SNP detection. Bioinformatics. 2009;25:6–13.
- Schatz M. Cloudburst: highly sensitive read mapping with mapreduce. Bioinformatics. 2009;25:1363–1369.
- Smith AD, et al. Using quality scores and longer reads improves accuracy of Solexa read mapping. BMC Bioinformatics. 2008;9:128.
- McKenna A, et al. Genome Res. 2010;20:1297–1303.
- Trapnell C, et al. Nat. Biotechnol. 2010;28:511–515.
- Langmead B, Hansen KD, Leek JT. Genome Biol. 2010;11:R83.
- Li H, Homer N. Brief. Bioinform. 2010;11:473–483.
- Ferragina P, Manzini G. Proc. 41st Annual Symposium on Foundations of Computer Science. IEEE Comput. Soc.; 2000. pp. 390–398.
- Langmead B, Trapnell C, Pop M, Salzberg SL. Genome Biol. 2009;10:R25.
- Lam T, et al. IEEE International Conference on Bioinformatics and Biomedicine. 2009:31–36.
- Li H, Durbin R. Bioinformatics. 2009;25:1754–1760.
- Li H, Durbin R. Bioinformatics. 2010;26:589–595.
- Li R, et al. Bioinformatics. 2009;25:1966–1967.
- Ajay SS, Parker SC, Ozel Abaan H, Fuentes Fajardo KV, Margulies EH. Genome Res. 2010;21:1498–1505.
- 1000 Genomes Project Consortium Nature. 2010;467:1061–1073.
- Rothberg JM, et al. Nature. 2011;475:348–352.
- Li H, et al. Bioinformatics. 2009;25:2078–2079.
```
:::

### pdb.mk

Retrieve a structure file from the Protein Data Bank.

```{note}
PDB identifiers are four-character alphanumerics such as _2HBS_. By conventional, alphabetic characters are in uppercase.
```

**{sc}`Parameters`**

- `ID`: a 4-character PDB identifier.

**{sc}`Example Usage`**

Download the SARS-CoV-2 spike glycoprotein with PDB ID `7FCD`.
```{code-cell} bash
make -f pdb.mk ID=7FCD
```
