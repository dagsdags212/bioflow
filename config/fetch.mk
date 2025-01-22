# Conda environment
ENV := bf-fetch

# Check if dependencies are installed
dependencies := entrez-direct sra-tools ncbi-datasets-cli pdb-tools

# === Parameters for the `sra` command ===
# Sequencing project accession.
PRJN ?=

# Sequecing read accession.
SRX ?=

# Number of spots to retrieve.
X ?=

# If true, separate into forward and reverse reads
PE ?=

# === Parameters for the `ref` command ===
ACC ?= 

# If true, include annotation file for the assembly
GFF ?= false

# If true, include genbank record for the sequence file
GB ?= false

# === Parameters for the `pdb` command ===
# Protein structure ID
PDB ?= 

# === Parameters for the `pubmed` command ===
# Query string of pubmed search
QUERY ?=
