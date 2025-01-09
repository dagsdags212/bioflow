#!/usr/bin/env bash

set -eux

# Accession for ASFV reference genome
REF=GCF_003047755.2

# Reads from a sequencing experiment on ASFV
SRR=SRR27644850

# Download ASFV reference genome
micromamba run -n bwf-fetch \
  make -f ../src/fetch.mk ref ACC=${REF} INCLUDE_GFF=true

# Download sequencing reads
micromamba run -n bwf-fetch \
  make -f ../src/fetch.mk sra SRR=${SRR} PE=true X=100000
