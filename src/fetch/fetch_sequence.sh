#!/usr/bin/env bash

usage() {
  echo "Retrieve a sequence file from NCBI"
  echo
  echo "Usage:"
  echo "  fetch_sequence.sh [-h] [-v] <acc>"
  echo
  echo "  <acc>  NCBI sequence accession"
  echo
  echo "Options:"
  echo "  -h  display this help message"
  echo "  -g  include the GFF file"
  echo "  -b  include GenBank record"
  echo "  -v  print tool version"
  exit
}

validate_accession() {
  local acc=$1
  local gene_regex="^[0-9]{8}$"
  local protein_regex="^[A-Z0-9]{4}$"
  local assembly_regex="^GC[AF]_[0-9]{9}\.[0-9]$"
  local dna_regex="^[AMNPX][ACFGTWZMR]_?[0-9]{6,9}\.[0-9]$"

  if [[ "${acc}" =~ ${dna_regex} ]]; then
    echo "Detected nuccore accession"
    declare -g ACC=${acc}
    declare -g SEQTYPE=dna
  elif [[ "${acc}" =~ ${gene_regex} ]]; then
    echo "Detected gene accession"
    declare -g ACC=${acc}
    declare -g SEQTYPE=gene
  elif [[ "${acc}" =~ ${protein_regex} ]]; then
    echo "Detected protein accession"
    declare -g ACC=${acc}
    declare -g SEQTYPE=protein
  elif [[ "${acc}" =~ ${assembly_regex} ]]; then
    echo "Detected assembly accession"
    declare -g ACC=${acc}
    declare -g SEQTYPE=assembly
  else
    echo "Error: provided accession is invalid"
    exit 1
  fi
}

download_dna() {
  local acc=$1
  local gb=$2

  echo "Downloading sequence file for ${acc}" 1>&2
  efetch -db nuccore -id ${acc} -format fasta >data/${acc%%.*}.fa && echo "Done!" 1>&2
  if [[ "${gb}" == true ]]; then
    echo "Including GenBank record"
    efetch -db nuccore -id ${acc} -format genbank >data/${acc%%.*}.gb && echo "Done!" 1>&2
  fi
}

download_gene_info() {
  local acc=$1
  echo "Retrieving genomic accession of ${acc}"
  gene_acc=$(efetch -db gene -id ${acc} -format tabular | cut -f12 | tail -n1)
  echo "Accession found: ${acc} > ${gene_acc}"
  echo "Downloading FASTA file for ${gene_acc}"
  efetch -db nuccore -id ${gene_acc} -format fasta >data/${acc%%.*}.fa && echo "Done!" 1>&2
}

download_assembly() {
  local acc=$1
  local gff=$2

  echo "Downloading assembly file for ${acc}"
  if [[ "${gff}" == true ]]; then
    echo "Including annotation file in GFF3 format"
    datasets download genome accession ${acc} --include genome,gff3 --filename ${acc%%.*}.zip
  else
    datasets download genome accession ${acc} --include genome --filename ${acc%%.*}.zip
  fi
  unzip ${acc%%.*}.zip -d ${acc%%.*}
  rm -f ${acc%%.*}.zip
  echo "Done!"
}

download_sequence() {
  local acc=$1
  local gff=$2
  local gb=$3
  mkdir -p data/

  case ${SEQTYPE} in
  dna) download_dna ${acc} ${gb} ;;
  gene) download_gene_info ${acc} ;;
  assembly) download_assembly ${acc} ${gff} ;;
  *) echo "Error: cannot retrieve data" && exit 1 ;;
  esac
}

declare -A params

# Set GFF to false by default
params[GFF]=false

# Set GB to false by default
params[GB]=false

optspec="hgb"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  g)
    params[GFF]=true
    ;;
  b)
    params[GB]=true
    ;;
  *)
    echo "Error: invalid option passed"
    usage
    ;;
  esac
done

main() {
  local acc=$1
  validate_accession ${acc}
  if [[ "${params[GFF]}" == true && "${params[GB]}" == true ]]; then
    download_sequence ${ACC} true true
  elif [[ "${params[GFF]}" == true && "${params[GB]}" == false ]]; then
    download_sequence ${ACC} true false
  elif [[ "${params[GFF]}" == false && "${params[GB]}" == true ]]; then
    download_sequence ${ACC} false true
  elif [[ "${params[GFF]}" == false && "${params[GB]}" == false ]]; then
    download_sequence ${ACC} false false
  else
    exit 10
  fi
}

main $1
