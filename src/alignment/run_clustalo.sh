#!/usr/bin/env bash

usage() {
  echo "Perform sequence alignment with ClustalO"
  echo
  echo "Usage:"
  echo "  clustalo.sh [-h] [options] <file>|<dir>"
  echo
  echo "  <file>  a FASTA file containg the sequences to align"
  echo "  <dir>   a directory containing FASTA files for alignment"
  echo
  echo "Options:"
  echo "  -t    number of threads"
  echo "  -p    number of iterations"
  echo "  -f    output format {fasta,clustal,msf,phylip,selex,stockholm,vienna}"
  echo "  -h    display this help message"
  exit
}

validate_fasta_file() {
  local fasta=$1
  # Check if input is a file
  if [[ ! -f "${fasta}" ]]; then
    echo "Error: input file not found"
    exit 1
  fi
  # Check for header line
  header=$(head -n 1 ${fasta})
  if [[ "${header:0:1}" != ">" ]]; then
    echo "Error: header line not found"
    exit 1
  fi
  # All checks passed
  declare -g FA=${fasta}
}

concatenate_fasta_files() {
  local dirpath=$1
  # Check if path points to a directory
  if [[ ! -d "${dirpath}" ]]; then
    echo "Error: not a directory"
    exit 1
  fi
  # Validate each fasta file
  fasta_files=$(find ${dirpath} -type f -name "*.fn?a\(sta\)?$")
  if [[ $(echo ${fasta_files} | wc -l) -gt 1 ]]; then
    echo "Detected multiple FASTA files" 1>&2
  fi

  for fa in ${fasta_files[@]}; do
    echo "Validating FASTA file: ${fa}" 1>&2
    validate_fasta_file ${fa}
  done
  # Concatenate all fasta files within directory
  cat ${dirpath}/* >/tmp/bioflow.seqs.fa
}

declare -A params

# Specify default parameters
params[format]="fasta"
params[iter]=0
params[threads]=4

# Parse command line arguments
optspec="hf:p:t:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  f)
    params[format]=${OPTARG}
    ;;
  p)
    params[iter]=${OPTARG}
    ;;
  t)
    params[threads]="${OPTARG}"
    ;;
  *)
    echo "Error: invalid option passed"
    usage
    ;;
  esac
done

main() {
  local input=$1

  clustalo_cmd="clustalo --threads=${params[threads]} --outfmt=${params[format]} "

  if [[ -f "${input}" ]]; then
    validate_fasta_file ${input}
    clustalo_cmd+="-i ${FA}"
  elif [[ -d "${input}" ]]; then
    concatenate_fasta_files ${input}
    clustalo_cmd+="-i /tmp/bioflow.seqs.fa"
  fi

  echo "Aligning sequences with Clustal-Omega" 1>&2
  eval ${clustalo_cmd}
  echo "Done!" 1>&2
}

main $1
