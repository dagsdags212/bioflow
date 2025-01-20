#!/usr/bin/env bash

usage() {
  echo "Perform sequence alignment with MAFFT"
  echo
  echo "Usage:"
  echo "  run_mafft.sh [-h] [options] <file>|<dir>"
  echo
  echo "  <file>  a FASTA file containg the sequences to align"
  echo "  <dir>   a directory containing FASTA files for alignment"
  echo
  echo "Options:"
  echo "  -t    number of threads"
  echo "  -o    gap opening penalty (default: 1.53)"
  echo "  -e    gap extension penality (default: 0.0)"
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
params[op]=1.53
params[ep]=0
params[maxiterate]=0
params[threads]=4

# Parse command line arguments
optspec="ho:e:t:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  o)
    params[op]=${OPTARG}
    ;;
  e)
    params[ep]=${OPTARG}
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

  mafft_cmd="mafft --maxiterate ${params[maxiterate]} --op ${params[op]} --ep ${params[ep]} --thread ${params[threads]} "

  if [[ -f "${input}" ]]; then
    validate_fasta_file ${input}
    mafft_cmd+="${FA}"
  elif [[ -d "${input}" ]]; then
    concatenate_fasta_files ${input}
    mafft_cmd+="/tmp/bioflow.seqs.fa"
  fi

  echo "Aligning sequences with MAFFT" 1>&2
  eval ${mafft_cmd}
  echo "Done!" 1>&2
}

main $1
