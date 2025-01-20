#!/usr/bin/env bash

usage() {
  echo "Perform sequence alignment with MUSCLE"
  echo
  echo "Usage:"
  echo "  run_mafft.sh [-h] [options] <file>|<dir>"
  echo
  echo "  <file>  a FASTA file containg the sequences to align"
  echo "  <dir>   a directory containing FASTA files for alignment"
  echo
  echo "Options:"
  echo "  -t    number of threads"
  echo "  -p    number of permutations"
  echo "  -s    set inital seed"
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
params[perm]=none
params[seed]=0
params[threads]=4

# Parse command line arguments
optspec="ho:p:s:t:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  p)
    params[perm]="${OPTARG}"
    ;;
  s)
    params[seed]="${OPTARG}"
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

  muscle_cmd="muscle "

  if [[ -f "${input}" ]]; then
    validate_fasta_file ${input}
    muscle_cmd+="-super5 ${FA}"
  elif [[ -d "${input}" ]]; then
    concatenate_fasta_files ${input}
    muscle_cmd+="-super5 /tmp/bioflow.seqs.fa"
  fi

  muscle_cmd+=" -output /tmp/bioflow.muscle.aln"

  echo "Aligning sequences with MUSCLE" 1>&2
  eval ${muscle_cmd}
  cat /tmp/bioflow.muscle.aln
  echo "Done!" 1>&2
}

main $1
