#!/usr/bin/env bash

usage() {
  echo "Retrieve sequencing reads from the SRA"
  echo
  echo "Usage: fetch_reads.sh [-h] [-v] [-a ACCESSION] [-x SPOTS]"
  echo "  -h  display this help message"
  echo "  -a  sequence project accession to fetch"
  echo "  -x  number of spots to download"
  echo "  -v  print tool version"
  exit
}

validate_read_accession() {
  local ACC=$1
  local SRA_REGEX=^[SED]R[PSXR].*

  if [[ "${ACC}" =~ ${SRA_REGEX} ]]; then
    declare -g ACC=${ACC}
  else
    echo "Error: invalid read accession"
    exit 1
  fi
}

download_reads_se() {
  local ACC=$1
  local X=$2
  local target=reads/${ACC}

  if [[ -n ${X} ]]; then
    echo "Downloading ${X} spots for ${ACC}"
    fastq-dump -X ${X} --origfmt -v -O ${target} ${ACC}
  else
    echo "Downloading ${ACC}"
    fastq-dump --origfmt -v -O ${target} ${ACC}
  fi
}

download_reads_pe() {
  local ACC=$1
  local X=$2
  local target=reads/${ACC}

  mkdir -p ${target}
  if [[ -n ${X} ]]; then
    echo "Downloading ${X} spots for ${ACC}"
    fastq-dump -X ${X} --split-3 --origfmt -v -O ${target} ${ACC}
  else
    echo "Downloading ${ACC}"
    fastq-dump --split-3 --origfmt -v -O ${target} ${ACC}
  fi
}

declare -A params

optspec="hpx:a:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  x)
    params[X]="$OPTARG"
    ;;
  a)
    params[ACC]="$OPTARG"
    ;;
  p)
    params[PE]=true
    ;;
  *)
    echo "Error: invalid option passed"
    usage
    ;;
  esac
done

main() {
  if [[ ! -n "${params[ACC]}" ]]; then
    echo "Error: accession not provided" && exit 1
  fi

  validate_read_accession ${params[ACC]}

  if [[ "${params[PE]}" == true ]]; then
    download_reads_pe ${ACC} ${params[X]}
  else
    download_reads_se ${ACC} ${params[X]}
  fi
}

main
