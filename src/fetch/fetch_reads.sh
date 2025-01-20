#!/usr/bin/env bash

usage() {
  echo "Retrieve sequencing reads from the SRA"
  echo
  echo "Usage:"
  echo "fetch_reads.sh [-h] [-x SPOTS] <acc>"
  echo
  echo "  <acc>  SRA read accession"
  echo
  echo "Options:"
  echo "  -h  display this help message"
  echo "  -x  number of spots to download"
  echo "  -v  print tool version"
  exit
}

validate_read_accession() {
  local acc=$1
  local sra_regex="^[SED]R[PSXR].*"

  if [[ "${acc}" =~ ${sra_regex} ]]; then
    declare -g ACC=${acc}
  else
    echo "Error: invalid read accession"
    exit 1
  fi
}

download_reads_se() {
  local acc=$1
  local x=$2
  local target=reads/${acc}

  if [[ "${x}" != "" ]]; then
    echo "Downloading ${x} spots for ${acc}"
    fastq-dump -X ${x} --origfmt -v -O ${target} ${acc}
  else
    echo "Downloading ${acc}"
    fastq-dump --origfmt -v -O ${target} ${acc}
  fi
}

download_reads_pe() {
  local acc=$1
  local x=$2
  local target=reads/${acc}

  mkdir -p ${target}
  if [[ "${x}" != "" ]]; then
    echo "Downloading ${x} spots for ${acc}"
    fastq-dump -X ${x} --split-3 --origfmt -v -O ${target} ${acc}
  else
    echo "Downloading ${acc}"
    fastq-dump --split-3 --origfmt -v -O ${target} ${acc}
  fi
}

declare -A params

optspec="hpx:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  x)
    params[X]="$OPTARG"
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
  local acc=$1
  if [[ ! -n "${acc}" ]]; then
    echo "Error: accession not provided" && exit 1
  fi

  validate_read_accession ${acc}

  if [[ "${params[PE]}" == true ]]; then
    if [[ "${params[X]}" != "" ]]; then
      download_reads_pe ${ACC} ${params[X]}
    else
      download_reads_pe ${ACC}
    fi
  else
    if [[ "${params[X]}" != "" ]]; then
      download_reads_se ${ACC} ${params[X]}
    else
      download_reads_se ${ACC}
    fi
  fi
}

main $1
