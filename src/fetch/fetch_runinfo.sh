#!/usr/bin/env bash

usage() {
  echo "Fetch metadata from a project accession"
  echo
  echo "Usage:"
  echo "  fetch_runinfo.sh [-h] <acc>"
  echo ""
  echo "  acc  sequencing project accession"
  echo
  echo "Options:"
  echo "  -h         display this help message"
  exit
}

mgnify2sra() {
  # Convert MGnify accession to SRA accession
  local id=$1
  local root=https://www.ebi.ac.uk/ena/browser/api/summary
  local url=${root}/${id}

  curl -X 'GET' ${url} -H 'accept: application/json' --silent |
    jq -r ".summaries[].project"
}

validate_project_id() {
  local mgnify_regex='^MGYS[0-9]{8}$'
  local sra_regex='^PRJ(NA|EB)[0-9]{4}$'

  if [[ "$1" =~ ${mgnify_regex} ]]; then
    declare -g ACC=$(mgnify2sra $1)
  elif [[ "$1" =~ ${sra_regex} ]]; then
    declare -g ACC=$1
  else
    echo "Error: invalid project accession"
    exit 1
  fi
}

optspec="h"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  *)
    echo "Error: invalid option passed"
    usage
    ;;
  esac
done

main() {
  validate_project_id $1
  echo "Fetching metadata for reads belonging to ${ACC}" 1>&2
  esearch -db sra -query ${ACC} | efetch -format runinfo
}

main $1
