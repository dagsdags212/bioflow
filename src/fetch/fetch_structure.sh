#!/usr/bin/env bash

usage() {
  echo "Retrieve a protein structure from PDB"
  echo
  echo "Usage: fetch_sequence.sh [-h] [-v] [ACCESSION]"
  echo
  echo "  accession  a four-character alphanumeric identifier"
  echo
  echo "  -h         display this help message"
  echo "  -v         print tool version"
  exit
}

validate_pdb_id() {
  local acc=$1
  local pdb_regex="^[[:alnum:]]{4}$"

  if [[ "${acc}" =~ ${pdb_regex} ]]; then
    declare -g ACC=${acc}
  else
    echo "Error: invalid PDB accession"
    exit 1
  fi
}

fetch_structure() {
  local acc=$1

  mkdir -p data/
  echo "Attempting to fetch structure file for ${acc}"
  pdb_fetch ${acc} >${acc}.pdb
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
  local acc=$1
  validate_pdb_id ${acc}
  fetch_structure ${ACC}
}

main $1
