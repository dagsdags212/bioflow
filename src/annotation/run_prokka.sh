#!/usr/bin/env bash

usage() {
  echo "Annotate prokaroyotic sequences with Prokka"
  echo
  echo "Usage:"
  echo "  run_prokka.sh [-h] [options] <file>"
  echo
  echo "  <file>  a FASTA file to annotate"
  echo
  echo "Options:"
  echo "  -t    number of threads"
  echo "  -o    specify output directory"
  echo "  -p    filename output prefix"
  echo "  -h    display this help message"
  exit
}

declare -A params

# Specify default parameters
params[outdir]=output/prokka
params[prefix]=""
params[threads]=4

# Parse command line arguments
optspec="ht:o:p:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  o)
    params[outdir]="${OPTARG}"
    ;;
  p)
    params[prefix]="${OPTARG}"
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
  local outdir=${params[outdir]}
  echo "Creating output directory: ${outdir}" 1>&2

  prokka_cmd="prokka --gffver 3 --force "
  # Specify output directory
  prokka_cmd+="--outdir ${outdir} "
  # Specify input FASTA file
  prokka_cmd+="${input}"

  echo "Annotating prokaroyotic sequence with Prokka" 1>&2
  eval ${prokka_cmd} && echo "Done!" 1>&2
}

main $1
