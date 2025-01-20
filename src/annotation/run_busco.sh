#!/usr/bin/env bash

usage() {
  echo "Perform sequence annotation with BUSCO"
  echo
  echo "Usage:"
  echo "  run_busco.sh [-h] [options] <file>"
  echo "  <file>  a FASTA file to annotate"
  echo
  echo "Options:"
  echo "  -t    number of threads"
  echo "  -m    BUSCO mode {genome,transcriptome,proteins}"
  echo "  -d    specify domain of sequence source {prokaryote,eukaryote}"
  echo "  -o    path to output directory"
  echo "  -h    display this help message"
  exit
}

declare -A params

# Specify default parameters
params[mode]="genome"
params[domain]="prokaryote"
params[outdir]="output/busco"
params[threads]=4

# Parse command line arguments
optspec="hm:d:o:t:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  m)
    params[mode]="${OPTARG}"
    ;;
  d)
    params[domain]="${OPTARG}"
    ;;
  o)
    params[outdir]="${OPTARG}"
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
  mkdir -p ${outdir}

  busco_cmd="busco -f -c ${params[threads]} "
  # Specify input file
  busco_cmd+="-i ${input} "
  # Specify BUCSO mode
  busco_cmd+="-m ${params[mode]} "
  # Indicate sequence domain
  if [[ "${params[domain]}" =~ ^[Pp]rok\(aryote\)$ ]]; then
    busco_cmd+="--auto-lineage-prok "
  elif [[ "${params[domain]}" =~ ^[Ee]uk\(aryote\)$ ]]; then
    busco_cmd+="--auto-lineage-euk "
  else
    busco_cmd+="--auto-lineage "
  fi
  # Specify output directory
  busco_cmd+="-o ${outdir}"

  echo "Performing genome annotation with BUSCO" 1>&2
  eval ${busco_cmd} && echo "Done!" 1>&2
}

main $1
