#!/usr/bin/env bash

usage() {
  echo "Perform de novo gene prediction with Prodigal"
  echo
  echo "Usage:"
  echo "  run_prodigal.sh [-h] [options] <file>"
  echo
  echo "  <file>  a FASTA file to annotate"
  echo
  echo "Options:"
  echo "  -t    number of threads"
  echo "  -o    filename for output predictions (file extension will be inferred)"
  echo "  -f    output format {gbk,gff,sco}"
  echo "  -h    display this help message"
  exit
}

declare -A params

# Specify default parameters
params[outname]="predicted"
params[outfmt]="gff"
params[threads]=4

# Parse command line arguments
optspec="hf:t:"
while getopts "${optspec}" optchar; do
  case "${optchar}" in
  h)
    usage
    ;;
  f)
    params[outfmt]="${OPTARG}"
    ;;
  o)
    params[outname]="${OPTARG}"
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
  local outdir=output/prodigal
  mkdir -p ${outdir}

  prodigal_cmd="prodigal "
  # Specify input file
  prodigal_cmd+="-i ${input} "
  # Write nucleotide sequence of genes
  prodigal_cmd+="-d ${outdir}/${params[outname]}.fna "
  # Write amino acid translation of genes
  prodigal_cmd+="-a ${outdir}/${params[outname]}.faa "
  # Specify format of prodigal output
  prodigal_cmd+="-f ${params[outfmt]} "
  # Write prodigal output in a specified format (default: gff)
  prodigal_cmd+="-o ${outdir}/${params[outname]}.${params[outfmt]}"

  echo "Performing ab initio gene prediction using Prodigal" 1>&2
  eval ${prodigal_cmd} && echo "Done!" 1>&2
}

main $1
