#!/usr/bin/env bash

# Path pointing to the bioflow repo.
BIOFLOW=~/.local/share/bioflow

# Set config directory.
[ -d ~/.config ] && CONFIG=~/.config/bioflow || CONFIG=~/.bioflow

# Empty target directory.
clean_target_path() {
  rm -rf ${BIOFLOW}
}

# Copy repository into target directory.
create_local_copy() {
  echo "Copying project files to ${BIOFLOW}"
  mkdir -p ${BIOFLOW}
  cp -r ./* ${BIOFLOW}
}

create_config_dir() {
  echo "Creating config directory"
  mkdir -p ${CONFIG}
}

# Export BIOFLOW variable.
export_bioflow_path() {
  local config_path=~/.zshrc
  local export_line="export BIOFLOW=${BIOFLOW}"
  echo "Exporting BIOFLOW path"
  grep "${export_line}" ${config_path} || echo ${export_line} >>${config_path}
}

set_aliases() {
  local alias_file=${CONFIG}/aliases.sh

  echo "Generating aliases"
  # Remove if already exists, and create a new file
  [ -f "${alias_file}" ] && rm -f ${alias_file} || touch ${alias_file}

  # Fetch module aliases.
  echo "# fetch module" >>${alias_file}
  echo 'alias bf-bioproject="make -f ${BIOFLOW}/src/fetch/bioproject.mk"' >>${alias_file}
  echo 'alias bf-sra="make -f ${BIOFLOW}/src/fetch/sra.mk"' >>${alias_file}
  echo 'alias bf-genbank="make -f ${BIOFLOW}/src/fetch/genbank.mk"' >>${alias_file}
  echo 'alias bf-pubmed="make -f ${BIOFLOW}/src/fetch/pubmed.mk"' >>${alias_file}
  echo 'alias bf-pdb="make -f ${BIOFLOW}/src/fetch/pdb.mk"' >>${alias_file}

  # QC module aliases.
  echo 'alias bf-fastqc="make -f ${BIOFLOW}/src/qc/fastqc.mk"' >>${alias_file}
  echo 'alias bf-multiqc="make -f ${BIOFLOW}/src/qc/multiqc.mk"' >>${alias_file}
  echo 'alias bf-fastp="make -f ${BIOFLOW}/src/qc/fastp.mk"' >>${alias_file}

  # Source alias file
  local config_path=~/.zshrc
  local export_line="source ${alias_file}"
  grep "${export_line}" ${config_path} || echo ${export_line} >>${config_path}
}

run_setup() {
  clean_target_path
  create_config_dir
  create_local_copy
  export_bioflow_path
  set_aliases
  sleep 1
  echo "Bioflow has been installed"
}

# Run complete setup
run_setup
