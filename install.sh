#!/usr/bin/env bash

# Path pointing to the bioflow repo.
BIOFLOW=~/.local/share/bioflow

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

# Export BIOFLOW variable.
export_bioflow_path() {
  local config_path=~/.zshrc
  local export_line="export BIOFLOW=${BIOFLOW}"
  echo "Exporting BIOFLOW path"
  grep "${export_line}" ${config_path} || echo ${export_line} >>${config_path}
}

run_setup() {
  clean_target_path
  create_local_copy
  export_bioflow_path
  sleep 1
  echo "Done!"
}

# Run complete setup
run_setup
