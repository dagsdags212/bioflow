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

run_setup() {
  clean_target_path
  create_local_copy
  sleep 1
  echo "Done!"
}

# Run complete setup
run_setup
