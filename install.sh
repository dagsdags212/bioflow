#!/usr/bin/env bash

# Path where bioflow will live in the local system
BIOFLOW_PREFIX=~/.local/share/bioflow

# Determine default shell
determine_shell() {
  if [[ "${SHELL}" == "/usr/bin/zsh" ]]; then
    declare -g RC_FILE=~/.zshrc
  elif [[ "${SHELL}" == "/usr/bin/bash" ]]; then
    declare -g RC_FILE=~/.bashrc
  fi
}

clean_target_path() {
  rm -rf ${BIOFLOW_PREFIX}
}

create_local_copy() {
  echo "Copying project files to ${BIOFLOW_PREFIX}"
  mkdir -p ${BIOFLOW_PREFIX}
  cp -r ./* ${BIOFLOW_PREFIX}
}

create_setup_file() {
  local target=${BIOFLOW_PREFIX}/setup.sh
  echo "Cleaning target directory"
  rm -f ${target}
  touch ${target}

  echo "Creating setup file"
  # export bioflow path
  echo "export BIOFLOW_PREFIX=~/.local/share/bioflow" >${target}
  # Treat each makefile an executable
  echo "Generating entry point for each module"

  echo "alias bf-fetch=\"make -f ${BIOFLOW_PREFIX}/src/fetch/Makefile\"" >>${target}
  echo "alias bf-align=\"make -f ${BIOFLOW_PREFIX}/src/alignment/Makefile\"" >>${target}
  echo "alias bf-annotate=\"make -f ${BIOFLOW_PREFIX}/src/annotation/Makefile\"" >>${target}
  echo "alias bf-assemble=\"make -f ${BIOFLOW_PREFIX}/src/assembly/Makefile\"" >>${target}
  echo "alias bf-map=\"make -f ${BIOFLOW_PREFIX}/src/mapping/Makefile\"" >>${target}

  # TODO: migrate to new module structure
  echo "alias bf-phylo=\"make -f ${BIOFLOW_PREFIX}/src/phylo.mk\"" >>${target}
  echo "alias bf-vc=\"make -f ${BIOFLOW_PREFIX}/src/variant_calling.mk\"" >>${target}
  echo "alias bf-qc=\"make -f ${BIOFLOW_PREFIX}/src/qc.mk\"" >>${target}
}

append_to_rc() {
  echo "Appending source file to RC config"
  grep -qxF "source \"/home/dagsdags/.local/share/bioflow/setup.sh\"" ~/.zshrc || echo "source \"${BIOFLOW_PREFIX}/setup.sh\"" >>${RC_FILE}
}

run_setup() {
  determine_shell
  clean_target_path
  create_local_copy
  create_setup_file
  append_to_rc

  sleep 1
  echo "Done!"
}

# Run complete setup
run_setup
