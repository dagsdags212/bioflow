#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $(realpath "$0"))

# Path where bioflow will live in the local system
LOCAL_DIR=~/.local/share/bioflow

RC_FILE=~/.zshrc

clean_target_path() {
  rm -rf ${LOCAL_DIR}
}

create_local_copy() {
  echo "Copying project files to ${LOCAL_DIR}"
  mkdir -p ${LOCAL_DIR}
  cp -r ./* ${LOCAL_DIR}
}

# Create a symbolic link to each makefile workflow
symlink_src() {
  local target=${LOCAL_DIR}/src
  mkdir -p ${target}
  for mk in $(find ${PROJECT_ROOT}/src -name "*.mk"); do
    ln -s -f ${mk} ${target}
  done
}

create_aliases() {
  echo "Setting up entry points to makefiles"
  local target=${LOCAL_DIR}/aliases.sh
  rm -f ${target}
  touch ${target}
  # Treat each makefile an executable
  echo "alias bf-align=\"make -f ${LOCAL_DIR}/src/alignment.mk\"" >${target}
  echo "alias bf-phylo=\"make -f ${LOCAL_DIR}/src/phylo.mk\"" >>${target}
  echo "alias bf-vc=\"make -f ${LOCAL_DIR}/src/variant_calling.mk\"" >>${target}
  echo "alias bf-map=\"make -f ${LOCAL_DIR}/src/mapping.mk\"" >>${target}
  echo "alias bf-qc=\"make -f ${LOCAL_DIR}/src/qc.mk\"" >>${target}
  echo "alias bf-assemble=\"make -f ${LOCAL_DIR}/src/assembly.mk\"" >>${target}
  echo "alias bf-fetch=\"make -f ${LOCAL_DIR}/src/fetch.mk\"" >>${target}
  echo "alias bf-annotate=\"make -f ${LOCAL_DIR}/src/annotation.mk\"" >>${target}
}

append_to_rc() {
  echo "Appending source file to RC config"
  grep -qxF "source \"/home/dagsdags/.local/share/bioflow/aliases.sh\"" ~/.zshrc || echo "source \"${LOCAL_DIR}/aliases.sh\"" >>${RC_FILE}
}

run_setup() {
  clean_target_path
  create_local_copy
  # symlink_src
  create_aliases
  append_to_rc

  sleep 1
  echo "Done!"
}

# Run complete setup
run_setup
