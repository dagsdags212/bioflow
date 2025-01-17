#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $(realpath "$0"))

# Path where bioflow will live in the local system
LOCAL_DIR=~/.local/share/bioflow

RC_FILE=~/.zshrc

# Create a symbolic link to each makefile workflow
symlink_src() {
  local target=${LOCAL_DIR}/bin
  mkdir -p ${target}
  for mk in $(find ${PROJECT_ROOT}/src -name "*.mk"); do
    ln -s -f ${mk} ${target}
  done
}

create_aliases() {
  local target=${LOCAL_DIR}/aliases.sh
  rm ${target}
  touch ${target}
  # Treat each makefile as a binary executable
  echo "alias bf-align=\"make -f ${LOCAL_DIR}/bin/alignment.mk\"" >${target}
  echo "alias bf-phylo=\"make -f ${LOCAL_DIR}/bin/phylo.mk\"" >>${target}
  echo "alias bf-vc=\"make -f ${LOCAL_DIR}/bin/variant_calling.mk\"" >>${target}
  echo "alias bf-map=\"make -f ${LOCAL_DIR}/bin/mapping.mk\"" >>${target}
  echo "alias bf-qc=\"make -f ${LOCAL_DIR}/bin/qc.mk\"" >>${target}
  echo "alias bf-assemble=\"make -f ${LOCAL_DIR}/bin/assembly.mk\"" >>${target}
  echo "alias bf-fetch=\"make -f ${LOCAL_DIR}/bin/fetch.mk\"" >>${target}
  echo "alias bf-annotate=\"make -f ${LOCAL_DIR}/bin/annotation.mk\"" >>${target}
}

append_to_rc() {
  grep -qxF "source \"/home/dagsdags/.local/share/bioflow/aliases.sh\"" ~/.zshrc || echo "source \"${LOCAL_DIR}/aliases.sh\"" >>${RC_FILE}
}

run_setup() {
  symlink_src
  create_aliases
  append_to_rc
}

# Run complete setup
run_setup
