#!/usr/bin/env bash

PROJECT_ROOT=$(dirname $(realpath "$0"))

# ln -s -f ${PROJECT_ROOT}/src/fetch.mk fetch

LOCAL_DIR=~/.local/share/bioflow

symlink_src() {
  mkdir -p ${LOCAL_DIR}/bin
  for mk in $(find ${PROJECT_ROOT}/src -name "*.mk"); do
    ln -s -f ${mk} ${LOCAL_DIR}/bin
  done
}

create_aliases() {
  local target=${LOCAL_DIR}/aliases.sh
  rm ${target}
  touch ${target}
  echo "alias bf-align=\"make -f ${LOCAL_DIR}/bin/alignment.mk\"" >${target}
  echo "alias bf-phylo=\"make -f ${LOCAL_DIR}/bin/phylo.mk\"" >>${target}
  echo "alias bf-vc=\"make -f ${LOCAL_DIR}/bin/variant_calling.mk\"" >>${target}
  echo "alias bf-map=\"make -f ${LOCAL_DIR}/bin/mapping.mk\"" >>${target}
  echo "alias bf-qc=\"make -f ${LOCAL_DIR}/bin/qc.mk\"" >>${target}
  echo "alias bf-assemble=\"make -f ${LOCAL_DIR}/bin/assembly.mk\"" >>${target}
  echo "alias bf-fetch=\"make -f ${LOCAL_DIR}/bin/fetch.mk\"" >>${target}
  echo "alias bf-annotate=\"make -f ${LOCAL_DIR}/bin/annotation.mk\"" >>${target}

}

run_setup() {
  symlink_src
  create_aliases
}

run_setup
