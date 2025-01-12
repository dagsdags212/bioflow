#!/usr/bin/env bash

ROOT=https://www.ebi.ac.uk/ena/browser/api/summary
ID=$1

URL=${ROOT}/${ID}

curl -X 'GET' ${URL} -H 'accept: application/json' --silent |
  jq -r ".summaries[].project"
