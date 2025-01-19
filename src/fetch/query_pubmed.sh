#!/usr/bin/env bash

main() {
  local query=$1
  esearch -db pubmed -query "${query}" |
    efetch -format docsum |
    # Extract PMID, publication year, and title of article
    xtract -pattern DocumentSummary -tab "|" -element Id -year PubDate -element Title |
    sort -t "|" -k2,2nr -k3,3 |
    column -t -s "|" -N PMID,YEAR,TITLE -W TITLE
}

main $1
