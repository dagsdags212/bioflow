#!/usr/bin/env bash

QUERY=$1

echoerr() { echo "$@" 1>&2; }

query_pubmed() {
  # display message prompt
  # echoerr "Querying PubMed with search term: '${QUERY}'"

  esearch -db pubmed -query "${QUERY}" |
    efetch -format docsum |
    # Extract PMID, publication year, and title of article
    xtract -pattern DocumentSummary -tab "|" -element Id -year PubDate -element Title |
    sort -t "|" -k2,2nr -k3,3 |
    column -t -s "|" -N PMID,YEAR,TITLE -W TITLE
}

query_pubmed
