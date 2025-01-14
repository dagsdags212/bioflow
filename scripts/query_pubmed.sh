#!/usr/bin/env bash

QUERY=$1

query_pubmed() {
  # display message prompt
  echo
  echo "Querying PubMed with search term: '${QUERY}'"
  echo

  esearch -db pubmed -query "${QUERY}" |
    efetch -format docsum |
    # Extract PMID, publication year, and title of article
    xtract -pattern DocumentSummary -sep "\t" -element Id -year PubDate -element Title |
    align-columns |
    # Most recent articles come first
    sort-table -k 2 -r
}

query_pubmed
