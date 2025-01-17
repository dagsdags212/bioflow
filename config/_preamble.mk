SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
.ONESHELL:
MAKEFLAGS += --warn-undefined-variables --no-builtin-rules
