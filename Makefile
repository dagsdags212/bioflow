# Absoluate path for the `bioflow` directory
PROJECT_ROOT := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))

.PHONY: serve config format

all:
	@echo "Installing bioflow..."
	. $(PROJECT_ROOT)/install.sh

serve:
	myst start

config:
	nvim ./myst.yml

format:
	ruff format scripts/
