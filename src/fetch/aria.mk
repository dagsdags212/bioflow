#
# Download data using the aria2c protocol.
#

# Makefile preamble.
SHELL := bash
.DELETE_ON_ERROR:
.ONESHELL:
.SHELLFLAGS := -eu -o pipefail -c
MAKEFLAGS += --warn-undefined-variables --no-print-directory

# Micromamba environment.
ENV = bf-fetch-aria

# Run command within environment.
ENV_RUN = micromamba run -n $(ENV)

# The target URL to be downloaded.
URL ?= https://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.fa.gz

# The directory to store the data.
OUT ?= data

# The file destination.
FILE ?= $(DIR)/$(notdir $(URL))

# Aria2c flags.
FLAGS := -x 5 -c --summary-interval=10

# Print usage.
help::
	@echo ""
	@echo "aria.mk: download data from a link"
	@echo ""
	@echo "Usage:"
	@echo "  bf-aria URL=<URL> [options] run"
	@echo ""
	@echo "Options:"
	@echo "  URL        an FTP/URL/torrent link"
	@echo "  OUT        a directory path for storing downloads"
	@echo ""

# Download the file with aria2c.
${FILE}:
	# Create output directory.
	mkdir -p $(dir ${FILE})

	# Download the file
	$(ENV_RUN) aria2c ${FLAGS} -o ${FILE} ${URL}

run: ${FILE}
	@ls -lh $<

# Delete the file.
clean:
	rm -f $(FILE)

# Alternative rule for clean.
run!: clean

# Run test suite.
test: clean run

DEPS := aria2
# Display dependencies.
install::
	micromamba create -n ${ENV}
	${ENV_RUN} micromamba install ${DEPS}

.PHONY: help run run! test install

