# Installation

Create a local copy of the github repository:
```bash
git clone https://github.com/dagsdags212/bioworkflows.git bioworkflows
cd bioworkflows
```

## Managing Environments

The project supports the use of `conda`, `mamba`, and `micromamba` for creating isolated environment to manage workflow dependencies. Choose a workflow and install its dependencies by invoking the `init` subcommand.

For example, we could create a micromamba environment for the **fetch.mk** workflow by running:
```bash
make -f src/fetch.mk init ENV_MANAGER=micromamba
```

Make sure to indicate the `ENV_MANAGER` parameter when creating a new environment, otherwise installation will fail.
