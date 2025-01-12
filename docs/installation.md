# Installation

Create a local copy of the github repository:
```bash
git clone https://github.com/dagsdags212/bioflow.git
cd bioworkflows
```

:::{figure} ./assets/recordings/install.gif
Install bioworkflows using `git`.
:::

## Managing Environments

The project supports the use of [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html), [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html), and [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for creating isolated environment to manage workflow dependencies. Choose a workflow and install its dependencies by invoking the `init` subcommand.

For example, we could create a micromamba environment for the **fetch.mk** workflow by running:
```bash
make -f src/fetch.mk init ENV_MANAGER=micromamba
```

:::{figure} ./assets/recordings/install_deps.gif
Install workflow dependencies.
:::

Make sure to indicate the `ENV_MANAGER` parameter when creating a new environment, otherwise installation will fail. By default, the manager is set to `micromamba`.
