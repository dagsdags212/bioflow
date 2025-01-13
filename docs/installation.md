# Getting Started

## Installation

Create a local copy of the github repository:
```bash
git clone https://github.com/dagsdags212/bioflow.git
cd bioworkflows
```

:::{figure} ./assets/recordings/install.gif
Install bioworkflows using `git`.
:::

## Installing GNU Make

Prior to running a workflow, ensure that [GNU Make](https://www.gnu.org/software/make/) is installed in your system. Check if the program is installed:
```bash
make -v
```

If not, the software is available thru APT (`apt-get`) in UNIX-based systems, Homebrew (`brew`) in MacOS, and Chocolatey (`choco`) in Windows.

To learn more about GNU Make, read the official manual [here](https://www.gnu.org/software/make/manual/make.html).

## Running a Workflow

Once `make` is installed in your system, run a supported pipeline by invoking the following command:
```bash
make -f src/<workflow> [arguments]
```

`<workflow>` is the Makefile which stores the workflow and `[arguments]` is a set of optional/required parameters used for a specific workflow. All workflows are stored under the `src` directory. 

To get a complete list of arguments for a workflow, invoke the following command:
```bash
make -f src/<workflow> help
```

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
