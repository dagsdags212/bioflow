# Getting Started

## Installation

Create a local copy of the github repository:
```bash
git clone https://github.com/dagsdags212/bioflow.git
cd bioflow
```

A bash script is provided to install the entry points of each workflow to your system. Inside the `bioflow` directory, run the `install.sh` script in the context of the current shell:
```bash
./install.sh
```

:::{figure} ./assets/recordings/install.gif
Install bioworkflows using `git`.
:::

After restarting the shell, the `BIOFLOW` variable which points to the project root of bioflow will be globally available. This can be used by the user for invoking the recipes in each module.

## Installing GNU Make

Prior to running a workflow, ensure that [GNU Make](https://www.gnu.org/software/make/) is installed in your system. Check if the program is installed:
```bash
make -v
```

If not, the software is available thru APT (`apt-get`) in UNIX-based systems, Homebrew (`brew`) in MacOS, and Chocolatey (`choco`) in Windows.

To learn more about GNU Make, read the official manual [here](https://www.gnu.org/software/make/manual/make.html).

## Running a Workflow

Upon running the install script, all makefiles are accessible within the `$HOME/.local/share/bioflow/bin/` directory. We could run them as a normal makefile by invoking `make`, specifying the path of the workflow using the `-f` flag:
```bash
make -f <workflow> [arguments]
```

`<workflow>` is the Makefile which stores the workflow and `[arguments]` is a set of optional/required parameters used for a specific workflow. All workflows are stored under the `src` directory. However, this is **not** the recommended approach.

Each workflow is provided with an entry point that can be treated as an executable. All commands are prefixed with `bf-` and a complete list is provided below:

- `bf-align` for sequence alignment
- `bf-annotate` for gene prediction and genome annotation
- `bf-assembly` for sequence assembly and draft genome assessment
- `bf-fetch` for retrieving biological data
- `bf-map` for read mapping
- `bf-phylo` for tree inference
- `bf-qc` for quality control, adapter trimming, and quality filtering
- `bf-vc` for variant calling

As an example, we could download all sequencing reads associated with a PRJNA accession by running the following command:
```bash
bf-fetch sra PRJNA=PRJNA1066786
```

To get a complete list of arguments for the `bf-fetch` workflow, invoke the following command:
```bash
bf-fetch help
```

## Managing Environments

The project supports the use of [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html), [mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html), and [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html) for creating isolated environment to manage workflow dependencies. Choose a workflow and install its dependencies by invoking the `init` subcommand.

For example, we could create a micromamba environment for the **fetch.mk** workflow by running:
```bash
bf-fetch init ENV_MANAGER=micromamba
```

:::{figure} ./assets/recordings/install_deps.gif
Install workflow dependencies.
:::

Make sure to indicate the `ENV_MANAGER` parameter when creating a new environment, otherwise installation will fail. By default, the manager is set to `micromamba`.
