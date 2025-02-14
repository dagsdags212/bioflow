# Welcome to bioflow!

The `bioflow` project is a collection of Makefiles for orchestrating tasks used in genomic analysis. 
It provides a simple interface for running command-line programs in bioinformatics.
Each tool if configured with a reasonable set of defaults, removing the overhead of manual parameterization for budding bioinformaticians.
However, the user if free to edit the configuration for each tool, adding, removing, or replacing parameters as one sees fit for the analysis.

## The `BIOFLOW` variable

The full install path of bioflow is given by the `BIOFLOW` global variable. 

## Modules

`bioflow` is organized into modules, with each module containing recipes for running tools commonly used together in a workflow. 
This provides the user with a set of modular components for running individual parts of their pipeline. 
Use the `BIOFLOW` variable to easily access a recipe. For a list of parameters and commands, invoke the `help` command:

```bash
make -f $BIOFLOW/src/<module>/<recipe> help
```

As an example, we can run the `bwa` recipe (from the **mapping** module) for mapping FASTQ reads against a reference genome as follows:
```bash
make -f $BIOFLOW/src/mapping/bwa.mk REF=genome.fa R1=reads_1.fq R2=reads_2.fq run
```

However, we recommend that you use the provided alias for each recipe. All aliases follow the pattern: `bf-<recipe>`.
Using an alias, we can run the `bwa` recipe as follows:
```bash
bf-bwa REF=genome.fa R1=reads_1.fq R2=reads_2.fq run
```

A complete list of modules can be seen in @bioflow-modules.

:::{table} bioflow modules.
:label: bioflow-modules
:align: center

| Module | Functionality |
| ---- | ------------ |
| `alignment` | pairwise alignment, multiple sequence alignment |
| `annotation` | gene prediction, assembly annotation |
| `assembly` | sequence assembly, assembly evaluation |
| `fetch` | download biological data |
| `mapping` | read mapping, mapping assessment |
| `phylo` | maximum likelihood inference, bootstrapping phylogenies |
| `qc` | produce FASTQ summary reports, adapter trimming, quality filtering |
| `variant_calling` | detecting and filtering variants |

:::

## Environments

For each module, `bioflow` maintains a separate environment to avoid dependency incompatibility.
It uses the `conda`-family of environment managers (i.e., `conda`, `mamba`, `micromamba`) to download programs from a wide variety of package repositories.

The environment for a specific recipe can be intialized by running the `install` command:
```bash
bf-bwa install
```

Refer to the section on [environments](./environments.md) for a complete list of conda environment files.
