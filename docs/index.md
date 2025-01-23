# Welcome to bioflow!

The `bioflow` project is a collection of Makefiles for orchestrating tasks used in genomic analysis. It provides a simple
interface for running command-line programs in bioinformatics. Each tool if configured with a reasonble set of defaults, 
removing the overhead of manual paramterization for budding bioinformaticians. However, the user if free to edit the configuration for each tool, adding, removing, or replacing parameters as one sees fit for the analysis.


## Modules

`bioflow` is organized into modules, with each module containing recipes for running tools commonly used together in a workflow. This provides the user with a set of modular components for running individual parts of their pipeline. Recipes (or submodules) can be invoked simply by running `make -f <recipe>`, where `<recipe>` is the path to a submodule. For a list of parameters, invoke the help command for a specific recipe:

```bash
make -f <recipe> help
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

For each module, `bioflow` maintains a separate environment to avoid dependency incompatibility. It uses the `conda`-family of environment managers (i.e., `conda`, `mamba`, `micromamba`) to download programs from a wide variety of package repositories.

All environments can be locally installed by downloading them from a YAML config. Refer to the section on [environments](./environments.md) for a complete list of conda environment files.
