# Scripts root

This is where all scripts for the pipeline live. It includes python and
rmarkdown, and hopefully will never include bash ;)

# Development

This directory is segregated by language and further by environment.

All common code is located in modules under `*/common`. Everything not under
`common` corresponds to script that run under one environment corresponding to
the `env.yml` file (which links back to `../envs` so they can be called by
snakemake rules easily).

Each of these environments includes `snakefmt` and `black` which will allow
editing/linting the snakemake files themselves. They also include
`mypy`/`flake8` (python) or `r-styler`/`r-lintr` (rmarkdown) depending on which
files are meant to be run in the environment. In the case of python, `mypy`
should be satisfied as all runtime dependencies are already included in each
environment.
