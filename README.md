# giab-ai-ebm

## deployment

Install the environment to run in snakemake by running this command at the root
of this repo:

```
mamba env create -f env.yml
```

For development packages (`black`, `flake8`, etc) run this (must be done after
creating the environment:

```
mamba install --file dev.txt -n snakemake-ebm -c conda-forge -c bioconda
```

## running the pipeline

This pipeline uses `snakemake` to process the data and build the models and
`dvc` to store the results in a hopefully-sane way.

`dvc` is implemented using a single stage pipeline which called the `snakemake`
pipeline. It's sole input is the `config/dynamic.yaml` file (which fully
describes an experiment) and its outputs are the contents of
`results/annotated_inputs` and `results/ebm` (which contain the input dataframes
and model results respectively).

Results can then be retrieved for a specific experiment by referencing the tag
(or commit if not available) of the version of this git repo used to generate
it.

### auto

Run the entire pipeline and store results in the cloud:

```
dvc repro
dvc push
```

### manual

Run the entire pipeline via snakemake with the following:

```
snakemake -p --verbose -r -j 1 --use-conda --rerun-incomplete
```

Store results in the cloud:

```
dvc commit
dvc push
```

### slurm

Submit a job to run the pipeline via slurm using this [script](snakemake-slurm):

```
sbatch snakemake-slurm
```

This does not invoke `dvc` so `dvc commit/push` will need to be run manually
to put the results in the cloud.

## retrieving results

Assuming that `dvc commit/push` was properly invoked on the results in question,
data can be accessed for any given commit/tag. To list the files for a given
tag:

```
dvc list --rev <tag> https://gitlab.nist.gov/gitlab/njd2/giab-ai-ebm.git --dvc-only -R
```

To pull the data from the `results/ebm` folder (eg the model output) from a
given tag/commit to a local file path:

```
dvc get --rev <tag> https://gitlab.nist.gov/gitlab/njd2/giab-ai-ebm.git results/ebm -o <local_path>
```

To track in a local `dvc` repo (eg for rigorous analysis) use `import` instead
of `get` in the above command (must be done in a `dvc` repo, which can be
initialized with `git init && dvc init`)

## configuring the pipeline

The configuration is split into [static](config/static.yml) and
[dynamic](config/dynamic.yml) components. The former is for downloaded resources
and parameters used to generate the annotated dataframes. The latter is for
selecting various features within the annotated dataframes, applying
transformations, and selecting hyperparameters when training the EBM models.

## pipeline output

Each entry under `ebm_runs` in the dynamic config corresponds to one EBM run
with its corresponding features and settings to be used. After running the
pipeline, each run should have a directory under `results/ebm` names like
`<git_tag>_<run_entry_name>` where `<git_tag>` is the current tag of the repo
(or the commit if there is none) and `<run_entry_name>` is the key under
`ebm_runs` in the [dynamic config](config/dynamic.yml).

Each directory will contain the input tsv of data used to train the EBM, a
config yml file with all settings used to train the EBM, and python pickles for
the X/Y train/test datasets as well as a pickle for the final model itself.
