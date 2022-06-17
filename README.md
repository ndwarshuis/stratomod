# giab-ai-ebm

The main pipeline to run EBM experiments using `snakemake` and `dvc`.

## Workflow

### Snakemake

Each experiment is 'built' using `snakemake` which runs all commands to retrieve
data, wrangle the input dataframes, and train the models.

This repo is designed such that (for the most part) the only flags required to
run `snakemake` are `--configfile` and `--profile`. This enables the experiement
configuration and runtime respectively to be tracked in git and easily
implemented in higher-level frameworks if needed.

### Profiles

In order to have snakemake run reproducibly on many different
machines/architectures, several profiles are provided at
`workflow/profiles/<profile_name`.

For now, there are two profiles:
- nisaba: for running on the NIST Nisaba cluster (specifically with slurm)
- local: for running on a local machine

### Configurations

Each experiment is configured using `config/dynamic.yml`. By convention, this
file is not tracked in any branch except for experiment branches (see below)
and needs to be created and tracked for each individual experiment.

The file at `config/dynamic-testing.yml` is a small-scale "experiment" used for
testing.

### DVC

While `snakemake` is used to run the actual experiments, `dvc` is deployed as a
wrapper around `snakemake` to store the results in a (hopefully) sane manner.
Note that `dvc` is only necessary/useful in the context of running experiments;
for testing and development it is easier to simply run `snakemake` directly.

The `dvc` 'pipeline' (see `dvc.yaml`) consistes of one stage whose sole purpose
is to run `snakemake`. Generally `dvc.yaml` doesn't need to be edited.

The behavior of `dvc` is controlled using `config/dvc-params.yml`. Here the
configuration file and the profile for `snakemake` (as described above) can be
set. This should be modified for each experiment.

Upon running `dvc repro` (see below) the `dvc.lock` file will link the
config/profile paramaters, the configuration file itself, and model results to a
specific git commit. These can then be pushed to an external data store (S3)
using `dvc push` and then imported somewhere else using `dvc import` and the
desired git commit which holds the `dvc.lock` data.

## Deployment

Install the environment to run in snakemake by running this command at the root
of this repo:

```
mamba env create -f env.yml
```

For development packages (`black`, `flake8`, etc) run this after creating the
environment:

```
mamba install --file dev.txt -n snakemake-ebm -c conda-forge -c bioconda
```

## Development and Experiment Workflow

Both development and experiments are tracked using git branches. There are two
main branches: `master` and `develop`

### Adding a New Feature

Make a new branch off of develop for the new feature, then merge into develop
when done (note `--no-ff`).

```
git checkout develop
git branch -n <new_feature>
git checkout <new_feature>

# do a bunch of stuff...

git checkout develop
git merge --no-ff <new_feature>
```

After feature(s) have been added and all tests have succeeded, update changelog,
add tag, and merge into master. Use semantic versioning for tags.

```
# update changelog
vim CHANGELOG.md

git commit

git tag vX.Y.Z
git checkout master
git merge --no-ff vX.Y.Z
```

NOTE: do not add an experiment-specific configuration to `master` or `develop`.
The yml files in `config` for these branches are used for testing. See below
for how to add an experiment.

### Adding an Experiment

Experiments are added by branching off master using a specific revision. By
convention, experiment branches should be prefixed with `x_`.

```
git checkout master
git branch x_examine_penguin_genes
git checkout x_examine_penguin_genes
```

The only modification to make on these branches is creating/editing the
`config/dynamic.yml` file which holds the configuration for the experiment. Once
this is modified to the desired state, commit it and run the pipeline as
described below.

If the master branch is updated to a new version, either merge that tag into the
experiment branch or create an entirely new experiment.

Don't merge any experiment branches back into master.

## Running the Pipeline

### Auto (dvc + snakemake)

Edit the dvc params with your favorite text editor to change the profile and
config as desired. If running on Nisaba, set the profile to "nisaba" (otherwise
local):

```
notepad config/dvc-params.yml
```

Run the entire pipeline and store results in the cloud:

```
dvc repro
dvc push
```

`dvc repro` will block the terminal so it is recommended to run it in the
background with `&` or use your favorite multiplexer (which is `tmux`).

### Manual (snakemake only)

Run the entire pipeline via `snakemake` with the following (substitute any
options as desired in the place of `--profile`):

```
snakemake --profile workflow/profiles/<profname> --configfile=config/<confname.yml>
```

Store results in the cloud:

```
dvc commit
dvc push
```

### Manual on Nisaba (snakemake only)

Use the nisaba profile:

```
snakemake --configfile config/dynamic.yml --profile workflow/profiles/nisaba
```

See the [profile](workflow/profiles/nisaba/config.yaml) for slurm/snakemake
options that are set.

The slurm logs will be found in `cluster_logs`, partitioned by each rule.

Note that this command will block the terminal so it is recommended to run it in
the background with `&` or use a multiplexer like `tmux`.

Commit and push as desired:

```
dvc commit
dvc push
```

## Retrieving Results

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

## Pipeline Output

Each entry under `ebm_runs` in the dynamic config corresponds to one EBM run
with its corresponding features and settings to be used. After running the
pipeline, each run should have a directory under `results/ebm` names like
`<git_tag>_<run_entry_name>` where `<git_tag>` is the current tag of the repo
(or the commit if there is none) and `<run_entry_name>` is the key under
`ebm_runs` in the [dynamic config](config/dynamic.yml).

Each directory will contain the input tsv of data used to train the EBM, a
config yml file with all settings used to train the EBM, and python pickles for
the X/Y train/test datasets as well as a pickle for the final model itself.

## Configuring the Pipeline

The configuration is split into [static](config/static.yml) and
[dynamic](config/dynamic.yml) components. The former is for downloaded resources
and parameters used to generate the annotated dataframes.

The latter is for selecting various features within the annotated dataframes,
applying transformations, and selecting hyperparameters when training the EBM
models.
