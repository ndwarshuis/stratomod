# giab-ai-ebm

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

NOTE: `config/dynamic.yml` can technically be any name one desires, but it is
currently hardcoded in `dvc.yaml` so it is best to leave as is.

## Running the Pipeline

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

NOTE: `config/dynamic.yml` should only be created on experimental branches (see
above), and also must be specified manually on the command line via
`--configfile`.

### Auto

Run the entire pipeline and store results in the cloud:

```
dvc repro
dvc push
```

### Manual

Run the entire pipeline via snakemake with the following:

```
snakemake -p --verbose -r -j 1 --use-conda --rerun-incomplete --configfile config/dynamic.yml
```

Store results in the cloud:

```
dvc commit
dvc push
```

### Slurm/Nisaba

This repository has a profile to run the pipeline using slurm on the NIST Nisaba
cluster.

Run it with this:

```
snakemake --configfile config/dynamic.yml --profile workflow/profiles/nisaba --cores 4
```

See the [profile](workflow/profiles/nisaba/config.yaml) for slurm/snakemake
options that are set.

The `cores` option will only affect the number of cores used for EBM training;
snakemake will submit up to 500 jobs to run rules in parallel. The slurm logs
will be found in `cluster_logs`, partitioned by each rule.

Note that this command will block the terminal so it is recommended to run it in
the background with `&` or use a multiplexer like `tmux`.

This does not invoke `dvc` so `dvc commit/push` will need to be run manually
to put the results in the cloud.

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
