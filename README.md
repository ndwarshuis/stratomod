# giab-ai-ebm

## deployment

Install the environment to run in snakemake by running this command at the root
of this repo:

```
conda env create -f env.yml
```

## running the pipeline

### manual

Run the entire pipeline via snakemake with the following:

```
snakemake -p --verbose -r -j 1 --use-conda
```

### slurm

Submit a job to run the pipeline via slurm using this [script](snakemake-slurm):

```
sbatch snakemake-slurm
```

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
