# StratoMod

A model-based tool to quantify the difficulty of calling a variant given genomic
context.

## Background

Intuitively we understand that acurately calling variants in a genome can be
more or less difficult depending on the context of that variant. For example,
many sequencing technologies have higher error rates in homopolymers, and this
error rate generally increases as homopolymers get longer. However, precisely
quantifying the relationship between these errors, the length of the
homopolymer, and the impact on the resulting variant call remain challenging.
Analogous arguments can be drawn for other "repetitive" regions in the genome,
such as tandem repeats, segmental duplications, transposable elements, and
diffitult-to-map regions.

The solution we present here is to use an interpretable modeling framework
called [explainable boosting machines](https://github.com/interpretml/interpret)
to predict variant calling errors as a function of genomic features (eg, whether
or not the variant in a tandem repeat, homopolymer, etc). The interpretability
of the model is important for allowing end users to understand the relationship
each feature has to the prediction, which facilitates understanding (for
example) at what lengths of homopolymers the likelihood of incorrectly calling a
variant drastically increases. This precision is an improvement over [existing
methods](https://github.com/ndwarshuis/giab-strats-smk) we have developed for
statifying the genome by difficulty into discrete bins. Furthermore, this
modeling framework allows understanding of interactions between different
genomic contexts, which is important as many repetative characterstics do not
exist in isolation.

We anticipate `StratoMod` would be useful for both method developers and
clinicians who wish to better understand variant calling error modalities. In
the case of method development, `StratoMod` can be used to accurately compare
the error modalities of different sequencing technologies. For clinicians, this
can be used for determining in which regions/genes (which may be clinically
interesting for a given study) variant errors are likely to occur, which may in
turn inform which technologies should be employed and/or other mitigation
strategies should be used.

Further information can be found in our
[preprint](https://www.biorxiv.org/content/10.1101/2023.01.20.524401v1).

## User Guide

### Pipeline steps

1. Compare user-supplied query vcf with GIAB benchmark vcf to produce labels
   (true positive, false positive, false negative). The labels comprise the
   dependent variable used in model training downstream.

2. Intersect comparison output labels with genomic features to produce the
   featurs (independent variables) used in model training.

3. Train the EBM model with random holdout for testing

4. If desired, test the model on other query vcfs (which may or may not also be
   labelled with a benchmark comparison).
   
5. Inspect the output features (plots showing the profile of each feature and
   its effect on the label).
   
NOTE: currently only two labels can be compared at once given that we used a
binary classifier. This means either one of the three labels must be omitted or
two need to be combined into one label.

### Data Inputs

The only mandatory user-supplied data requrired to run is a query vcf.
Optionally one can supply other vcfs for testing the model.

Unless one is using esoteric references or benchmarks, the pipeline is
preconfigured to retrieve commonly-used data defined by flags in the
configuration. This includes:
- a GIAB benchmark, including the vcf, bed, and reference fasta
- reference-specific bed files which will provide "contextual" features for each
  variant call, including:
  - difficult-to-map regions (GIAB stratification bed file)
  - segmental duplications (UCSC superdups database)
  - tandem repeats (UCSC simple repeats database)
  - transposable elements (USCS Repeat Masker)

### Installation

This assumes the user has a working `conda` or `mamba` installation.

Run the following to set up the runtime environment.

```
mamba env create -f env.yml
```

### Configuration

A sample configuration file may be found in `config/dynamic-testing.yml` which
may be copies as a starting point. This file is heavily annotated to explain all
the options/flags and their purpose.

### Running

Execute the pipeline using snakemake:


```
snakemake --use-conda -c <num_cores> --rerun-incomplete --configfile=config/<confname.yml> all
```

## Pipeline Output

Each entry under `ebm_runs` in the dynamic config corresponds to one EBM run
with its corresponding features and settings to be used. After running the
pipeline, each run should have a directory under `results/model`.

Each directory will contain the input tsv of data used to train the EBM, a
config yml file with all settings used to train the EBM, and python pickles for
the X/Y train/test datasets as well as a pickle for the final model itself.

Within the run directory will also be a `test` directory which will contain all
test runs (eg the results of the model test and the input data used for the
test).

## Developer Guide

"Development" refers to adding features/altering code. For configuring and
committing experiments, see above.

### Environments

By convention, the conda environment specified by `env.yml` only has runtime
dependencies for the pipeline itself.

To install development environments, run the following:

```
./setup_dev.sh
```

In addition to creating new environments, this script will update existing
ones if they are changed during development.

Note that scripts in the pipeline are segregated by environment in order to
prevent dependency hell while maintaining reproducible builds. When editing, one
will need to switch between environments in the IDE in order to benefit from the
features they provide. Further details on which environments correspond to which
files can be found in `workflow/scripts`.

Note that this will only install environments necessary for running scripts (eg
rules with a `script` directive).

### Linting

All python code should be error free when finalizing any new features. Linting
will be performed automatically as part of the CI/CD pipeline, but to run it
manually, invoke the following:

```
./lint.sh
```

This assumes all development environments are installed (see above).

### New Feature Workflow

There are two main development branches: `master` and `develop`.

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
