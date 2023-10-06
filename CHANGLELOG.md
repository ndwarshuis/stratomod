# 8.1.1

- fix missing bcftools

# 8.1.0

- add biallelic split option
- fix sloppy env definitions
- fix float bug in truncation annotation for plots

# 8.0.5

- fix missing bsdtar tool

# 8.0.4

- fix TR_length not being a valid feature

# 8.0.3

- fix plot display bug in train reports

# 8.0.2

- fix bug when using only one limit in plotting
- fix absolute paths when using rmd outputs in subworkflows

# 8.0.1

- fix bug in test model specification and variable set check

# 8.0.0

- made configuration file more intuitive/powerful
  - defaulted all unnecessary blocks
  - added lots of documentation
  - merged "static" and "dynamic" configs (static things are part of pydantic
    class)
  - remove run_key level of configuration in each model
  - can now specify repeat_masker features on a per-reference basis
- vcf standardization and correction are now in a single rule
- simplify/speed up vcf parser
- simplify output directory layout

# 7.1.0

- clean up repository and present it more as a tool
- add lots of documentation
- remove all domain specific references to NIST-specific infrastructure

# 7.0.1

- fix chromosome X/Y name bugs

# 7.0.0

- add support for GRCh37
- add support for assigning variables to individual vcfs
- fix lots of bugs
- use pydantic for config validation (and end-user error messages)
- CI/CD improvements (linting)
- cleaner config syntax (split labeled and unlabeled queries, use reference
- sets for filtering, etc)

# 6.2.0

- remove tags from output dirs
- add useful info at top of reports

# 6.1.2

- update snakemake version

# 6.1.1

- allow user to run train without any test datasets

# 6.1.0

- add column for raw absolute index in annotated dataframes which will be
  carried through the train/test split step; this will allow full mapping of
  each row in the EBM test/train dataframes back to the original input dataframe

# 6.0.1

- fix off-by-one error in vcf -> bed parser

# 6.0.0

- simplify homopolymers
  - use individual instead of combined bases
  - only add 1bp slop
- increase memory for summarization steps

# 5.6.0

- fix max ref/alt filter (before it didn't listen to max ref)
- fix continuous vs continuous interaction plots
- add lots of diagnostic outputs to train report
  - global calibration plot
  - fraction of true positive plots
  - error profiles for each feature
- test indel length vs homopol length as an interaction

# 5.5.1

- remove all structural variants by default

# 5.5.0

- save coordinates in ebm input dataframes to facilitate genome browser
  extravaganzas
- banish csv files from this project (and only permit tsv files, the obviously
  superior alternative)
  remove MHC region entirely from analysis (before, anything in the MHC region
  would have been flagged a false positive with most benchmarks)

# 5.4.3

- fix parse errors in MHC region for HG003/4/6/7

# 5.4.2

- fix flipped bed and vcf files

# 5.4.1

- fix typo
- fix benchmark assert rule

# 5.4.0

- add all v4.2.1 giab benchmarks to static config

# 5.3.0

- add intercept to model output/report
- add interpretations for intercept and global scores

# 5.2.0

- use gzip compression (almost) everywhere (according to the testing pipeline
  times it only slows down the pipeline by 4-5%)

# 5.1.2

- raise memory requirements for vcf parsing

# 5.1.1

- rename nisaba resource rule (and actually allocate truckload of memory)

# 5.1.0

- store logs in dvc repo

# 5.0.0

features added
- add options to parse/not parse DP/VAF/GT/GQ fields in the FORMAT column
- add logic to 'train' a model and 'test' it on other VCFs (including those
  without benchmarks/labels)
- ensure all chromosomes are standardized to 1-24 (with X/Y being 23/24)

bugfixes
- use REF instead of ALT for region length
- fix missed HG002 vcf parse fix
- restrict all wildcards to exclude slashes (eg they can only match a single
  dir/filename)

# 4.6.2

- use log10 for transforms instead of natural log

# 4.6.1

- fix super silly bug that scrambled the tandem repeat length feature
- don't include tandem repeats whose unit size is 1bp (these are homopolymers)

# 4.6.0

- add 2.7 xy benchmark to static configuration
- remove tbi_url config options (we calculate tbi on the fly now)

# 4.5.0

- make dvc repro actually nice to use and configure

# 4.4.2

- fix vcfeval tmp dir clash

# 4.4.1

- update pipeline to work with dockerized CI/CD pipeline

# 4.4.0

- add option to split continuous EBM tree plots by missing and non-missing
- make all EBM plots the same width

# 4.3.4

- fix vcfeval memory limit
- fix thread typo

# 4.3.3

- fix vcfeval mem management for large datasets
- specify threads manually for vcfeval and ebm train

# 4.3.2

- fix slurm out of memory errors
- fix quote errors when making cluster output files with weird characters

# 4.3.1

- dynamically allocate memory for large jobs to keep slurm happy

# 4.3.0

- remove old slurm script
- add profile for running pipeline on clusters (which really means Nisaba)
- update docs for running on clusters
- add resource annotations to memory-guzzling snakemake rules to prevent cluster
  admins from booting us

# 4.2.0

- add option to set truncation in config
- add option to set plot type in config
- add captions for truncated ranges
- add vcf input key legend at the top of each EBM report
- fix flipped legend in bivariate cat/cat plots

# 4.1.1

- fix segdup feature errors (off by one columns)
- add benchmarking directives to key rules for debugging

# 4.1.0

- use config to control feature naming (mostly)
- add checks for feature names in initial validation steps

# 4.0.0

- use dynamic-testing.yml only for testing
- improve CI pipeline to test snakemake global and rule conda envs as well as
  resources and results computation for chr21
- remove superfluous reference download for homopolymers
- make curl commands quieter
- add target to download all resources (mostly for testing)

# 3.2.0

- add working CI pipeline with super lame test

# 3.1.0

- add option to change name of features in config file
- add 'binary' feature transformation option

# 3.0.0

- rename features to have standard prefixes per group
- rename `include_filtered` to `filtered_are_candidates` (config option)
- make mappability features independent (eg no overlap)
- remove `consensusSize` feature
- make static config annotations dependent on reference (futureproofing)

# 2.5.1

- make Rmd files run faster with large files

# 2.5.0

- add logging to all (important steps)
- add filtering to annotations code to make testing not take insanely long

# 2.4.1

- add features overview file

# 2.4.0

- add AT/GC gap length/frac homopolymer features
- only sort/filter the massive homopolymer raw dataframe once

# 2.3.0

- add AT/GC perfect homopoly fraction
- fix annotations using -1 for missing values in some cases

# 2.2.2

- add option to include all interactions with one feature
- fix (almost) all things regarding interaction plots

# 2.2.1

- update dvcignore

# 2.2.0

- add error bars to output plots
- remove 'config' version tag and use pure semantic versioning for the pipeline

# 2.1.0-0

- add QUAL score
- add config option to control if filtered entries are counted in the query
- add draft benchmark to static config

# 2.0.0-0

- add option to concatenate multiple input dataframes as needed for EBM training
  set

# 1.5.0-3

- use DVC to track results

# 1.4.0-3

- include FN and FP and include logic for mapping these to negative class

# 1.3.0-3

- make ebm reporting more modular

# 1.2.0-3

- add transformed univariate model plots

# 1.1.2-3

- fix html output for models

# 1.1.1-3

- remove all log transform stuff from config (we probably don't need this)
- remove all intersections from config (too complex at this point)
- add precision FDA VCF as input to compare with the Illumna VCF

# 1.1.1-2

- fix colors/ranges on correlation plots

# 1.1.0-2

- add input dataframe summaries
  - tables with key feature statistics
  - information ranking
  - true colocalization matrix with asymmetric jaccard index
  - correlation plots for features that perfectly colocalize
  - feature plots with transforms
- add model performance summaries
- new features
  - combined AT and GC percents (min/max/med) from tandem repeats
  - LINE families
