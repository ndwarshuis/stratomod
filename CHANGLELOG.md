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
