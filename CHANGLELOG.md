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
