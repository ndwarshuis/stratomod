# 1.1.0-3

- remove all log transform stuff from config (we probably don't need this)
- remove all intersections from config (too complex at this point)
- add precision FDA VCF as input to compare with the Illumna VCF

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
