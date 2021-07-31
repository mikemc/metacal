# metacal (development version)

* `build_matrix()` now works correctly on grouped tibbles

* `pairwise_ratios()` now correctly handles phyloseq objects with just 1 sample
or taxon, and is properly exported

* `estimate_bias()` now allows `observed` to have extra samples and taxa not in `actual`, by automatically subsetting to those in `actual`

# metacal 0.2.0

## New features

### New estimation and calibration functions

`estimate_bias()` and `calibrate()` provide easy-to-use high-level interfaces to the original metacal bias-estimation and calibration method.
Their use is illustrated in the new tutorial.

### New utility functions

* `pairwise_ratios()` allows computing ratios between pairs of taxa and/or samples.

* `perturb()` applies a compositional perturbation to all observations in an abundance matrix (or phyloseq object).

### Support for phyloseq objects

The new functions `estimate_bias()`, `perturb()`, `calibrate()`, and `pairwise_ratios()` all work with objects from the [phyloseq](https://joey711.github.io/phyloseq/) package as well as plain matrices. 
Phyloseq is now a required dependency, though may be made optional in the future.

### New tutorial

The [new tutorial](https://mikemc.github.io/metacal/articles/tutorial.html) demonstrates the new estimation and calibration functions on a new dataset from [Leopold and Busby (2020)](https://doi.org/10.1016/j.cub.2020.06.011).

## Minor fixes

* Fixed failure of `coocurrence_matrix()` when names were missing (#5).

# metacal 0.1.0

Direct implementation of the methods described in [McLaren MR, Willis AD, Callahan BJ (2019)](https://elifesciences.org/articles/46923) and used for [the analysis of that paper](https://github.com/mikemc/mgs-bias-manuscript).
