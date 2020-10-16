# metacal

The metacal package provides tools for bias estimation and calibration in
marker-gene and metagenomics sequencing experiments. It implements the methods
described in [McLaren MR, Willis AD, Callahan BJ
(2019)](https://elifesciences.org/articles/46923) and is used for the
analysis associated with that manuscript, available at the [manuscript's
repository](https://github.com/mikemc/mgs-bias-manuscript).

## Installation

Install the development version of metacal from from GitHub,

``` r
# install.packages("devtools")
devtools::install_github("mikemc/metacal")
```

## Usage

See the [package tutorial](https://mikemc.github.io/metacal/articles/tutorial.html)
for a demonstration of how to estimate bias from control samples with known
composition (i.e., mock community samples), and how to calibrate the relative
abundances in unknown samples of the taxa that were in the controls.

The primary utility of this package is quantitatively estimating the bias of
protocols in quality control experiments, where samples with known composition
are measured or samples with unknown composition are measured by multiple
protocols.

It is currently not possible to calibrate the composition of a natural
community without making strong and untested assumptions about bias being the
same for constructed and natural samples and about the efficiencies of taxa not
in the controls (e.g., approximating them by that of the closest relative or
the average efficiency). For this and other limitations described in the
Discussion of our manuscript, calibration as a practical method to obtain
quantitatively accurate composition measurements is not currently feasible
using this or any package. However, calibration using a hypothesized bias
(perhaps partially informed by experimental measurement) can still be useful to
analyze the sensitivity of downstream results to bias, a use case we will
illustrate in a future vignette.
