context("Methods for computing compositional means")

# Setup -----------------------------------------------------------------------

library(tibble)
library(tidyr)
library(dplyr)

set.seed(1)
# Bias for 5 fake taxa:
K <- 5
taxa <- paste0("T", seq(K))
bias <- tibble(Taxon = taxa, Bias = exp(rnorm(K, 0, 2))) %>%
    mutate(
        Bias = center_elts(Bias),
    )
# Fake experimental dataset 1
N <- 10
samples <- paste0("S", seq(N))
exp1 <- crossing(Sample = samples, Taxon = taxa) %>%
    left_join(bias, by = "Taxon")
exp1 <- exp1 %>%
    mutate(
        Actual = rexp(n(), 1),
        Noise = exp(rnorm(n(), mean = log(1), sd = 0.5)),
        Observed = Actual * Bias * Noise
    ) %>%
    group_by(Sample) %>%
    # mutate_at(vars(Actual, Noise, Observed), center_elts) %>%
    ungroup
# Fake experimental dataset 2
exp2 <- exp1 %>%
    mutate(
        Actual = Actual * purrr::rbernoulli(n(), p = 0.66),
        Observed = Actual * Bias * Noise
    )
# Fake experimental dataset 3
exp3 <- exp1 %>%
    mutate(
        Actual = case_when(
            Sample <= "S5" & Taxon <= "T3" ~ Actual,
            Sample <= "S5" & Taxon > "T3" ~ 0,
            Sample > "S5" & Taxon <= "T3" ~ 0,
            Sample > "S5" & Taxon > "T3" ~ Actual,
        ),
        Observed = Actual * Bias * Noise
    )
# Get compositional errors in matrix form
d1 <- exp1 %>%
    transmute(Sample, Taxon, Error = Observed / Actual) %>%
    build_matrix(Sample, Taxon, Error)
d2 <- exp2 %>%
    transmute(Sample, Taxon, Error = Observed / Actual) %>%
    build_matrix(Sample, Taxon, Error)
d3 <- exp3 %>%
    transmute(Sample, Taxon, Error = Observed / Actual) %>%
    build_matrix(Sample, Taxon, Error)
# Example set of weights
w <- rmultinom(1, N, rep(1, N))[,1]

# Deterministic datasets
d4 <- exp1 %>%
    build_matrix(Sample, Taxon, Bias)
d5 <- exp2 %>%
    transmute(Sample, Taxon, Error = Actual * Bias / Actual) %>%
    build_matrix(Sample, Taxon, Error)

# TODO: Dataset w/ covariance in the errors

# Mean ------------------------------------------------------------------------

test_that("all methods agree when all samples have all taxa", {
    expect_equal(center(d1, method = "proj"), center(d1, method = "gm"))
    expect_equal(center(d1, method = "proj"), center(d1, method = "rss"),
        tolerance = 1e-5)
})

test_that("the 'proj' and 'rss' methods agree when samples have different taxa, but sufficient overlap", {
    expect_equal(center(d2, method = "proj"), center(d2, method = "rss"),
        tolerance = 1e-5)
})

test_that("the 'gm' method returns all NaNs when some taxa-sample observations are missing", {
    expect_equal(center(d2, method = "gm"), 
        rlang::rep_named(colnames(d2), NaN))
})

test_that("the estimate does not depend on first normalizing the `weights`", {
    expect_equal(
        center(d2, weights = w, method = "rss"), 
        center(d2, weights = w / sum(w), method = "rss"), 
        tolerance = 1e-5)
    expect_equal(
        center(d2, weights = w, method = "proj"), 
        center(d2, weights = w / sum(w), method = "proj"))
    expect_equal(
        center(d2, weights = w, method = "gm"), 
        center(d2, weights = w / sum(w), method = "gm"))
})

test_that("sampling rows or using `weights` gives an equal estimate", {
    rows <- sample(N, N, replace = TRUE)
    w.rows <- table(factor(rows, seq(N))) %>% c
    expect_equal(
        center(d1[rows,], method = "gm"), 
        center(d1, weights = w.rows, method = "gm"))
    expect_equal(
        center(d2[rows,], method = "proj"), 
        center(d2, weights = w.rows, method = "proj"))
    expect_equal(
        center(d2[rows,], method = "rss"), 
        center(d2, weights = w.rows, method = "rss"))
})


# Additional tests to add
# - Check various in scale and out scale combinations
# - Check behavior in the exp3 case (different cliques of taxa)
# - Check that the projection method gives the right answer in the
# deterministic missing observations case
