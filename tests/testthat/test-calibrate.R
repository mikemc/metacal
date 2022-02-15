# run just this file:
# devtools::test_active_file(here::here("tests", "testthat", "test-calibrate.R"))
#
# Covers estimate_bias(), perturb(), calibrate()

# Setup -----------------------------------------------------------------------

test_xy <- function(n_taxa, n_samples) {
  set.seed(42)
  # random abundance matrix with zeros
  x <- sample(0:5, n_taxa*n_samples, replace = TRUE) %>% 
    matrix(nrow = n_taxa)
  rownames(x) <- paste0("t", seq(n_taxa))
  colnames(x) <- paste0("s", seq(n_samples))
  # perturbation vector
  y <- seq(n_taxa)

  list(x = x, y = y)
}

l <- test_xy(5, 10)
x <- l$x; y <- l$y
otu <- otu_table(x, taxa_are_rows = TRUE)

# Tests -----------------------------------------------------------------------

test_that("`estimate_bias()` correctly recovers a deterministic perturbation", {
  actual <- otu
  bias <- rlang::set_names(center_elts(y), phyloseq::taxa_names(otu))
  observed <- perturb(otu, bias, norm = "keep")
  # Mix up the sample and taxa order
  observed <- observed[rownames(observed) %>% rev, colnames(observed) %>% rev]

  expect_equal(bias, estimate_bias(observed, actual) %>% coef)
  expect_equal(bias, estimate_bias(observed, actual %>% t) %>% coef)

  # non-zero values is observed that are zero in actual should be automatically
  # zeroed with a message
  observed1 <- otu_table(observed + 23 * (observed == 0), taxa_are_rows = TRUE)
  expect_message(fit <- estimate_bias(observed1, actual, boot = TRUE, times = 2))
  expect_equal(bias, coef(fit))

  # Matrices stored in `fit` should have samples as rows, and so work without
  # error in perturb() with `margin = 1`.
  perturb(fit$actual, fit$estimate, margin = 1, norm = "close")

  # Check perturb on mc_bias_fit's
  x <- perturb(fit, 1/fit$estimate)
  expect_equal(x$estimate / x$estimate, apply(x$bootreps, 2, mean))
  expect_error(perturb(fit, 1:3))

  # NOTE: doesn't quite fit this test title
  # Test that calibrate works on mc_bias_fit objects
  # fit <- estimate_bias(observed, actual)
  expect_equal(
    calibrate(observed, fit),
    calibrate(observed, coef(fit))
  )

  # Test mean efficiency computation
  me1 <- apply(fit$actual, 1, function(x) weighted.mean(coef(fit), x)) 
  me2 <- mean_efficiency(fit$actual, coef(fit), margin = 1, type = 'actual')
  me3 <- mean_efficiency(fit)
  expect_equal(me1, me3)
  expect_equal(me2, me3)
  me4 <- actual %>%
    phyloseq::transform_sample_counts(close_elts) %>%
    perturb(coef(fit), norm = 'none') %>%
    phyloseq::sample_sums() %>%
    .[names(me1)]
  expect_equal(me1, me4)

  # can also compute within `calibrate()`
  sam <- data.frame(
    dummy_var = seq(phyloseq::nsamples(observed)),
    row.names = sample_names(observed)
    ) %>%
    phyloseq::sample_data()
  cal <- phyloseq::phyloseq(observed, sam) %>%
    calibrate(coef(fit), mean_name = 'silly_name')
  me5 <- sample_data(cal)$silly_name
  names(me5) <- sample_names(cal)
  me5 <- me5[names(me1)]
  expect_equal(me1, me4)
})

test_that("`estimate_bias()` silently drops taxa and samples not in 'actual'", {
  actual <- otu %>%
    phyloseq::prune_taxa(taxa_names(.)[-3], .) %>%
    phyloseq::prune_samples(sample_names(.)[-4], .)
  bias <- rlang::set_names(center_elts(y), phyloseq::taxa_names(otu))
  observed <- perturb(otu, bias, norm = "keep")
  # Mix up the sample and taxa order
  observed <- observed[rownames(observed) %>% rev, colnames(observed) %>% rev]

  expect_equal(
    bias[-3] %>% center_elts, 
    estimate_bias(observed, actual) %>% coef
  )
})

test_that("`calibrate()` and `perturb()` are inverse operations", {
  expect_equal(
    otu,
    otu %>% perturb(y, norm = "none") %>% calibrate(y, norm = "none")
  )
  mat <- otu %>% as("matrix")
  expect_equal(
    mat,
    mat %>% perturb(y, 2, norm = "none") %>% calibrate(y, 2, norm = "none")
  )
  # Should also work if y is named and has a different order
  y1 <- rlang::set_names(y, taxa_names(otu)) %>% rev
  expect_equal(
    otu,
    otu %>% perturb(y1, norm = "none") %>% calibrate(y, norm = "none")
  )
})

test_that("`calibrate()` (and `perturb()`) normalizations work", {
  expect_equal(
    rlang::set_names(rep(1, phyloseq::nsamples(otu)), phyloseq::sample_names(otu)),
    calibrate(otu, y, norm = "close") %>% phyloseq::sample_sums()
  )
  expect_equal(
    otu %>% sample_sums,
    calibrate(otu, y, norm = "keep") %>% phyloseq::sample_sums()
  )
})

test_that("taxa disagreements are properly handled", {
  y1 <- rlang::set_names(y, phyloseq::taxa_names(otu)) %>% rev
  y2 <- rlang::set_names(y, letters[seq_along(y)])
  # calibrate should work as long as some taxa names are shared with `bias`
  cal <- calibrate(otu, y1[1:4])
  expect_identical(taxa_names(cal), intersect(taxa_names(otu), names(y1)[1:4]))
  expect_error(calibrate(otu, y[1:4]))
  expect_error(calibrate(otu, y2))
  expect_error(calibrate(otu, y2[1:4]))
})
