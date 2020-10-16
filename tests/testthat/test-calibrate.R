# run just this file:
# devtools::test_file(here::here("tests", "testthat", "test-calibrate.R"))
#
# Covers estimate_bias(), perturb(), calibrate()

# Setup -----------------------------------------------------------------------

test_xy <- function(n_taxa, n_samples) {
  set.seed(42)
  # random abundance matrix
  x <- sample(100, n_taxa*n_samples, replace = TRUE) %>% 
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

  expect_equal(bias, estimate_bias(observed, actual))
  expect_equal(bias, estimate_bias(observed, actual %>% t))
})

test_that("`calibrate()` and `perturb()` are inverse operations", {
  expect_equal(
    otu,
    otu %>% perturb(y, norm = "none") %>% calibrate(y, norm = "none")
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
