# run just this file:
# devtools::test_active_file(here::here("tests", "testthat", "test-ratios.R"))

data("enterotype", package = "phyloseq")

# Test for issue #18 fix
test_that("Maintains name with just 1 item on opposing margin; throws error if 1 item on given margin", {
  p0 <- enterotype %>%
    phyloseq::prune_taxa(taxa_names(.)[3], .) %>%
    phyloseq::prune_samples(sample_names(.)[1:5], .)
  p1 <- p0 %>% pairwise_ratios("samples")
  expect_equal(taxa_names(p1), taxa_names(p0))
  expect_error(p0 %>% pairwise_ratios("taxa"))

  p0 <- enterotype %>%
    phyloseq::prune_taxa(taxa_names(.)[3:5], .) %>%
    phyloseq::prune_samples(sample_names(.)[5], .)
  p1 <- p0 %>% pairwise_ratios("taxa")
  expect_equal(sample_names(p1), sample_names(p0))
  expect_error(p0 %>% pairwise_ratios("samples"))
})
