# run just this file:
# devtools::test_active_file(here::here("tests", "testthat", "test-utilities.R"))

# Setup: Load Brooks files
x <- system.file("extdata", "brooks2015-actual.csv", package = "metacal") %>%
  readr::read_csv() %>%
  dplyr::arrange(Sample) # Not strictly needed since already in this order
sam <- system.file("extdata", "brooks2015-sample-data.csv", package = "metacal") %>%
  readr::read_csv() 

test_that("`as_matrix()` and `build_matrix()` work as inverses", {
  y <- x %>%
    as_matrix(Sample)
  z <- x %>%
    tidyr::pivot_longer(-Sample) %>%
    build_matrix(Sample, name, value)
  expect_equal(y, z)
})

test_that("`build_matrix()` works on grouped variables", {
  y <- x %>%
    dplyr::left_join(sam, by = "Sample") %>%
    tidyr::pivot_longer(-colnames(sam), names_to = ".otu")
  z <- y %>%
    dplyr::group_by(Plate, .otu) %>%
    dplyr::summarize(across(value, sum), .groups = "drop_last") %>%
    build_matrix(Plate, .otu, value)
  expect_equal(
    rownames(z), 
    sam$Plate %>% unique %>% as.character %>% sort
  )
  expect_equal(
    colnames(z), 
    y$.otu %>% unique %>% sort
  )
})
