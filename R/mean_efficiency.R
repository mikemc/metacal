#' Compute sample mean efficiencies
#'
#' Given a matrix of 'observed' or 'actual' (e.g. calibrated) abundance
#' profiles and a vector of relative efficiencies, compute the estimated mean
#' efficiency for each sample.
#'
#' The mean efficiency of a sample equals `weighted.mean(bias, y)`, where
#' `bias` is the vector of taxa efficiencies and `y` is the vector of taxa
#' proportions.
#'
#' The `type %in% c('actual', 'observed')` specifies whether the data is
#' 'actual' (i.e. nominally known or calibrated abundances) or 'observed' (i.e.
#' uncalibrated) data. 'Observed' data is calibrated to determine `y` prior to
#' computing the mean efficiency.
#'
#' If `x` is an `mc_bias_fit` object, then the efficiencies will be extracted
#' with `coef(x)`. An abundance matrix can optionally be given with `newdata`;
#' otherwise, the mean efficiency will be computed for the control samples in
#' `x`.
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param x A matrix, a phyloseq object with an otu_table, or an object of
#'   class 'mc_bias_fit'
#' @param bias A (possibly named) numeric vector of relative efficiencies
#' @param newdata NULL or an abundance matrix
#' @param margin The margin containing the samples (observations)
#' @param type 'actual' or 'observed', indicating the type of abundance
#'   profiles in the matrix `x` or `newdata`
#' 
#' @seealso [estimate_bias()]
#' 
#' @rdname mean_efficiency
#' @export
# mean_efficiency <- function(bias, actual, observed, margin) {
#   UseMethod("mean_efficiency")
# }
mean_efficiency <- function(x, ...) {
  UseMethod("mean_efficiency")
}

#' @rdname mean_efficiency
#' @method mean_efficiency matrix
#' @export
mean_efficiency.matrix <- function(x,
                                   bias,
                                   margin,
                                   type) {
  # First, get to point where `x` is a matrix with the 'actual' (unbiased)
  # profiles with samples as rows, and the taxa (columns) are properly aligned
  # with `bias`.
  if (!margin %in% 1:2)
    stop("`margin` must be 1 or 2.")
  if (margin == 2) {
    x <- t(x)
    margin <- 1
  }
  taxa_margin <- 2
  if (type == 'observed') {
    # Note: Consider allowing taxa to be dropped from 'observed' with a
    # message, in the case where 'bias' is named.
    stopifnot(identical(dim(x)[[taxa_margin]], length(bias)))
    x <- calibrate(x, bias, margin = margin, norm = 'none')
  } else if (type != 'actual') {
    stop("`type` must be 'actual' or 'observed'.")
  }
  # TODO: Add ability to extract efficiencies from mc_bias_fit objects
  # Align taxa order between x and bias
  if (!is.null(names(bias))) {
    if (!setequal(dimnames(x)[[taxa_margin]], names(bias)))
      stop("Taxa names in `bias` and `actual` do not match")
    bias <- bias[dimnames(x)[[taxa_margin]]]
  }

  # Can also do the same thing slightly slower with weighted.mean()
  # apply(bias$actual, 1, function(x) weighted.mean(bias, x)) 
 
  x %>% 
    # normalize to proportions
    sweep(., 1, rowSums(.), `/`) %>%
    # multiply taxa by their efficiency
    perturb(bias, margin = 1, norm = 'none') %>%
    # sum up within samples
    rowSums
}

#' @rdname mean_efficiency
#' @method mean_efficiency otu_table
#' @export
mean_efficiency.otu_table <- function(x, bias, type) {
  mean_efficiency(
    as(x, "matrix"), 
    bias,
    margin = 1 + taxa_are_rows(x),
    type = type
  )
}

#' @rdname mean_efficiency
#' @method mean_efficiency phyloseq
#' @export
mean_efficiency.phyloseq <- function(x, bias, type) {
  mean_efficiency(otu_table(x), bias, type = type)
}

#' @rdname mean_efficiency
#' @method mean_efficiency mc_bias_fit
#' @export
mean_efficiency.mc_bias_fit <- function(x,
                                        newdata = NULL, 
                                        margin = NULL,
                                        type = NULL) {
  bias = coef(x)

  if (is.null(newdata)) {
    stopifnot(is.null(margin))
    x <- x$actual
    type <- 'actual'
    margin <- 1
  } else {
    if (! type %in% c('actual', 'observed')) {
      stop("`type` must be 'actual' or 'observed'.")
    }
    if (! margin %in% 1:2) {
      stop("`margin` must be 1 or 2.")
    }
  }

  mean_efficiency(
    x = x,
    bias = bias,
    type = type,
    margin = margin
  )
}

