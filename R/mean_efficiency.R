#' Compute sample mean efficiencies
#'
#' Given an `mc_bias_fit` object or a vector of taxa efficiencies, compute the
#' mean efficiency in the control samples in the `mc_bias_fit` object or a set
#' of user-supplied samples. This function currently only works on matrices,
#' not phyloseq objects.
#'
#' If `bias` is an `mc_bias_fit` object and `actual` and `observed` are both
#' missing, then the estimated mean efficiency for the control samples used in
#' the `mc_bias_fit` object will be computed. New samples should be given in
#' either `actual` or `observed`, depending on whether they are nominal or
#' calibrated abundances or observed (i.e. uncalibrated) data; observed data
#' will be calibrated first. 
#'
#' The mean efficiencies in a sample with is equal to `weighted.mean(bias, x))`
#' where `bias` is the vector of taxa efficiencies and `x` is the vector of
#' taxa proportions.
#'
#'
#' `r lifecycle::badge("experimental")`
#'
#' @param bias A numeric vector of relative efficiencies or an object of class
#'   'mc_bias_fit' (from which efficiencies will be extracted)
#' @param actual A matrix or NULL
#' @param observed A matrix or NULL
#' @param margin The margin containing the samples/observations
#' 
#' @seealso [estimate_bias()]
#' 
#' @rdname mean_efficiency
#' @export
mean_efficiency <- function(bias, actual, observed, margin) {
  UseMethod("mean_efficiency")
}

#' @rdname mean_efficiency
#' @method mean_efficiency mc_bias_fit
#' @export
mean_efficiency.mc_bias_fit <- function(bias,
                                        actual = NULL, 
                                        observed = NULL, 
                                        margin = NULL) {
  if (!is.null(actual) && !is.null(observed)) {
    stop('At most one of `actual` and `observed` should be non-null.')
  } else if (is.null(actual) && is.null(observed)) {
    stopifnot(is.null(margin))
    actual <- bias$actual
    margin <- 1
  }

  mean_efficiency(
    coef(bias), 
    actual = actual, 
    observed = observed, 
    margin = margin
  )
}

#' @rdname mean_efficiency
#' @method mean_efficiency numeric
#' @export
mean_efficiency.numeric <- function(bias,
                                    actual = NULL, 
                                    observed = NULL, 
                                    margin) {
  # First we need to set `actual` depending on the supplied arguments
  if (!is.null(actual) && !is.null(observed)) {
    stop('Only one of `actual` and `observed` should be non-null.')
  } else if (is.null(actual) && is.null(observed)) {
    stop('Either `actual` or `observed` must be non-null.')
  } else {
    stopifnot(margin %in% 1:2)
    taxa_margin <- setdiff(1:2, margin)
  }
  if (!is.null(observed)) {
    # Note: Consider allowing taxa to be dropped from 'observed' with a
    # message, in the case where 'bias' is named.
    stopifnot(identical(dim(observed)[[taxa_margin]], length(bias)))
    actual <- calibrate(observed, bias, margin = margin)
  }
  stopifnot(identical(dim(actual)[[taxa_margin]], length(bias)))

  # Align taxa order between actual and bias
  if (!is.null(names(bias))) {
    if (!setequal(dimnames(actual)[[taxa_margin]], names(bias)))
      stop("Taxa names in `bias` and `actual` do not match")
    bias <- bias[dimnames(actual)[[taxa_margin]]]
  }

  # Can also do the same thing slightly slower with weighted.mean()
  # apply(bias$actual, 1, function(x) weighted.mean(bias, x)) 
 
  actual %>% 
    # normalize to proportions
    sweep(., 1, rowSums(.), `/`) %>%
    # multiply taxa by their efficiency
    perturb(bias, margin = 1, norm = 'none') %>%
    # sum up within samples
    rowSums
}
