#' Estimate bias from control measurements
#'
#' Estimate bias using the compositional least-squares approach described in
#' McLaren, Willis, and Callahan (2019).
#'
#' Bias is estimated by applying [center()] to the compositional error matrix
#' defined by `observed/actual`, which requires that `observed` and `actual`
#' are non-zero for the same sample-taxa pairs. For convenience, this
#' function will automatically set values in `observed` to 0 whose
#' corresponding entries are 0 in `actual`, but it is up to you to replace 0
#' values in `observed` with a non-zero value (such as a pseudocount).
#'
#' Name requirements for `observed` and `actual`: The row and column names (for
#' matrices) or taxa and sample names (for phyloseq objects) must match, but
#' can be in different orders.
#'
#' @param observed Abundance matrix of observed compositions.
#' @param actual Abundance matrix of actual or reference compositions for the
#'   same samples and taxa in `observed`.
#' @param margin Matrix margin that corresponds to observations (samples); 
#'   `1` for rows, `2` for columns.
#' @param ... Arguments passed to the matrix method.
#' @param boot Whether to perform bootstrapping.
#' @param times Number of bootstrap replicates.
#'
#' @return A `mc_bias_fit` object with [coef()], [fitted()], [residuals()], and
#'   [summary()] methods.
#' 
#' @seealso [center()] [calibrate()] 
#' 
#' @rdname estimate_bias
#' @export
#'
#' @examples
#' # Load data from the cellular mock communities of Brooks et al 2015
#' dr <- system.file("extdata", package = "metacal")
#' list.files(dr)
#' actual <- file.path(dr, "brooks2015-actual.csv") |>
#'   read.csv(row.names = "Sample") |>
#'   as("matrix")
#' observed <- file.path(dr, "brooks2015-observed.csv") |>
#'   read.csv(row.names = "Sample") |>
#'   subset(select = - Other) |>
#'   as("matrix")
#' sam <- file.path(dr, "brooks2015-sample-data.csv") |> read.csv()
#'
#' # Estimate bias with bootstrapping for error estimation
#' mc_fit <- estimate_bias(observed, actual, margin = 1, boot = TRUE)
#' summary(mc_fit)
estimate_bias <- function(observed, actual, ...) {
  UseMethod("estimate_bias")
}

#' @rdname estimate_bias
#' @method estimate_bias matrix
#' @export
estimate_bias.matrix <- function(observed, 
                                 actual, 
                                 margin, 
                                 boot = FALSE, 
                                 times = 1000) {
  stopifnot(margin %in% 1:2)
  stopifnot(identical(dim(observed), dim(actual)))

  # Align samples and taxa
  stopifnot(setequal(rownames(observed), rownames(actual)))
  stopifnot(setequal(colnames(observed), colnames(actual)))
  observed <- observed[rownames(actual), colnames(actual)]
  # Standardize on samples as rows
  if (margin == 2) {
    observed <- t(observed)
    actual <- t(actual)
  }
  # Handle spurious non-zero observations
  n <- sum(observed[actual == 0] > 0)
  if (n > 0) {
    message(paste("Zeroing", n, "values in `observed`"))
    observed[actual == 0] <- 0
  }

  error <- observed / actual
  estimate <- center(error)

  if (boot)
    bootreps <- bootrep_center(error, R = times)
  else 
    bootreps <- NULL

  structure(class = "mc_bias_fit",
    list(
      observed = observed,
      actual = actual,
      estimate = estimate,
      bootreps = bootreps
    )
  )
}

#' @rdname estimate_bias
#' @method estimate_bias otu_table
#' @export
estimate_bias.otu_table <- function(observed, actual, ...) {
  stopifnot(setequal(taxa_names(observed), taxa_names(actual)))
  stopifnot(setequal(sample_names(observed), sample_names(actual)))

  # Orient both tables with samples as rows before passing to matrix method,
  # which will adjust row and column order.
  if (taxa_are_rows(observed)) observed <- t(observed)
  if (taxa_are_rows(actual)) actual <- t(actual)

  estimate_bias.matrix(
    as(observed, "matrix"), 
    as(actual, "matrix"),
    margin = 1,
    ...
  )
}

#' @rdname estimate_bias
#' @method estimate_bias phyloseq
#' @export
estimate_bias.phyloseq <- function(observed, actual, ...) {
  estimate_bias.otu_table(otu_table(observed), otu_table(actual), ...)
}

setMethod("estimate_bias", c("matrix", "matrix"), estimate_bias.matrix)
setMethod("estimate_bias", c("otu_table", "otu_table"), estimate_bias.otu_table)
setMethod("estimate_bias", c("phyloseq", "phyloseq"), estimate_bias.phyloseq)
setMethod("estimate_bias", c("otu_table", "phyloseq"), estimate_bias.phyloseq)
setMethod("estimate_bias", c("phyloseq", "otu_table"), estimate_bias.phyloseq)

#' @export
fitted.mc_bias_fit <- function(object) {
  # TODO: give a norm "keep" option that keeps the totals in _observed_
  # Or, could be a "scale" param with a "counts" or "original" or "observed"
  # option
  perturb(object$actual, object$estimate, margin = 1, norm = "close")
}

#' @export
residuals.mc_bias_fit <- function(object) {
  # For now, give residuals on the proportion scale; in future, give options
  # for other types.
  observed <- object$observed %>% sweep(., 1, rowSums(.), `/`)
  fitted <- fitted(object)
  observed - fitted
}

#' @export
coef.mc_bias_fit <- function(object) {
  object$estimate
}

#' @export
print.mc_bias_fit <- function(x) {
  cat('A metacal bias fit.', fill = TRUE)
  cat(fill = TRUE)
  cat('Estimated relative efficiencies:', fill = TRUE)
  print(x$estimate)
  if (!is.null(x$bootreps)) {
    cat(fill = TRUE)
    cat('Contains', nrow(x$bootreps), 'bootstrap replicates.', fill = TRUE)
  }
  invisible(x)
}

#' @export
summary.mc_bias_fit <- function(object) {
  s <- structure(list(), class = "mc_bias_fit_summary")

  # If no taxon names, use the integer index
  taxon_names <- names(object$estimate)
  if (is.null(taxon_names))
    taxon_names <- seq_along(object$estimate)

  # Create a tibble with the estimated coefficients, with means and standard
  # errors from the bootreps (if they exist)
  s$coefficients <- tibble::tibble(
    taxon = taxon_names,
    estimate = object$estimate,
  )
  if (!is.null(object$bootreps)) {
    s$coefficients <- s$coefficients %>% 
      tibble::add_column(
        gm_mean = object$bootreps %>% apply(2, gm_mean),
        gm_se = object$bootreps %>% apply(2, gm_sd)
      )
    s$n_bootreps <- nrow(object$bootreps)
  }

  s
}

#' @export
print.mc_bias_fit_summary <- function(x) {
  cat('Summary of a metacal bias fit.', fill = TRUE)
  cat(fill = TRUE)
  cat('Estimated relative efficiencies:', fill = TRUE)
  print(x$coefficients)
  if (!is.null(x$n_bootreps)) {
    cat(fill = TRUE)
    cat('Geometric standard error estimated from', x$n_bootreps, 
      'bootstrap replicates.', 
      fill = TRUE)
  }
  invisible(x)
}
