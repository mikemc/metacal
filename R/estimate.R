#' Estimate bias from control measurements
#'
#' Estimate bias using the compositional least-squares approach described in
#' McLaren, Willis, and Callahan (2019).
#'
#' Name requirements for `observed` and `actual`: The row and column names (for
#' matrices) or taxa and sample names (for phyloseq objects) must match, but
#' can be in different orders.
#'
#' @param observed Abundance matrix of observed compositions
#' @param actual Abundance matrix of actual or reference compositions for the
#'   same samples and taxa in `observed`
#' @param margin Matrix margin that corresponds to observations (samples); 
#'   `1` for rows, `2` for columns
#' @param ... Arguments passed to [center()]
#' 
#' @seealso [center()] [calibrate()] 
#' 
#' @rdname estimate_bias
#' @export
estimate_bias <- function(observed, actual, ...) {
  UseMethod("estimate_bias")
}

#' @rdname estimate_bias
#' @method estimate_bias matrix
#' @export
estimate_bias.matrix <- function(observed, actual, margin, ...) {
  stopifnot(margin %in% 1:2)
  stopifnot(identical(dim(observed), dim(actual)))

  # Align samples and features (taxa)
  stopifnot(setequal(rownames(observed), rownames(actual)))
  stopifnot(setequal(colnames(observed), colnames(actual)))
  observed <- observed[rownames(actual), colnames(actual)]

  error <- observed / actual
  # `center()` requires `error` is oriented with samples as rows
  if (margin == 2) 
    error <- t(error)
  center(error, ...)
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
  # Note: In the case where both are oriented with taxa as rows, we could
  # reduce the total number of transpose operations by 1 by keeping the
  # original orientation

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
