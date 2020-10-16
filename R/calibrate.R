# perturb ---------------------------------------------------------------------

#' Compositionally perturb a relative-abundance matrix
#'
#' The perturbation of two compositional vectors `v` and `y` is given (up to
#' compositional equivalence) by the elementwise product `v * y`.
#' 
#' Let `w = v * y`, where `v` is a row or column of `x`. The normalization
#' options specified by the `norm` argument are
#'
#' * "close" Return the compositional closure of w into the simplex, 
#'   `w / sum(w)`
#' * "keep" Keep the same total abundance as the original vector by returning
#'   `w * sum(v) / sum(w)`
#' * "none" Return `w` without any normalization
#'
#' @param x An abundance matrix or phyloseq object containing one
#' @param y A numeric vector with which to perturb the observations in x
#' @param margin Matrix margin that corresponds to observations (samples); 
#'   `1` for rows, `2` for columns
#' @param norm String specifying how to normalize the perturbed observations;
#'   see Details.
#' 
#' @rdname perturb
#' @export
perturb <- function(x, y, ...) {
  UseMethod("perturb")
}

#' @rdname perturb
#' @method perturb matrix
#' @export
perturb.matrix <- function(x, y, margin, norm = "close") {
  if (!is.null(names(y)) & !identical(rownames(x), names(y)))
    stop("`rownames(x)` and `names(y)` do not match")
  stopifnot(margin %in% 1:2)
  stopifnot(norm %in% c("close", "keep", "none"))

  if (margin == 1) {
    z <- diag(y) %*% x
    rownames(z) <- rownames(x)
    if (norm %in% c("close", "keep"))
      z <- sweep(z, 2, colSums(z), `/`)
    if (norm == "keep")
      z <- sweep(z, 2, colSums(x), `*`)
  } else {
    z <- x %*% diag(y)
    colnames(z) <- colnames(x)
    if (norm %in% c("close", "keep"))
      z <- sweep(z, 1, rowSums(z), `/`)
    if (norm == "keep")
      z <- sweep(z, 1, rowSums(x), `*`)
  }

  z
}

#' @rdname perturb
#' @method perturb otu_table
#' @export
perturb.otu_table <- function(x, y, norm = "close") {
  if (!is.null(names(y)) & !identical(taxa_names(x), names(y)))
    stop("`taxa_names(x)` and `names(y)` do not match")

  margin <- 2 - taxa_are_rows(x)
  z <- perturb.matrix(as(x, "matrix"), y, margin = margin, norm = norm)
  otu_table(z, taxa_are_rows = taxa_are_rows(x))
}

#' @rdname perturb
#' @method perturb phyloseq
#' @export
perturb.phyloseq <- function(x, y, norm = "close") {
  otu_table(x) <- perturb(otu_table(x), y, norm = norm)
  x
}

setMethod("perturb", c("matrix", "numeric"), perturb.matrix)
setMethod("perturb", c("otu_table", "numeric"), perturb.otu_table)
setMethod("perturb", c("phyloseq", "numeric"), perturb.phyloseq)

# calibrate -------------------------------------------------------------------

#' Calibrate a relative-abundance matrix by a bias vector
#'
#' Calibration via the simple deterministic procedure described in McLaren,
#' Willis, and Callahan (2019), simply involved dividing the observed vector of
#' relative abundances by the estimated bias vector and (optionally)
#' normalizing the result to sum to 1 or some other chosen value.
#'
#' Normalization options specified by `norm`:
#'
#' * "close": Divide the calibrated abundance vector by its sum, so that it
#'   sums to 1
#' * "keep": Keep the same total abundance as the original observation
#' * "none": Return the calibrated abundances without any normalization
#'
#' @param observed An abundance matrix or phyloseq object containing one
#' @param bias A numeric vector of relative efficiencies
#' @param margin Matrix margin that corresponds to observations (samples); 
#'   `1` for rows, `2` for columns
#' @param norm String specifying how to normalize the calibrated observations;
#'   see Details.
#' 
#' @seealso [perturb()]
#' 
#' @rdname calibrate
#' @export
calibrate <- function(observed, bias, ...) {
  UseMethod("calibrate")
}

#' @rdname calibrate
#' @method calibrate matrix
#' @export
calibrate.matrix <- function(observed, bias, margin, norm = "close") {
  perturb.matrix(observed, 1 / bias, margin = margin, norm = norm)
}

#' @rdname calibrate
#' @method calibrate otu_table
#' @export
calibrate.otu_table <- function(observed, bias, norm = "close") {
  perturb.otu_table(observed, 1 / bias, norm = norm)
}

#' @rdname calibrate
#' @method calibrate phyloseq
#' @export
calibrate.phyloseq <- function(observed, bias, norm = "close") {
  otu_table(observed) <- calibrate(otu_table(observed), bias, norm = norm)
  observed
}

setMethod("calibrate", c("matrix", "numeric"), calibrate.matrix)
setMethod("calibrate", c("otu_table", "numeric"), calibrate.otu_table)
setMethod("calibrate", c("phyloseq", "numeric"), calibrate.phyloseq)
