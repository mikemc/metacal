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
#' If `y` is named, then the names must agree with the taxa names in `x` and
#' will be used to reorder `y` to match the taxa order in `x`.
#'
#' @param x An abundance matrix or phyloseq object containing one
#' @param y A numeric vector with which to perturb the observations in x
#' @param margin Matrix margin that corresponds to observations (samples); 
#'   `1` for rows, `2` for columns
#' @param norm String specifying how to normalize the perturbed observations;
#'   see Details.
#'
#' @seealso [calibrate()]
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
  stopifnot(margin %in% 1:2)
  stopifnot(norm %in% c("close", "keep", "none"))
  taxa_margin <- setdiff(1:2, margin)
  stopifnot(identical(dim(x)[[taxa_margin]], length(y)))

  # Align taxa
  if (!is.null(names(y))) {
    if (!setequal(dimnames(x)[[taxa_margin]], names(y)))
      stop("Taxa names in `x` and `y` do not match")
    y <- y[dimnames(x)[[taxa_margin]]]
  }

  if (margin == 2) {
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
  margin <- 1 + taxa_are_rows(x)
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

#' @rdname perturb
#' @method perturb mc_bias_fit
#' @export
perturb.mc_bias_fit <- function(x, y) {
  stopifnot(identical(length(x$estimate), length(y)))
  x$estimate <- x$estimate * y
  if (!is.null(x$bootreps))
    x$bootreps <- sweep(x$bootreps, 2, y, `*`)
  x
}

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
#' If `bias` is named, then the names must agree with the taxa names in
#' `observed` and will be used to reorder `bias` to match the taxa order in
#' `observed`.
#'
#' @param observed An abundance matrix or phyloseq object containing one
#' @param bias A numeric vector of relative efficiencies or an object of class
#'   'mc_bias_fit' (from which efficiencies will be extracted)
#' @param margin Matrix margin that corresponds to observations (samples); 
#'   `1` for rows, `2` for columns
#' @param norm String specifying how to normalize the calibrated observations;
#'   see Details.
#' @param mean_name Character vector or NULL. Name of the column in the sample
#'   data in which to store the mean efficiency, or NULL to skip.
#' 
#' @seealso [perturb()] [estimate_bias()]
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
  if (inherits(bias, 'mc_bias_fit'))
    bias <- coef(bias)
  perturb.matrix(observed, 1 / bias, margin = margin, norm = norm)
}

#' @rdname calibrate
#' @method calibrate otu_table
#' @export
calibrate.otu_table <- function(observed, bias, norm = "close") {
  if (inherits(bias, 'mc_bias_fit'))
    bias <- coef(bias)
  # The otu-table method will try to subset to the shared taxa in observed and
  # bias
  if (!identical(ntaxa(observed), length(bias))) {
    if (is.null(names(bias)))
      stop("`bias` must be named if `!identical(ntaxa(observed), length(bias))`")
    taxa <- intersect(taxa_names(observed), names(bias))
    if (!length(taxa) > 0)
      stop("`intersect(taxa_names(observed), names(bias))` is empty")
    observed <- phyloseq::prune_taxa(taxa, observed)
    bias <- bias[taxa]
  } # else: perturb() will make sure the taxa names match the names in bias
  perturb.otu_table(observed, 1 / bias, norm = norm)
}

#' @rdname calibrate
#' @method calibrate phyloseq
#' @export
calibrate.phyloseq <- function(observed, 
                               bias, 
                               norm = "close", 
                               mean_name = ".mean_efficiency") {
  otu_table(observed) <- calibrate(otu_table(observed), bias, norm = norm)
  if (!is.null(mean_name) & !is.null(phyloseq::access(observed, 'sam_data'))) {
    if (!is.character(mean_name))
      stop('`mean_name` must be a character vector')
    me <- mean_efficiency(
      observed, 
      bias, 
      # `observed` is now calibrated so type = 'actual'
      type = 'actual'
    )
    sample_data(observed)[, mean_name] <- me
  }
  observed
}
