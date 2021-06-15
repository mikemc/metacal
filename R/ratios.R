# Functions for computing ratios

# TODO: ensure that comp_vars default won't pick group_vars; add ability to use
# symbols as arguments; Add ability to specify the taxon_var

#' Compute taxon ratios from a tidy data frame
#'
#' Assumes the taxon variable name is `Taxon`
#'
#' @param drop Whether to drop the individual taxon quantities
#' @export
compute_ratios <- function(.data, 
    # taxon_var = "Taxon",
    group_vars = c("Sample"),
    comp_vars = setdiff(names(dplyr::select_if(.data, is.numeric)), 
        group_vars),
    drop = TRUE) {

    .data <- .data %>%
        dplyr::select(Taxon, group_vars, comp_vars)

    # Get all pairs of taxa and join the abundances.
    tb <- tidyr::expand(.data, tidyr::nesting(!!!rlang::syms(group_vars)), 
        Taxon.x = Taxon, Taxon.y = Taxon) %>%
        # mutate(Pair = paste(Taxon.x, Taxon.y, sep = ":")) %>%
        dplyr::left_join(.data, by = c(group_vars, "Taxon.x" = "Taxon")) %>%
        dplyr::left_join(.data, by = c(group_vars, "Taxon.y" = "Taxon"))
    # Each var in `comp_vars` now appears twice, with a ".x" and a ".y" suffix. 

    # Compute the ratios for each var in `comp_vars`:
    for (v in comp_vars) {
        tb[v] <- tb[paste0(v, ".x")] / tb[paste0(v, ".y")]
    }

    if (drop == TRUE) {
        tb <- tb %>%
            dplyr::select(Taxon.x, Taxon.y, group_vars, comp_vars)
    }

    tb
}

# Pairwise ratio vectors and matrices -----------------------------------------

#' Pairwise ratios of vector elements and matrix rows or columns
#'
#' For a matrix, computes the ratios of the supplied margin (1 = rows, 2 =
#' cols). Otherwise, computes the ratios of the elements.
#' 
#' @param x A vector or matrix
#' @param margin Margin (1 = rows, 2 = cols) of matrix x to compute ratios of
#' @param f Vectorized binary function applied to numerator and denominator
#' @param filter Whether to filter out redundant pairs
#' @param set_names Whether to set names for the pairs
#' @param sep Seperator for pair names
#'
#' @name pairwise_ratios
#' @export
#' @examples
#'   mat <- seq(10) %>%
#'     matrix(nrow = 2)
#'   rownames(mat) <- paste0('r', seq(2))
#'   colnames(mat) <- paste0('c', seq(5))
#'   
#'   pairwise_ratios(mat, 2)
#'   pairwise_ratios(mat[1,])
pairwise_ratios <- function(x, 
                            ..., 
                            f = `/`, 
                            filter = TRUE, 
                            set_names = TRUE, 
                            sep = ":") {
  UseMethod("pairwise_ratios")
}

setGeneric("pairwise_ratios")

#' @rdname pairwise_ratios
pairwise_ratios.default <- function(x,
                                    f = `/`, 
                                    filter = TRUE, 
                                    set_names = TRUE, 
                                    sep = ":") {
  idx <- seq_along(x)
  tb <- tidyr::crossing(i = idx, j = idx)
  if (filter)
    tb <- tb %>% dplyr::filter(i < j)
  ratios <- purrr::pmap_dbl(tb, function(i, j) f(x[i], x[j]))
  if (set_names) {
    nms <- names(x)
    if (is.null(nms)) 
      nms <- idx
    names(ratios) <- purrr::pmap_chr(
      tb, 
      function(i, j) paste(nms[i], nms[j], sep = sep)
    )
  }
  ratios
}
# NOTE: Could also be done with `outer(x, x, `/`) %>% as.vector`, but this way
# is ~20x faster and easily extends to work with matrices

#' @rdname pairwise_ratios
pairwise_ratios.matrix <- function(x, 
                                   margin, 
                                   f = `/`, 
                                   filter = TRUE, 
                                   set_names = TRUE, 
                                   sep = ":") {
  # How are the names of the other margin normally carried over?
  # Currently they are lost when there is only one element in that dimension
  idx <- seq(dim(x)[[margin]])
  tb <- tidyr::crossing(i = idx, j = idx)
  if (filter)
    tb <- tb %>% dplyr::filter(i < j)
  if (margin == 1) {
    ratios <- purrr::pmap(tb, 
      function(i, j) f(x[i, , drop = FALSE], x[j, , drop = FALSE])
    )
    ratio_mat <- do.call(rbind, ratios)
  } else if (margin == 2) {
    ratios <- purrr::pmap(tb, 
      function(i, j) f(x[, i, drop = FALSE], x[, j, drop = FALSE])
    )
    ratio_mat <- do.call(cbind, ratios)
  } else
    stop("margin must be 1 or 2")
  if (set_names) {
    nms <- dimnames(x)[[margin]]
    if (is.null(nms)) 
      nms <- idx
    dimnames(ratio_mat)[[margin]] <- purrr::pmap_chr(
      tb, 
      function(i, j) paste(nms[i], nms[j], sep = sep)
    )
  }
  ratio_mat
}

# Version that works with phyloseq otu_table objects
# TODO: check that phyloseq is installed and throw and error if not

#' @rdname pairwise_ratios
pairwise_ratios.otu_table <- function(x, 
                                      margin = "taxa", 
                                      f = `/`, 
                                      filter = TRUE, 
                                      sep = ":") {
  if (!requireNamespace("phyloseq", quietly = TRUE))
    stop("Phyloseq package required for this function")
  stopifnot(margin %in% c("taxa", "samples"))
  taxa_are_rows <- phyloseq::taxa_are_rows(x)
  xmat <- as(x, "matrix")
  # Ensure that taxa/features are columns
  if (taxa_are_rows)
    xmat <- t(xmat)
  # Compute ratios of taxa (default) or of samples
  if (margin == "taxa") {
    mat_margin <- 2
  } else if (margin == "samples") {
    mat_margin <- 1
  }
  ratios <- pairwise_ratios.matrix(
    xmat,
    mat_margin,
    f = f,
    filter = filter,
    set_names = TRUE,
    sep = sep
  )
  # Return as an otu_table in original orientation
  if (taxa_are_rows)
    ratios <- t(ratios)
  phyloseq::otu_table(ratios, taxa_are_rows)
}

setMethod("pairwise_ratios", "otu_table", pairwise_ratios.otu_table)

#' @importClassesFrom phyloseq phyloseq
#' @rdname pairwise_ratios
pairwise_ratios.phyloseq <- function(x, 
                                     margin = "taxa", 
                                     f = `/`, 
                                     filter = TRUE, 
                                     sep = ":") {
  if (!requireNamespace("phyloseq", quietly = TRUE))
    stop("Phyloseq package required for this function")
  # New otu table
  new_otu <- pairwise_ratios.otu_table(
    phyloseq::otu_table(x), 
    margin = margin,
    f = f,
    filter = filter,
    sep = sep
  )
  if (margin == "taxa") {
    return(phyloseq::phyloseq(new_otu, phyloseq::sample_data(x)))
  }
  # New sample data needed when margin = "samples". In this case, should have
  # columns "sample.1" and "sample.2" and the data for both.
  if (margin == "samples") {
    sam <- phyloseq::sample_data(x)
    samtb <- sam %>%
      as("data.frame") %>%
      tibble::as_tibble(rownames = "sample")
    tb <- tidyr::crossing(
      sample.1 = phyloseq::sample_names(x),
      sample.2 = phyloseq::sample_names(x)
      )
    tb <- tb %>%
      dplyr::left_join(samtb, by = c(sample.1 = "sample")) %>%
      dplyr::left_join(samtb, by = c(sample.2 = "sample"), 
        suffix = c(".1", ".2")) %>%
      dplyr::mutate(pair = paste(sample.1, sample.2, sep = sep))
    new_sam <- tb %>%
      tibble::column_to_rownames("pair") %>%
      phyloseq::sample_data()
    # Call to phyloseq() will filter sample pairs that are missing in new_otu
    # if `filter == TRUE`
    return(phyloseq::phyloseq(new_otu, new_sam))
  }
}
# TODO: check that this method of filtering samples is ok; and see which order
# they end up in. Might instead do to exactly match the otu_table order
#

setMethod("pairwise_ratios", "phyloseq", pairwise_ratios.phyloseq)

# TODO: preserve taxonomy and refseq information?

# idx <- seq(dim(x)[[margin]])
# tb <- tidyr::crossing(i = idx, j = idx)
# if (filter)
#   tb <- tb %>% dplyr::filter(i < j)
